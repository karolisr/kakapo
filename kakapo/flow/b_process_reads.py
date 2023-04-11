"""Kakapo workflow: Process Reads."""

import fileinput
import pickle
# import random
import re

from collections import OrderedDict
from copy import deepcopy
from os import remove
from os import stat
from os.path import basename
from os.path import commonprefix
from os.path import exists as ope
from os.path import isfile
from os.path import join as opj
from os.path import sep as ops
from os.path import splitext
from shutil import copyfile
from shutil import move
from shutil import rmtree
from time import sleep

from kakapo.tools.bioio import read_fasta
from kakapo.tools.bioio import seq_records_to_dict
from kakapo.tools.bioio import write_fasta
from kakapo.tools.blast import make_blast_db
from kakapo.tools.bowtie2 import build_bt2_index, run_bowtie2_se, run_bowtie2_pe
from kakapo.tools.config import PICKLE_PROTOCOL
from kakapo.tools.eutils import accs as accs_eutil
from kakapo.tools.eutils import search as search_eutil
from kakapo.tools.eutils import seqs as dnld_ncbi_seqs
from kakapo.tools.eutils import sra_run_info
from kakapo.tools.kraken import run_kraken_filters
from kakapo.tools.rcorrector import filter_unc_se, filter_unc_pe
from kakapo.tools.rcorrector import run_rcorrector_se, run_rcorrector_pe
from kakapo.tools.seq import SEQ_TYPE_NT
from kakapo.tools.seqtk import seqtk_fq_to_fa
from kakapo.tools.transl_tables import TranslationTable
from kakapo.tools.trimmomatic import trimmomatic_se, trimmomatic_pe
from kakapo.utils.logging import Log
from kakapo.utils.misc import make_dirs
from kakapo.utils.misc import plain_or_gzip
from kakapo.utils.misc import rename_fq_seqs
from kakapo.utils.misc import splitext_gz
from kakapo.utils.misc import avg_read_len_fq
from kakapo.utils.subp import run


MT = 'mitochondrion'
PT = 'plastid'

MAX_READ_LEN_ILLUMINA = 1023
MIN_READ_COUNT_SRA = 10000


def dnld_sra_info(sras, dir_cache_prj):

    sra_runs_info = {}
    sras_acceptable = []

    if len(sras) > 0:
        Log.inf('Downloading SRA run information.')
    else:
        return sra_runs_info, sras_acceptable

    __ = opj(dir_cache_prj, 'sra_runs_info_cache')

    if ope(__):
        with open(__, 'rb') as f:
            sra_runs_info = pickle.load(f)

    sras_local = [k for k in sra_runs_info.keys()]
    sras_to_dnld = set(sras).difference(set(sras_local))
    if len(sras_to_dnld) > 0:
        temp = sra_run_info(list(sras_to_dnld))
        new_sra_runs_info = {i['Run']: i for i in temp}
        sra_runs_info.update(new_sra_runs_info)

    for sra in sras:

        if sra in sra_runs_info:

            info = sra_runs_info[sra]

            sra_lib_layout = info['LibraryLayout'].lower()
            sra_lib_source = info['LibrarySource'].lower()
            sra_lib_strategy = info['LibraryStrategy']
            sra_seq_platform = info['Platform'].lower().capitalize()
            sra_seq_platform_model = info['Model']

            sra_species = info['ScientificName'].replace('/', '_')
            sra_runs_info[sra]['ScientificName'] = sra_species

            sra_taxid = info['TaxID']
            sra_spots = int(info['spots'])
            sra_spots_with_mates = int(info['spots_with_mates'])

            sample_base_name = (sra_species.replace(' ', '_') + '_' +
                                sra_taxid + '_' + sra)

            sra_runs_info[sra]['KakapoSampleBaseName'] = sample_base_name

            src_check = sra_lib_source.lower()
            strategy_check = sra_lib_strategy.lower()

            if not ('transcript' in src_check or
                    'rna' in src_check or
                    'rna' in strategy_check):

                sra_info_str = ('{sra}: the SRA library source type "{ltype}" '
                                'or library strategy "{strategy}" '
                                'is not supported.').format(
                                    sra=sra, ltype=sra_lib_source,
                                    strategy=sra_lib_strategy)

                Log.err(sra_info_str, 'Skipping.')

            elif sra_seq_platform != 'Illumina':
                sra_info_str = ('{sra}: the SRA library sequencing platform '
                                '"{plat}" is not supported').format(
                                    sra=sra, plat=sra_seq_platform)

                Log.err(sra_info_str, 'Skipping.')

            else:

                Log.msg('{sra}:'.format(sra=sra),
                        '{strategy} {layout}-end library ({source}).'.format(
                            strategy=sra_lib_strategy,
                            layout=sra_lib_layout,
                            source=sra_lib_source.title()))

                sra_runs_info[sra]['KakapoLibraryLayout'] = \
                    sra_runs_info[sra]['LibraryLayout']

                if sra_lib_layout == 'paired' and sra_spots_with_mates == 0:
                    sra_runs_info[sra]['KakapoLibraryLayout'] = 'SINGLE'
                    Log.wrn('      Note:', 'Listed as containing '
                            'paired-end reads, but only a single set of reads '
                            'is available. Treating as single-ended.')

                elif (sra_lib_layout == 'paired' and
                      sra_spots != sra_spots_with_mates):
                    sra_runs_info[sra]['KakapoLibraryLayout'] = 'PAIRED_UNP'
                    Log.wrn('      Note:', 'Listed as containing '
                            'paired-end reads, but not all reads are paired.')

                Log.msg('    Source:',
                        '{species} (TaxID: {txid}).'.format(
                            species=sra_species,
                            txid=sra_taxid), False)
                Log.msg('Technology:',
                        '{platform} platform / {model}.'.format(
                            platform=sra_seq_platform,
                            model=sra_seq_platform_model), False)

                sras_acceptable.append(sra)

    with open(__, 'wb') as f:
        pickle.dump(sra_runs_info, f, protocol=PICKLE_PROTOCOL)

    return sra_runs_info, sras_acceptable


def dnld_sra_fastq_files(sras, sra_runs_info, dir_fq_data, fasterq_dump,
                         threads, dir_temp):

    if len(sras) > 0:
        if fasterq_dump is None:
            Log.err('fasterq-dump from SRA Toolkit is not available. ' +
                    'Cannot continue. Exiting.')
            exit(0)

        print()
        Log.inf('Downloading SRA read data.')

    se_fastq_files = {}
    pe_fastq_files = {}

    for sra in sras:
        run_info = sra_runs_info[sra]
        sra_lib_layout = run_info['LibraryLayout'].lower()
        sra_lib_layout_k = run_info['KakapoLibraryLayout'].lower()
        sample_base_name = run_info['KakapoSampleBaseName']
        sra_taxid = int(run_info['TaxID'])
        avg_len = int(run_info['avgLength'])

        sra_dnld_needed = False

        se_file = None

        pe_file_1 = None
        pe_file_2 = None
        pe_file_3 = None

        pe_file_1_renamed = None
        pe_file_2_renamed = None
        pe_file_3_renamed = None

        if sra_lib_layout == 'single' or sra_lib_layout_k == 'single':
            se_file = opj(dir_fq_data, sra + '.fastq')
            se_file_gz = se_file + '.gz'

            se_file_blacklisted_gz = opj(dir_fq_data, sra + '_blacklisted.fastq.gz')
            if ope(se_file_blacklisted_gz):
                Log.wrn('Ignoring blacklisted SRA:', sra)
                continue

            se_fastq_files[sample_base_name] = {'path': se_file_gz}
            se_fastq_files[sample_base_name]['src'] = 'sra'
            se_fastq_files[sample_base_name]['avg_len'] = avg_len
            se_fastq_files[sample_base_name]['tax_id'] = sra_taxid
            se_fastq_files[sample_base_name]['organelle'] = None
            if not ope(se_file_gz):
                sra_dnld_needed = True

        elif sra_lib_layout == 'paired':
            pe_file_1 = opj(dir_fq_data, sra + '_1.fastq')
            pe_file_2 = opj(dir_fq_data, sra + '_2.fastq')
            pe_file_1_renamed = opj(dir_fq_data, sra + '_R1.fastq')
            pe_file_2_renamed = opj(dir_fq_data, sra + '_R2.fastq')

            pe_file_1_renamed_gz = pe_file_1_renamed + '.gz'
            pe_file_2_renamed_gz = pe_file_2_renamed + '.gz'

            pe_file_1_renamed_blacklisted_gz = opj(
                dir_fq_data, sra + '_R1_blacklisted.fastq.gz')

            if ope(pe_file_1_renamed_blacklisted_gz):
                Log.wrn('Ignoring blacklisted SRA:', sra)
                continue

            pe_fastq_files[sample_base_name] = {'path': [pe_file_1_renamed_gz,
                                                         pe_file_2_renamed_gz]}
            pe_fastq_files[sample_base_name]['src'] = 'sra'
            pe_fastq_files[sample_base_name]['avg_len'] = avg_len // 2
            pe_fastq_files[sample_base_name]['tax_id'] = sra_taxid
            pe_fastq_files[sample_base_name]['organelle'] = None
            if sra_lib_layout_k == 'paired_unp':
                pe_file_3 = opj(dir_fq_data, sra + '.fastq')
                pe_file_3_renamed = opj(dir_fq_data, sra + '_R3.fastq')
                pe_file_3_renamed_gz = pe_file_3_renamed + '.gz'
                pe_fastq_files[sample_base_name]['path'].append(
                    pe_file_3_renamed_gz)
            if not ope(pe_file_1_renamed_gz) or not ope(pe_file_2_renamed_gz):
                sra_dnld_needed = True

        if not sra_dnld_needed:
            Log.msg('FASTQ reads are available locally:', sample_base_name)

        retry_count = 0
        while sra_dnld_needed:

            if retry_count > 50:
                Log.err('Download failed. Exiting.')
                rmtree(dir_temp)
                exit(1)

            elif retry_count > 0:
                Log.wrn('Download failed. Retrying.')
                sleep(2)

            retry_count += 1

            Log.msg('Downloading FASTQ reads for:', sample_base_name)

            _threads = min(6, threads)
            _bufsize = 100
            # ToDo: Check the line below if _mem is NOT above RAM available.
            _mem = _threads * 2 * 100

            cmd = [fasterq_dump,
                   '--threads', str(_threads),
                   '--split-3',
                   '--bufsize', str(_bufsize) + 'M',
                   '--mem', str(_mem) + 'M',
                   '--outdir', dir_fq_data,
                   '--temp', dir_temp, sra]

            run(cmd, do_not_raise=True)

            if sra_lib_layout == 'single' or sra_lib_layout_k == 'single':
                if not ope(se_file):
                    continue

            elif sra_lib_layout == 'paired':
                if not ope(pe_file_1) or not ope(pe_file_2):
                    continue
                else:
                    move(pe_file_1, pe_file_1_renamed)
                    move(pe_file_2, pe_file_2_renamed)

                if sra_lib_layout_k == 'paired_unp':
                    if not ope(pe_file_3):
                        continue
                    else:
                        move(pe_file_3, pe_file_3_renamed)

            sra_dnld_needed = False

            if sra_lib_layout == 'single' or sra_lib_layout_k == 'single':
                if ope(se_file):
                    Log.msg('Renaming FASTQ reads in:', se_file)
                    read_count = rename_fq_seqs(se_file, sra, '1:N:0')

                    if read_count < MIN_READ_COUNT_SRA:
                        Log.wrn('FASTQ file contains less than ' +
                                str(MIN_READ_COUNT_SRA) +
                                ' reads, blacklisting:', sra)
                        del se_fastq_files[sample_base_name]
                        move(se_file_gz, se_file_gz.replace(
                            '.fastq',
                            '_blacklisted.fastq'))

            elif sra_lib_layout == 'paired':

                read_1_count = 0
                read_3_count = 0

                if ope(pe_file_1_renamed):
                    Log.msg('Renaming FASTQ reads in:', pe_file_1_renamed)
                    read_1_count = rename_fq_seqs(pe_file_1_renamed, sra,
                                                  '1:N:0')

                if ope(pe_file_2_renamed):
                    Log.msg('Renaming FASTQ reads in:', pe_file_2_renamed)
                    read_2_count = rename_fq_seqs(pe_file_2_renamed, sra,
                                                  '2:N:0')
                    assert read_1_count == read_2_count

                if sra_lib_layout_k == 'paired_unp':
                    if ope(pe_file_3_renamed):
                        Log.msg('Renaming FASTQ reads in:', pe_file_3_renamed)
                        read_3_count = rename_fq_seqs(pe_file_3_renamed, sra +
                                                      '_unpaired', '1:N:0')

                if (read_1_count + read_3_count) < MIN_READ_COUNT_SRA:
                    Log.wrn('FASTQ files contain less than ' +
                            str(MIN_READ_COUNT_SRA) + ' reads, blacklisting:',
                            sra)

                    del pe_fastq_files[sample_base_name]
                    move(pe_file_1_renamed_gz,
                         pe_file_1_renamed_gz.replace('.fastq',
                                                      '_blacklisted.fastq'))
                    move(pe_file_2_renamed_gz,
                         pe_file_2_renamed_gz.replace('.fastq',
                                                      '_blacklisted.fastq'))
                    if ope(pe_file_3_renamed_gz):
                        move(pe_file_3_renamed_gz,
                             pe_file_3_renamed_gz.replace('.fastq',
                                                          '_blacklisted.fastq'))

    return se_fastq_files, pe_fastq_files


def user_fastq_files(fq_se, fq_pe):
    if len(fq_se) > 0 or len(fq_pe) > 0:
        Log.inf('Preparing user provided FASTQ files.')

    se_fastq_files = {}
    pe_fastq_files = {}

    fq_type_1_regex = r'(.*)_L\d\d\d(_R.)_\d\d\d(.*)'

    for se in fq_se:
        tax_id = se[0]
        path = se[1]
        base = basename(path)
        if plain_or_gzip(base)[4] != '':
            base = splitext(base)[0]
        base = splitext(base)[0]
        fq_type_1_match = re.findall(fq_type_1_regex, base)
        if len(fq_type_1_match) > 0 and len(fq_type_1_match[0]) == 3:
            base = fq_type_1_match[0][0]
        sample_base_name = base.rstrip('_').rstrip('-')
        se_fastq_files[sample_base_name] = {'path': path}
        se_fastq_files[sample_base_name]['src'] = 'usr'
        se_fastq_files[sample_base_name]['avg_len'] = None
        se_fastq_files[sample_base_name]['tax_id'] = tax_id
        se_fastq_files[sample_base_name]['organelle'] = None
        Log.msg(sample_base_name + ':', basename(path))

    for pe in fq_pe:
        tax_id = pe[0]
        path = pe[1]
        base = basename(path[0])
        if plain_or_gzip(base)[4] != '':
            base = splitext(base)[0]
        base = splitext(base)[0]
        fq_type_1_match = re.findall(fq_type_1_regex, base)
        if len(fq_type_1_match) > 0 and len(fq_type_1_match[0]) == 3:
            base = fq_type_1_match[0][0]
        else:
            base = basename(commonprefix(path)).rstrip('_- R')
        sample_base_name = base.rstrip('_').rstrip('-')
        pe_fastq_files[sample_base_name] = {'path': path}
        pe_fastq_files[sample_base_name]['src'] = 'usr'
        pe_fastq_files[sample_base_name]['avg_len'] = None
        pe_fastq_files[sample_base_name]['tax_id'] = tax_id
        pe_fastq_files[sample_base_name]['organelle'] = None
        Log.msg(sample_base_name + ':', basename(path[0]) + '\n' +
                ' ' * (len(sample_base_name) + 2) + basename(path[1]))

    return se_fastq_files, pe_fastq_files


def min_accept_read_len(se_fastq_files, pe_fastq_files, dir_temp,
                        dir_cache_fq_minlen):
    # lowest allowable
    low = 25

    if len(se_fastq_files) > 0 or len(pe_fastq_files) > 0:
        print()
        Log.inf('Calculating minimum acceptable read length.')
    else:
        return None

    __ = opj(dir_cache_fq_minlen, 'minlen')

    pickled = {}

    if ope(__):
        with open(__, 'rb') as f:
            pickled = pickle.load(f)

    queue = []

    for se in se_fastq_files:
        # src = se_fastq_files[se]['src']
        # avg_len = se_fastq_files[se]['avg_len']
        # if src == 'sra':
        #     ml = max(avg_len // 3, low)
        #     se_fastq_files[se]['min_acc_len'] = ml
        #     Log.msg(str(ml) + ' nt:', se)
        #     continue

        # fq_path = se_fastq_files[se]['path']
        stats_file = opj(dir_temp, se + '_stats.txt')
        queue.append([se, [se_fastq_files[se]['path']], stats_file, 'se'])

    for pe in pe_fastq_files:
        # src = pe_fastq_files[pe]['src']
        # avg_len = pe_fastq_files[pe]['avg_len']
        # if src == 'sra':
        #     ml = max(avg_len // 3, low)
        #     pe_fastq_files[pe]['min_acc_len'] = ml
        #     Log.msg(str(ml) + ' nt:', pe)
        #     continue

        # fq_path = pe_fastq_files[pe]['path'][0]
        stats_file = opj(dir_temp, pe + '_stats.txt')
        queue.append([pe, pe_fastq_files[pe]['path'], stats_file, 'pe'])

    for x in queue:

        if x[0] in pickled:
            ml = pickled[x[0]]
        else:
            ml = avg_read_len_fq(x[1][0])

            if ml < MAX_READ_LEN_ILLUMINA:
                ml = max(int(ml) // 5, low)
            else:
                ml = MAX_READ_LEN_ILLUMINA

            pickled[x[0]] = ml

        if ml is not None:
            if ml < MAX_READ_LEN_ILLUMINA:
                Log.msg(str(ml) + ' nt:', x[0])
            else:
                Log.wrn('FASTQ file(s) contain reads longer than ' +
                        str(MAX_READ_LEN_ILLUMINA) + ' bp, blacklisting:',
                        x[0])
                for f_to_blacklist in x[1]:
                    if ope(f_to_blacklist):
                        move(f_to_blacklist, f_to_blacklist.replace(
                            '.fastq', '_blacklisted.fastq'))
        else:
            Log.msg(' ?' + ' nt:', x[0])
            ml = low

        if x[3] == 'se':

            if ml >= MAX_READ_LEN_ILLUMINA:
                del se_fastq_files[x[0]]
                continue

            se_fastq_files[x[0]]['min_acc_len'] = ml

        elif x[3] == 'pe':

            if ml >= MAX_READ_LEN_ILLUMINA:
                del pe_fastq_files[x[0]]
                continue

            pe_fastq_files[x[0]]['min_acc_len'] = ml

        with open(__, 'wb') as f:
            pickle.dump(pickled, f, protocol=PICKLE_PROTOCOL)


def file_name_patterns():
    pe_trim_pair_1_sfx = '_paired_1'
    pe_trim_pair_2_sfx = '_paired_2'
    pe_trim_unpr_1_sfx = '_unpaired_1'
    pe_trim_unpr_2_sfx = '_unpaired_2'

    pe_trim_suffixes = [pe_trim_pair_1_sfx, pe_trim_pair_2_sfx,
                        pe_trim_unpr_1_sfx, pe_trim_unpr_2_sfx]

    pe_file_pattern = opj('@D@', '@N@')

    pe_trim_fq_file_patterns = list(
        zip([pe_file_pattern] * 4, pe_trim_suffixes, ['.fastq'] * 4))
    pe_trim_fq_file_patterns = [''.join(x) for x in pe_trim_fq_file_patterns]

    pe_trim_fa_file_patterns = [x.replace('.fastq', '.fasta') for x in
                                pe_trim_fq_file_patterns]

    pe_blast_db_file_patterns = list(
        zip([pe_file_pattern] * 4, pe_trim_suffixes))
    pe_blast_db_file_patterns = [''.join(x) for x in pe_blast_db_file_patterns]

    pe_blast_results_file_patterns = [x.replace('.fastq', '__@Q@.txt') for x in
                                      pe_trim_fq_file_patterns]

    pe_vsearch_results_file_patterns = pe_blast_results_file_patterns

    return (pe_trim_fq_file_patterns, pe_trim_fa_file_patterns,
            pe_blast_db_file_patterns, pe_blast_results_file_patterns,
            pe_vsearch_results_file_patterns)


def run_trimmomatic(se_fastq_files, pe_fastq_files, dir_fq_trim_data,
                    trimmomatic, adapters, fpatt, threads,
                    rcorrector_before_trimmomatic, gz_out=True):
    if len(se_fastq_files) > 0 or len(pe_fastq_files) > 0:
        print()
        Log.inf('Running Trimmomatic.')
        if trimmomatic is None:
            Log.err('trimmomatic is not available. Cannot continue. Exiting.')
            exit(0)
    for se in se_fastq_files:
        dir_fq_trim_data_sample = opj(dir_fq_trim_data, se)

        if rcorrector_before_trimmomatic is True:
            fq_path = se_fastq_files[se]['cor_path_fq']
        else:
            fq_path = se_fastq_files[se]['path']

        _, _, _, _, ext_in = plain_or_gzip(fq_path)

        ext_out = ext_in
        if gz_out is True:
            ext_out = '.gz'

        min_acc_len = se_fastq_files[se]['min_acc_len']
        stats_f = opj(dir_fq_trim_data_sample, se + '.txt')
        out_f = opj(dir_fq_trim_data_sample, se + '.fastq' + ext_out)
        se_fastq_files[se]['trim_path_fq'] = out_f

        if ope(dir_fq_trim_data_sample):
            Log.msg('Trimmed FASTQ file already exists:', se)
        else:
            make_dirs(dir_fq_trim_data_sample)
            Log.msg('SE mode:', se)
            trimmomatic_se(
                trimmomatic=trimmomatic,
                adapters=adapters,
                in_file=fq_path,
                out_file=out_f,
                stats_file=stats_f,
                threads=threads,
                minlen=min_acc_len)

    for pe in pe_fastq_files:
        dir_fq_trim_data_sample = opj(dir_fq_trim_data, pe)

        if rcorrector_before_trimmomatic is True:
            pe_key = 'cor_path_fq'
        else:
            pe_key = 'path'

        fq_path_1 = pe_fastq_files[pe][pe_key][0]
        fq_path_2 = pe_fastq_files[pe][pe_key][1]
        fq_path_3 = None
        if len(pe_fastq_files[pe][pe_key]) == 3:
            fq_path_3 = pe_fastq_files[pe][pe_key][2]

        _, _, _, _, ext_in = plain_or_gzip(fq_path_1)

        ext_out = ext_in
        if gz_out is True:
            ext_out = '.gz'

        min_acc_len = pe_fastq_files[pe]['min_acc_len']
        stats_f = opj(dir_fq_trim_data_sample, pe + '.txt')
        out_fs = [x.replace('@D@', dir_fq_trim_data_sample) for x in fpatt]
        out_fs = [x.replace('@N@', pe) for x in out_fs]
        out_fs = [x + ext_out for x in out_fs]
        pe_fastq_files[pe]['trim_path_fq'] = out_fs

        if ope(dir_fq_trim_data_sample):
            Log.msg('Trimmed FASTQ files already exist:', pe)
        else:
            make_dirs(dir_fq_trim_data_sample)
            Log.msg('PE mode:', pe)
            trimmomatic_pe(
                trimmomatic=trimmomatic,
                adapters=adapters,
                in_file_1=fq_path_1,
                in_file_2=fq_path_2,
                out_file_paired_1=out_fs[0],
                out_file_paired_2=out_fs[1],
                out_file_unpaired_1=out_fs[2],
                out_file_unpaired_2=out_fs[3],
                stats_file=stats_f,
                threads=threads,
                minlen=min_acc_len)

            if fq_path_3 is not None:

                out_f = opj(dir_fq_trim_data_sample, 'unpaired.fastq' + ext_out)
                stats_f = opj(dir_fq_trim_data_sample, pe + '_unpaired.txt')

                Log.msg(
                    'SE mode (Paired-read SRA run contains unpaired reads):',
                    pe)

                trimmomatic_se(
                    trimmomatic=trimmomatic,
                    adapters=adapters,
                    in_file=fq_path_3,
                    out_file=out_f,
                    stats_file=stats_f,
                    threads=threads,
                    minlen=min_acc_len)

                p_temp = opj(dir_fq_trim_data_sample, 'temp.fastq' + ext_out)
                _, w_mode, _, fqopen, _ = plain_or_gzip(p_temp)
                f_temp = fqopen(p_temp, w_mode)
                with fileinput.FileInput(
                        files=[out_fs[2], out_f],
                        openhook=fileinput.hook_compressed) as f:
                    for line in f:
                        f_temp.write(line)
                f_temp.close()

                remove(out_fs[2])
                remove(out_f)
                copyfile(p_temp, out_fs[2])
                remove(p_temp)


def run_rcorrector(se_fastq_files, pe_fastq_files, dir_fq_cor_data, rcorrector,
                   threads, dir_temp, fpatt, should_run,
                   rcorrector_before_trimmomatic, gz_out=True):
    if len(se_fastq_files) > 0 or len(pe_fastq_files) > 0:
        print()
        if should_run is False:
            Log.wrn('Skipping Rcorrector as requested.')
        else:
            Log.inf('Running Rcorrector.')

            if rcorrector is None:
                Log.err('Rcorrector is not available. Cannot continue. '
                        'Exiting.')
                exit(0)

    for se in se_fastq_files:
        dir_fq_cor_data_sample = opj(dir_fq_cor_data, se)

        if rcorrector_before_trimmomatic is True:
            fq_path = se_fastq_files[se]['path']
        else:
            fq_path = se_fastq_files[se]['trim_path_fq']

        _, _, _, _, ext_in = plain_or_gzip(fq_path)

        ext_out = ext_in
        if gz_out is True:
            ext_out = '.gz'

        log_f = opj(dir_fq_cor_data_sample, se + '.txt')
        out_f = opj(dir_fq_cor_data_sample, se + '.fastq' + ext_out)

        se_fastq_files[se]['cor_path_fq'] = out_f

        if should_run is False:
            se_fastq_files[se]['cor_path_fq'] = fq_path
            continue

        if ope(dir_fq_cor_data_sample):
            Log.msg('Corrected FASTQ file already exists:', se)
        else:
            make_dirs(dir_fq_cor_data_sample)
            Log.msg('SE mode:', se)
            run_rcorrector_se(rcorrector=rcorrector,
                              in_file=fq_path,
                              out_dir=dir_fq_cor_data_sample,
                              threads=threads,
                              dir_temp=dir_temp)

            fq_base_path = opj(dir_fq_cor_data_sample, basename(fq_path))
            fq_cor_path = splitext_gz(fq_base_path)[0] + '.cor.fq' + ext_in

            filter_unc_se(in_file=fq_cor_path, out_file=out_f, log_file=log_f)

            remove(fq_cor_path)

    for pe in pe_fastq_files:
        dir_fq_cor_data_sample = opj(dir_fq_cor_data, pe)

        if rcorrector_before_trimmomatic is True:
            fq_path_1 = pe_fastq_files[pe]['path'][0]
            fq_path_2 = pe_fastq_files[pe]['path'][1]
            fq_path_3 = None
            if len(pe_fastq_files[pe]['path']) == 3:
                fq_path_3 = pe_fastq_files[pe]['path'][2]
            fq_path_4 = None
        else:
            fq_path_1 = pe_fastq_files[pe]['trim_path_fq'][0]
            fq_path_2 = pe_fastq_files[pe]['trim_path_fq'][1]
            fq_path_3 = pe_fastq_files[pe]['trim_path_fq'][2]
            fq_path_4 = pe_fastq_files[pe]['trim_path_fq'][3]

        _, _, _, _, ext_in = plain_or_gzip(fq_path_1)

        ext_out = ext_in
        if gz_out is True:
            ext_out = '.gz'

        log_f = opj(dir_fq_cor_data_sample, pe + '_paired.txt')

        out_fs = [x.replace('@D@', dir_fq_cor_data_sample) for x in fpatt]
        out_fs = [x.replace('@N@', pe) for x in out_fs]
        out_fs = [x + ext_out for x in out_fs]

        if fq_path_3 is None:
            out_fs = [out_fs[0], out_fs[1]]

        elif fq_path_4 is None:
            out_fs[2] = out_fs[2].replace('_unpaired_1', '_unpaired')
            out_fs = [out_fs[0], out_fs[1], out_fs[2]]

        pe_fastq_files[pe]['cor_path_fq'] = out_fs

        if should_run is False:
            if fq_path_3 is None:
                pe_fastq_files[pe]['cor_path_fq'] = [fq_path_1, fq_path_2]
            elif fq_path_4 is None:
                pe_fastq_files[pe]['cor_path_fq'] = [fq_path_1, fq_path_2,
                                                     fq_path_3]
            else:
                pe_fastq_files[pe]['cor_path_fq'] = [fq_path_1, fq_path_2,
                                                     fq_path_3, fq_path_4]
            continue

        if ope(dir_fq_cor_data_sample):
            Log.msg('Corrected FASTQ files already exist:', pe)
        else:
            make_dirs(dir_fq_cor_data_sample)
            Log.msg('PE mode:', pe)
            run_rcorrector_pe(rcorrector=rcorrector,
                              in_file_1=fq_path_1,
                              in_file_2=fq_path_2,
                              out_dir=dir_fq_cor_data_sample,
                              threads=threads,
                              dir_temp=dir_temp)

            fq_base_path_1 = opj(dir_fq_cor_data_sample, basename(fq_path_1))
            fq_cor_path_1 = splitext_gz(fq_base_path_1)[0] + '.cor.fq' + ext_in
            fq_base_path_2 = opj(dir_fq_cor_data_sample, basename(fq_path_2))
            fq_cor_path_2 = splitext_gz(fq_base_path_2)[0] + '.cor.fq' + ext_in

            filter_unc_pe(in_file_1=fq_cor_path_1,
                          in_file_2=fq_cor_path_2,
                          out_file_1=out_fs[0],
                          out_file_2=out_fs[1],
                          log_file=log_f)

            remove(fq_cor_path_1)
            remove(fq_cor_path_2)

            # unpaired 1
            if fq_path_3 is not None and ope(fq_path_3) and \
                    stat(fq_path_3).st_size > 512:

                run_rcorrector_se(rcorrector=rcorrector,
                                  in_file=fq_path_3,
                                  out_dir=dir_fq_cor_data_sample,
                                  threads=threads,
                                  dir_temp=dir_temp)

                fq_base_path_3 = opj(dir_fq_cor_data_sample,
                                     basename(fq_path_3))
                fq_cor_path_3 = splitext_gz(fq_base_path_3)[0] + \
                    '.cor.fq' + ext_in

                if fq_path_4 is None:
                    log_f_3 = opj(dir_fq_cor_data_sample, pe +
                                  '_unpaired.txt')
                else:
                    log_f_3 = opj(dir_fq_cor_data_sample, pe +
                                  '_unpaired_1.txt')

                filter_unc_se(in_file=fq_cor_path_3, out_file=out_fs[2],
                              log_file=log_f_3)

                remove(fq_cor_path_3)
            else:
                if rcorrector_before_trimmomatic is False:
                    with open(out_fs[2], 'w') as f:
                        f.write('')

            # unpaired 2
            if fq_path_4 is not None and ope(fq_path_4) and \
                    stat(fq_path_4).st_size > 512:

                run_rcorrector_se(rcorrector=rcorrector,
                                  in_file=fq_path_4,
                                  out_dir=dir_fq_cor_data_sample,
                                  threads=threads,
                                  dir_temp=dir_temp)

                fq_base_path_4 = opj(dir_fq_cor_data_sample,
                                     basename(fq_path_4))
                fq_cor_path_4 = splitext_gz(fq_base_path_4)[0] + \
                    '.cor.fq' + ext_in
                log_f_4 = opj(dir_fq_cor_data_sample,
                              pe + '_unpaired_2.txt')

                filter_unc_se(in_file=fq_cor_path_4, out_file=out_fs[3],
                              log_file=log_f_4)
                remove(fq_cor_path_4)

            else:
                if rcorrector_before_trimmomatic is False:
                    with open(out_fs[3], 'w') as f:
                        f.write('')


def dnld_refseqs_for_taxid(taxid, filter_term, taxonomy, dir_cache_refseqs,
                           query='', db='nuccore'):

    if filter_term == 'plastid':
        ft = '("chloroplast"[filter] OR "plastid"[filter])'
    else:
        ft = '("' + filter_term + '"[filter])'

    tax_terms = tuple(reversed(taxonomy.lineage_for_taxid(taxid)['names']))
    tax_term = ''
    accs = set()
    for tax_term in tax_terms:
        if tax_term is None:
            tax_term = taxonomy.scientific_name_for_taxid(taxid)
        term = '"RefSeq"[Keyword] AND "{}"[Primary Organism] AND {}'
        term = query + term.format(tax_term, ft)
        accs = set(accs_eutil(search_eutil(db, term)))
        if len(accs) > 0:
            plural = 'sequences'
            if len(accs) == 1:
                plural = 'sequence'
            Log.msg('Found {} RefSeq {} {} for'.format(len(accs), filter_term,
                                                       plural), tax_term)
            # ToDo: Using a random sample of X RefSeq sequences.
            # Random sample ###################################################
            # if len(accs) > 10:
            #     Log.wrn('Using a random sample of ten RefSeq sequences.')
            #     random.seed(a=len(accs), version=2)
            #     accs = set(random.sample(population=accs, k=10))
            ###################################################################
            break
        else:
            Log.wrn('No RefSeq {} sequences were found for'.format(
                filter_term), tax_term)

    cache_path = opj(dir_cache_refseqs, filter_term + '__' +
                     tax_term.replace(' ', '_') + '.fasta')
    parsed_fasta_cache = {}
    if ope(cache_path):
        parsed_fasta_cache = read_fasta(cache_path, seq_type=SEQ_TYPE_NT,
                                        def_to_first_space=True)
        parsed_fasta_cache = seq_records_to_dict(parsed_fasta_cache)
        for acc in parsed_fasta_cache:
            if acc in accs:
                accs.remove(acc)
    if len(accs) > 0:
        parsed_fasta = dnld_ncbi_seqs(db, list(accs))
        parsed_fasta = seq_records_to_dict(parsed_fasta, prepend_acc=True)
        parsed_fasta.update(parsed_fasta_cache)
        write_fasta(parsed_fasta, cache_path)

    return cache_path


def _should_run_bt2(taxid, taxonomy, bt2_order, bowtie2, bowtie2_build):

    dbs = OrderedDict()

    for x in bt2_order:
        db_path_ok = False

        if x == MT:
            if taxonomy.is_eukaryote(taxid) is True:
                if bt2_order[MT] == '':
                    dbs[MT] = MT
                    db_path_ok = True

        elif x == PT:
            if taxonomy.is_eukaryote(taxid) is True:
                if taxonomy.contains_plastid(taxid) is True:
                    if bt2_order[PT] == '':
                        dbs[PT] = PT
                        db_path_ok = True

        if db_path_ok is False:
            db_path = bt2_order[x]
            if ope(db_path) and isfile(db_path):
                dbs[x] = db_path
            else:
                Log.err('File not found:', db_path)
                exit(1)

    if len(dbs) > 0:

        if bowtie2 is None:
            Log.err('bowtie2 is not available. ' +
                    'Cannot continue. Exiting.')
            exit(0)

        if bowtie2_build is None:
            Log.err('bowtie2-build is not available. ' +
                    'Cannot continue. Exiting.')
            exit(0)

    return dbs


def run_bt2_fq(se_fastq_files, pe_fastq_files, dir_fq_filter_data,
               bowtie2, bowtie2_build, threads, dir_temp, bt2_order,
               fpatt, taxonomy, dir_cache_refseqs,
               rcorrector_before_trimmomatic, gz_out=True):

    new_se_fastq_files = dict()
    new_pe_fastq_files = dict()

    msg_printed = False

    ext_out = ''
    if gz_out is True:
        ext_out = '.gz'

    # SE
    for se in se_fastq_files:

        taxid = se_fastq_files[se]['tax_id']
        dbs = _should_run_bt2(taxid, taxonomy, bt2_order, bowtie2,
                              bowtie2_build)

        if rcorrector_before_trimmomatic is True:
            in_f = se_fastq_files[se]['trim_path_fq']
        else:
            in_f = se_fastq_files[se]['cor_path_fq']

        in_f_orig = in_f

        if len(dbs) == 0:
            se_fastq_files[se]['filter_path_fq'] = in_f
            continue

        if msg_printed is False:
            print()
            Log.inf('Running Bowtie2.')
            msg_printed = True

        dir_fq_bt_data_sample_un = opj(dir_fq_filter_data, se)

        for i, db in enumerate(dbs):

            db_path = dbs[db]

            dir_fq_bt_data_sample = opj(dir_fq_filter_data, se, db)

            new_se = se + '_' + db

            out_f = opj(dir_fq_bt_data_sample, new_se + '.fastq' + ext_out)

            out_f_un = opj(dir_temp, new_se + '_bt2_unaligned' + '.fastq' + ext_out)

            sam_f = opj(dir_fq_bt_data_sample, new_se + '.sam')
            new_se_fastq_files[new_se] = deepcopy(se_fastq_files[se])
            new_se_fastq_files[new_se]['path'] = None
            new_se_fastq_files[new_se]['cor_path_fq'] = None
            new_se_fastq_files[new_se]['trim_path_fq'] = None
            taxid = new_se_fastq_files[new_se]['tax_id']
            gc = new_se_fastq_files[new_se]['gc_id']
            new_se_fastq_files[new_se]['organelle'] = None
            if db == MT:
                gc = taxonomy.mito_genetic_code_for_taxid(taxid)
                new_se_fastq_files[new_se]['gc_id'] = gc
                new_se_fastq_files[new_se]['organelle'] = MT
            elif db == PT:
                gc = taxonomy.plastid_genetic_code_for_taxid(taxid)
                new_se_fastq_files[new_se]['gc_id'] = gc
                new_se_fastq_files[new_se]['organelle'] = PT
            new_se_fastq_files[new_se]['gc_tt'] = TranslationTable(gc)
            new_se_fastq_files[new_se]['filter_path_fq'] = out_f
            if ope(dir_fq_bt_data_sample):
                Log.msg('Bowtie2 filtered FASTQ file already exists:', new_se)
                in_f = opj(dir_fq_bt_data_sample_un, se + '.fastq' + ext_out)
            else:
                Log.msg('SE mode:', new_se)
                make_dirs(dir_fq_bt_data_sample)

                if db_path in (MT, PT):
                    db_fasta_path = dnld_refseqs_for_taxid(
                        taxid, db, taxonomy, dir_cache_refseqs, query='',
                        db='nuccore')
                    bt2_idx_path = splitext(db_fasta_path)[0]
                else:
                    db_fasta_path = db_path
                    bt2_idx_path = opj(dir_cache_refseqs,
                                       splitext(basename(db_fasta_path))[0])

                if not ope(bt2_idx_path + '.1.bt2'):
                    build_bt2_index(bowtie2_build, [db_fasta_path],
                                    bt2_idx_path, threads)

                run_bowtie2_se(bowtie2=bowtie2,
                               input_file=in_f,
                               output_file=out_f,
                               output_file_un=out_f_un,
                               sam_output_file=sam_f,
                               index=bt2_idx_path,
                               threads=threads,
                               dir_temp=dir_temp,
                               gz_out=gz_out)

                if i > 0:
                    remove(in_f)

                in_f = out_f_un

        out_f_un = opj(dir_fq_bt_data_sample_un, se + '.fastq' + ext_out)
        se_fastq_files[se]['filter_path_fq'] = out_f_un

        if in_f != in_f_orig:
            move(in_f, out_f_un)

    se_fastq_files.update(new_se_fastq_files)

    # PE
    for pe in pe_fastq_files:

        taxid = pe_fastq_files[pe]['tax_id']
        dbs = _should_run_bt2(taxid, taxonomy, bt2_order, bowtie2,
                              bowtie2_build)

        if rcorrector_before_trimmomatic is True:
            in_fs = pe_fastq_files[pe]['trim_path_fq']
        else:
            in_fs = pe_fastq_files[pe]['cor_path_fq']

        in_fs_orig = tuple(in_fs)

        if len(dbs) == 0:
            pe_fastq_files[pe]['filter_path_fq'] = in_fs
            continue

        if msg_printed is False:
            print()
            Log.inf('Running Bowtie2.')
            msg_printed = True

        dir_fq_bt_data_sample_un = opj(dir_fq_filter_data, pe)

        for i, db in enumerate(dbs):

            db_path = dbs[db]

            dir_fq_bt_data_sample = opj(dir_fq_filter_data, pe, db)

            new_pe = pe + '_' + db

            out_fs = [x.replace('@D@', dir_fq_bt_data_sample) for x in fpatt]
            out_fs = [x.replace('@N@', new_pe) for x in out_fs]
            out_fs = [x + ext_out for x in out_fs]

            out_fs_un = [x.replace('@D@', dir_temp) for x in fpatt]
            out_fs_un = [x.replace('@N@', new_pe + '_bt2_unaligned')
                         for x in out_fs_un]
            out_fs_un = [x + ext_out for x in out_fs_un]

            sam_f = opj(dir_fq_bt_data_sample, new_pe + '.sam')
            new_pe_fastq_files[new_pe] = deepcopy(pe_fastq_files[pe])
            new_pe_fastq_files[new_pe]['path'] = None
            new_pe_fastq_files[new_pe]['cor_path_fq'] = None
            new_pe_fastq_files[new_pe]['trim_path_fq'] = None
            taxid = new_pe_fastq_files[new_pe]['tax_id']
            gc = new_pe_fastq_files[new_pe]['gc_id']
            new_pe_fastq_files[new_pe]['organelle'] = None
            if db == MT:
                gc = taxonomy.mito_genetic_code_for_taxid(taxid)
                new_pe_fastq_files[new_pe]['gc_id'] = gc
                new_pe_fastq_files[new_pe]['organelle'] = MT
            elif db == PT:
                gc = taxonomy.plastid_genetic_code_for_taxid(taxid)
                new_pe_fastq_files[new_pe]['gc_id'] = gc
                new_pe_fastq_files[new_pe]['organelle'] = PT
            new_pe_fastq_files[new_pe]['gc_tt'] = TranslationTable(gc)
            new_pe_fastq_files[new_pe]['filter_path_fq'] = out_fs
            if ope(dir_fq_bt_data_sample):
                Log.msg('Bowtie2 filtered FASTQ files already exist:', new_pe)
                in_fs = [x.replace('@D@', dir_fq_bt_data_sample_un)
                         for x in fpatt]
                in_fs = [x.replace('@N@', pe) for x in in_fs]
                in_fs = [x + ext_out for x in in_fs]
            else:
                Log.msg('PE mode:', new_pe)
                make_dirs(dir_fq_bt_data_sample)

                if db_path in (MT, PT):
                    db_fasta_path = dnld_refseqs_for_taxid(
                        taxid, db, taxonomy, dir_cache_refseqs, query='',
                        db='nuccore')
                    bt2_idx_path = splitext(db_fasta_path)[0]
                else:
                    db_fasta_path = db_path
                    bt2_idx_path = opj(dir_cache_refseqs,
                                       splitext(basename(db_fasta_path))[0])

                if not ope(bt2_idx_path + '.1.bt2'):
                    build_bt2_index(bowtie2_build, [db_fasta_path],
                                    bt2_idx_path,
                                    threads)

                paired_out_pattern = out_fs[0].replace(
                    '_paired_1.fastq', '_paired_%.fastq')

                paired_out_pattern_un = out_fs_un[0].replace(
                    '_paired_1.fastq', '_paired_%.fastq')

                run_bowtie2_pe(bowtie2=bowtie2,
                               input_files=in_fs,
                               paired_out_pattern=paired_out_pattern,
                               paired_out_pattern_un=paired_out_pattern_un,
                               unpaired_out_1=out_fs[2],
                               unpaired_out_2=out_fs[3],
                               unpaired_out_1_un=out_fs_un[2],
                               unpaired_out_2_un=out_fs_un[3],
                               sam_output_file=sam_f,
                               index=bt2_idx_path,
                               threads=threads,
                               dir_temp=dir_temp,
                               gz_out=gz_out)

                if i > 0:
                    remove(in_fs[0])
                    remove(in_fs[1])
                    remove(in_fs[2])
                    remove(in_fs[3])

                in_fs = out_fs_un

        out_fs_un = [x.replace('@D@', dir_fq_bt_data_sample_un) for x in fpatt]
        out_fs_un = [x.replace('@N@', pe) for x in out_fs_un]
        out_fs_un = [x + ext_out for x in out_fs_un]
        pe_fastq_files[pe]['filter_path_fq'] = out_fs_un

        if tuple(in_fs) != in_fs_orig:
            move(in_fs[0], out_fs_un[0])
            move(in_fs[1], out_fs_un[1])
            move(in_fs[2], out_fs_un[2])
            move(in_fs[3], out_fs_un[3])

    pe_fastq_files.update(new_pe_fastq_files)


def run_kraken2(order, dbs, se_fastq_files, pe_fastq_files, dir_fq_filter_data,
                confidence, kraken2, threads, dir_temp, fpatt):

    if (len(se_fastq_files) > 0 or len(pe_fastq_files) > 0) and len(order) > 0:
        print()
        Log.inf('Running Kraken2.', 'Confidence: ' + str(confidence))
        if kraken2 is None:
            Log.err('kraken2 is not available. Cannot continue. Exiting.')
            exit(0)

    nuclear = None
    for nuc in order:
        if nuc[1] == 'nuclear':
            nuclear = nuc[0]
            break

    for se in se_fastq_files:

        if len(order) == 0:
            continue

        if se_fastq_files[se]['path'] is None:
            continue

        fq_path = se_fastq_files[se]['filter_path_fq']
        dir_fq_filter_data_sample = opj(dir_fq_filter_data, se)

        if nuclear is None:
            out_f = opj(dir_fq_filter_data_sample, se + '.fastq')
        else:
            out_f = opj(dir_fq_filter_data_sample, nuclear, se + '.fastq')

        se_fastq_files[se]['filter_path_fq'] = out_f

        if ope(dir_fq_filter_data_sample):
            Log.msg('Kraken2 filtered FASTQ files already exist:', se)
        else:
            make_dirs(dir_fq_filter_data_sample)
            print()
            Log.msg('SE mode:', se)
            run_kraken_filters(
                order=order,
                dbs=dbs,
                base_name=se,
                in_files=fq_path,
                dir_out=dir_fq_filter_data_sample,
                confidence=confidence,
                kraken2=kraken2,
                threads=threads,
                dir_temp=dir_temp)

    for pe in pe_fastq_files:

        if len(order) == 0:
            continue

        if pe_fastq_files[pe]['path'] is None:
            continue

        fq_path = pe_fastq_files[pe]['filter_path_fq']
        dir_fq_filter_data_sample = opj(dir_fq_filter_data, pe)

        if nuclear is None:
            dir_name_nuclear = dir_fq_filter_data_sample
        else:
            dir_name_nuclear = dir_fq_filter_data_sample + ops + nuclear

        out_fs = [x.replace('@D@', dir_name_nuclear) for x in fpatt]
        out_fs = [x.replace('@N@', pe) for x in out_fs]

        pe_fastq_files[pe]['filter_path_fq'] = out_fs

        if ope(dir_fq_filter_data_sample):
            Log.msg('Kraken2 filtered FASTQ files already exist:', pe)
        else:
            make_dirs(dir_fq_filter_data_sample)
            print()
            Log.msg('PE mode:', pe)
            run_kraken_filters(
                order=order,
                dbs=dbs,
                base_name=pe,
                in_files=fq_path,
                dir_out=dir_fq_filter_data_sample,
                confidence=confidence,
                kraken2=kraken2,
                threads=threads,
                dir_temp=dir_temp)


def filtered_fq_to_fa(se_fastq_files, pe_fastq_files, dir_fa_trim_data, seqtk,
                      fpatt):
    if len(se_fastq_files) > 0 or len(pe_fastq_files) > 0:
        print()
        Log.inf('Converting FASTQ to FASTA using Seqtk.')
        if seqtk is None:
            Log.err('seqtk is not available. Cannot continue. Exiting.')
            exit(0)
    for se in se_fastq_files:
        dir_fa_trim_data_sample = opj(dir_fa_trim_data, se)
        fq_path = se_fastq_files[se]['filter_path_fq']
        out_f = opj(dir_fa_trim_data_sample, se + '.fasta')
        se_fastq_files[se]['filter_path_fa'] = out_f

        if ope(dir_fa_trim_data_sample):
            Log.msg('Filtered FASTA files already exist:', se)
        else:
            make_dirs(dir_fa_trim_data_sample)
            Log.msg(basename(fq_path))
            seqtk_fq_to_fa(seqtk, fq_path, out_f)

    for pe in pe_fastq_files:
        dir_fa_trim_data_sample = opj(dir_fa_trim_data, pe)
        fq_paths = pe_fastq_files[pe]['filter_path_fq']
        out_fs = [x.replace('@D@', dir_fa_trim_data_sample) for x in fpatt]
        out_fs = [x.replace('@N@', pe) for x in out_fs]
        pe_fastq_files[pe]['filter_path_fa'] = out_fs

        if ope(dir_fa_trim_data_sample):
            Log.msg('Filtered FASTA files already exist:', pe)
        else:
            make_dirs(dir_fa_trim_data_sample)
            pe_trim_files = zip(fq_paths, out_fs)
            for x in pe_trim_files:
                Log.msg(basename(x[0]))
                seqtk_fq_to_fa(seqtk, x[0], x[1])


def makeblastdb_fq(se_fastq_files, pe_fastq_files, dir_blast_fa_trim,
                   makeblastdb, fpatt):
    for se in se_fastq_files:
        dir_blast_fa_trim_sample = opj(dir_blast_fa_trim, se)
        fa_path = se_fastq_files[se]['filter_path_fa']
        out_f = opj(dir_blast_fa_trim_sample, se)
        se_fastq_files[se]['blast_db_path'] = out_f

        if ope(dir_blast_fa_trim_sample):
            Log.msg('BLAST database already exists:', se)
        else:
            make_dirs(dir_blast_fa_trim_sample)
            Log.msg(basename(fa_path))
            make_blast_db(
                exec_file=makeblastdb,
                in_file=fa_path,
                out_file=out_f,
                title=se,
                dbtype='nucl')

    for pe in pe_fastq_files:
        dir_blast_fa_trim_sample = opj(dir_blast_fa_trim, pe)
        fa_paths = pe_fastq_files[pe]['filter_path_fa']
        out_fs = [x.replace('@D@', dir_blast_fa_trim_sample) for x in fpatt]
        out_fs = [x.replace('@N@', pe) for x in out_fs]
        pe_fastq_files[pe]['blast_db_path'] = out_fs

        if ope(dir_blast_fa_trim_sample):
            Log.msg('BLAST database already exists:', pe)
        else:
            make_dirs(dir_blast_fa_trim_sample)
            pe_trim_files = zip(fa_paths, out_fs)
            for x in pe_trim_files:
                Log.msg(basename(x[0]))
                make_blast_db(
                    exec_file=makeblastdb,
                    in_file=x[0],
                    out_file=x[1],
                    title=basename(x[1]),
                    dbtype='nucl')
