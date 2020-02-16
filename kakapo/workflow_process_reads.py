#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Kakapo workflow: Process Reads."""

import fileinput
import pickle
import re

from copy import deepcopy
from os import remove as osremove
from os.path import basename
from os.path import commonprefix
from os.path import exists as ope
from os.path import join as opj
from os.path import sep as ops
from os.path import splitext
from shutil import copyfile
from shutil import rmtree
from time import sleep

from kakapo.bioio import read_fasta
from kakapo.bioio import write_fasta
from kakapo.blast import BLST_RES_COLS_1
from kakapo.blast import make_blast_db, run_blast
from kakapo.bowtie2 import build_bt2_index, run_bowtie2_se, run_bowtie2_pe
from kakapo.config import PICKLE_PROTOCOL
from kakapo.entrez import accessions as accessions_ncbi
from kakapo.entrez import dnld_seqs_fasta_format
from kakapo.entrez import esearch
from kakapo.entrez import sra_run_info
from kakapo.helpers import combine_text_files
from kakapo.helpers import keep_unique_lines_in_file
from kakapo.helpers import make_dir
from kakapo.helpers import plain_or_gzip
from kakapo.helpers import splitext_gz
from kakapo.kraken import run_kraken_filters
from kakapo.rcorrector import filter_unc_se, filter_unc_pe
from kakapo.rcorrector import run_rcorrector_se, run_rcorrector_pe
from kakapo.seqtk import seqtk_fq_to_fa, seqtk_extract_reads
from kakapo.shell import call
from kakapo.translation_tables import TranslationTable
from kakapo.trimmomatic import trimmomatic_se, trimmomatic_pe
from kakapo.vsearch import run_cluster_fast, run_vsearch


def dnld_sra_info(sras, dir_cache_prj, linfo=print):  # noqa
    sra_runs_info = {}
    sras_acceptable = []

    if len(sras) > 0:
        linfo('Downloading SRA run information')
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
            sra_species = info['ScientificName']
            sra_taxid = info['TaxID']
            sra_spots = int(info['spots'])
            sra_spots_with_mates = int(info['spots_with_mates'])

            sample_base_name = (sra_species.replace(' ', '_') + '_' +
                                sra_taxid + '_' + sra)

            sra_runs_info[sra]['KakapoSampleBaseName'] = sample_base_name

            if sra_lib_source != 'transcriptomic':
                sra_info_str = (
                    '{sra}: the SRA library source type "{ltype}" '
                    'is not supported').format(
                    sra=sra, ltype=sra_lib_source)

            elif sra_seq_platform != 'Illumina':
                sra_info_str = (
                    '{sra}: the SRA library sequencing platform "{plat}" '
                    'is not supported').format(
                    sra=sra, plat=sra_seq_platform)

            else:
                sra_info_str = ('SRA run {sra} {source} '
                                '{layout}-end library. '
                                'Sourced from {species} '
                                '(TaxID: {txid}). '
                                'Sequenced using {platform} platform on '
                                '{model}.').format(
                                    sra=sra,
                                    source=sra_lib_source.title(),
                                    strategy=sra_lib_strategy,
                                    layout=sra_lib_layout,
                                    platform=sra_seq_platform,
                                    model=sra_seq_platform_model,
                                    species=sra_species,
                                    txid=sra_taxid)

                sra_runs_info[sra]['KakapoLibraryLayout'] = \
                    sra_runs_info[sra]['LibraryLayout']

                if sra_lib_layout == 'paired' and sra_spots_with_mates == 0:
                    sra_runs_info[sra]['KakapoLibraryLayout'] = 'SINGLE'
                    sra_info_str = (
                        sra_info_str + ' Listed as containing '
                        'paired-end reads, but only a single set of reads '
                        'is available. Treating as single-ended.')

                elif (sra_lib_layout == 'paired' and
                      sra_spots != sra_spots_with_mates):
                    sra_runs_info[sra]['KakapoLibraryLayout'] = 'PAIRED_UNP'
                    sra_info_str = (
                        sra_info_str + ' Listed as containing '
                        'paired-end reads, but not all reads are paired.')

                sras_acceptable.append(sra)

            linfo(sra_info_str)

    with open(__, 'wb') as f:
        pickle.dump(sra_runs_info, f, protocol=PICKLE_PROTOCOL)

    return sra_runs_info, sras_acceptable


def dnld_sra_fastq_files(sras, sra_runs_info, dir_fq_data, fasterq_dump,
                         threads, dir_temp, linfo=print): # noqa

    if len(sras) > 0:
        if fasterq_dump is None:
            linfo('fasterq-dump from SRA Toolkit is not available. ' +
                  'Cannot continue. Exiting.')
            exit(0)

    se_fastq_files = {}
    pe_fastq_files = {}

    for sra in sras:
        sra_run_info = sra_runs_info[sra]
        sra_lib_layout = sra_run_info['LibraryLayout'].lower()
        sra_lib_layout_k = sra_run_info['KakapoLibraryLayout'].lower()
        sample_base_name = sra_run_info['KakapoSampleBaseName']
        sra_taxid = int(sra_run_info['TaxID'])
        avg_len = int(sra_run_info['avgLength'])

        sra_dnld_needed = False

        if sra_lib_layout == 'single' or sra_lib_layout_k == 'single':
            se_file = opj(dir_fq_data, sra + '.fastq')
            se_fastq_files[sample_base_name] = {'path': se_file}
            se_fastq_files[sample_base_name]['src'] = 'sra'
            se_fastq_files[sample_base_name]['avg_len'] = avg_len
            se_fastq_files[sample_base_name]['tax_id'] = sra_taxid
            if not ope(se_file):
                sra_dnld_needed = True

        elif sra_lib_layout == 'paired':
            pe_file_1 = opj(dir_fq_data, sra + '_1.fastq')
            pe_file_2 = opj(dir_fq_data, sra + '_2.fastq')
            pe_fastq_files[sample_base_name] = {'path': [pe_file_1, pe_file_2]}
            pe_fastq_files[sample_base_name]['src'] = 'sra'
            pe_fastq_files[sample_base_name]['avg_len'] = avg_len // 2
            pe_fastq_files[sample_base_name]['tax_id'] = sra_taxid
            if sra_lib_layout_k == 'paired_unp':
                pe_file_3 = opj(dir_fq_data, sra + '.fastq')
                pe_fastq_files[sample_base_name]['path'].append(pe_file_3)
            if not ope(pe_file_1) or not ope(pe_file_2):
                sra_dnld_needed = True

        if not sra_dnld_needed:
            linfo('FASTQ reads for the SRA run ' + sample_base_name +
                  ' are available locally')

        retry_count = 0
        while sra_dnld_needed:

            if retry_count > 50:
                linfo('Download failed. Exiting.')
                rmtree(dir_temp)
                exit(1)

            elif retry_count > 0:
                linfo('Download failed. Retrying.')
                sleep(2)

            retry_count += 1

            linfo('Downloading FASTQ reads for ' + sample_base_name)

            cmd = [fasterq_dump,
                   '--threads', str(threads * 4),
                   '--split-3',
                   '--bufsize', '819200',
                   '--outdir', dir_fq_data,
                   '--temp', dir_temp, sra]

            call(cmd)

            if sra_lib_layout == 'single' or sra_lib_layout_k == 'single':
                if not ope(se_file):
                    continue

            elif sra_lib_layout == 'paired':
                if not ope(pe_file_1) or not ope(pe_file_2):
                    continue

            sra_dnld_needed = False

    return se_fastq_files, pe_fastq_files, sra_runs_info


def user_fastq_files(fq_se, fq_pe, linfo=print): # noqa
    if len(fq_se) > 0 or len(fq_pe) > 0:
        linfo('Preparing user provided FASTQ files')

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
        sample_base_name = base
        se_fastq_files[sample_base_name] = {'path': path}
        se_fastq_files[sample_base_name]['src'] = 'usr'
        se_fastq_files[sample_base_name]['avg_len'] = None
        se_fastq_files[sample_base_name]['tax_id'] = tax_id
        linfo(sample_base_name + ': ' + path)

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
        sample_base_name = base
        pe_fastq_files[sample_base_name] = {'path': path}
        pe_fastq_files[sample_base_name]['src'] = 'usr'
        pe_fastq_files[sample_base_name]['avg_len'] = None
        pe_fastq_files[sample_base_name]['tax_id'] = tax_id
        linfo(sample_base_name + ': ' + path[0] + ', ' + path[1])

    return se_fastq_files, pe_fastq_files


def min_accept_read_len(se_fastq_files, pe_fastq_files, dir_temp,
                        dir_cache_fq_minlen, vsearch, linfo=print): # noqa
    # lowest allowable
    low = 35

    if len(se_fastq_files) > 0 or len(pe_fastq_files) > 0:
        linfo('Calculating minimum acceptable read length')
        if vsearch is None:
            linfo('vsearch is not available. ' +
                  'Cannot continue. Exiting.')
            exit(0)
    else:
        return None

    __ = opj(dir_cache_fq_minlen, 'minlen')

    pickled = {}

    if ope(__):
        with open(__, 'rb') as f:
            pickled = pickle.load(f)

    queue = []

    for se in se_fastq_files:
        src = se_fastq_files[se]['src']
        avg_len = se_fastq_files[se]['avg_len']
        if src == 'sra':
            ml = max(avg_len // 3, low)
            se_fastq_files[se]['min_acc_len'] = ml
            linfo(str(ml) + ' nt: ' + se)
            continue

        fq_path = se_fastq_files[se]['path']
        stats_file = opj(dir_temp, se + '_stats.txt')
        queue.append([se, fq_path, stats_file, 'se'])

    for pe in pe_fastq_files:
        src = pe_fastq_files[pe]['src']
        avg_len = pe_fastq_files[pe]['avg_len']
        if src == 'sra':
            ml = max(avg_len // 3, low)
            pe_fastq_files[pe]['min_acc_len'] = ml
            linfo(str(ml) + ' nt: ' + pe)
            continue

        fq_path = pe_fastq_files[pe]['path'][0]
        stats_file = opj(dir_temp, pe + '_stats.txt')
        queue.append([pe, fq_path, stats_file, 'pe'])

    for x in queue:

        if x[0] in pickled:
            ml = pickled[x[0]]

        else:
            cmd = [vsearch, '--fastq_stats', x[1], '--log', x[2]]
            call(cmd)

            with open(x[2]) as f:
                stats = f.read()

            osremove(x[2])

            ml = re.findall(r'>=\s+(\d+)', stats)

            if len(ml) != 0:
                ml = max(int(ml[0]) // 3, low)
            else:
                ml = None

            pickled[x[0]] = ml

        if ml is not None:
            linfo(str(ml) + ' nt: ' + x[0])
        else:
            linfo(' ?' + ' nt: ' + x[0])
            ml = low

        if x[3] == 'se':
            se_fastq_files[x[0]]['min_acc_len'] = ml

        elif x[3] == 'pe':
            pe_fastq_files[x[0]]['min_acc_len'] = ml

        with open(__, 'wb') as f:
            pickle.dump(pickled, f, protocol=PICKLE_PROTOCOL)


def run_rcorrector(se_fastq_files, pe_fastq_files, dir_fq_cor_data, rcorrector,
                   threads, dir_temp, linfo=print):  # noqa
    if len(se_fastq_files) > 0 or len(pe_fastq_files) > 0:
        if rcorrector is None:
            linfo('Rcorrector is not available. ' +
                  'Cannot continue. Exiting.')
            exit(0)
    for se in se_fastq_files:
        dir_fq_cor_data_sample = opj(dir_fq_cor_data, se)
        fq_path = se_fastq_files[se]['path']
        r_mode, w_mode, a_mode, fqopen, ext = plain_or_gzip(fq_path)
        log_f = opj(dir_fq_cor_data_sample, se + '.txt')
        out_f = opj(dir_fq_cor_data_sample, se + '.fastq' + ext)
        se_fastq_files[se]['cor_path_fq'] = out_f

        if ope(dir_fq_cor_data_sample):
            linfo('Corrected FASTQ file for sample ' + se + ' already exists')
        else:
            make_dir(dir_fq_cor_data_sample)
            linfo('Running Rcorrector in SE mode: ' + se)
            run_rcorrector_se(rcorrector=rcorrector,
                              in_file=fq_path,
                              out_dir=dir_fq_cor_data_sample,
                              threads=threads,
                              dir_temp=dir_temp)

            fq_base_path = opj(dir_fq_cor_data_sample, basename(fq_path))
            fq_cor_path = splitext_gz(fq_base_path)[0] + '.cor.fq' + ext

            filter_unc_se(in_file=fq_cor_path, out_file=out_f, log_file=log_f)

            osremove(fq_cor_path)

    for pe in pe_fastq_files:
        dir_fq_cor_data_sample = opj(dir_fq_cor_data, pe)
        fq_path_1 = pe_fastq_files[pe]['path'][0]
        fq_path_2 = pe_fastq_files[pe]['path'][1]
        fq_path_3 = None
        out_f_3 = None
        r_mode, w_mode, a_mode, fqopen, ext = plain_or_gzip(fq_path_1)
        log_f = opj(dir_fq_cor_data_sample, pe + '.txt')
        out_f_1 = opj(dir_fq_cor_data_sample, pe + '_R1.fastq' + ext)
        out_f_2 = opj(dir_fq_cor_data_sample, pe + '_R2.fastq' + ext)
        pe_fastq_files[pe]['cor_path_fq'] = [out_f_1, out_f_2]

        if len(pe_fastq_files[pe]['path']) == 3:
            fq_path_3 = pe_fastq_files[pe]['path'][2]
            out_f_3 = opj(dir_fq_cor_data_sample, pe + '_R3.fastq' + ext)
            pe_fastq_files[pe]['cor_path_fq'].append(out_f_3)

        if ope(dir_fq_cor_data_sample):
            linfo('Corrected FASTQ files for sample ' + pe + ' already exist')
        else:
            make_dir(dir_fq_cor_data_sample)
            linfo('Running Rcorrector in PE mode: ' + pe)
            run_rcorrector_pe(rcorrector=rcorrector,
                              in_file_1=fq_path_1,
                              in_file_2=fq_path_2,
                              out_dir=dir_fq_cor_data_sample,
                              threads=threads,
                              dir_temp=dir_temp)

            fq_base_path_1 = opj(dir_fq_cor_data_sample, basename(fq_path_1))
            fq_cor_path_1 = splitext_gz(fq_base_path_1)[0] + '.cor.fq' + ext
            fq_base_path_2 = opj(dir_fq_cor_data_sample, basename(fq_path_2))
            fq_cor_path_2 = splitext_gz(fq_base_path_2)[0] + '.cor.fq' + ext

            filter_unc_pe(in_file_1=fq_cor_path_1,
                          in_file_2=fq_cor_path_2,
                          out_file_1=out_f_1,
                          out_file_2=out_f_2,
                          log_file=log_f)

            osremove(fq_cor_path_1)
            osremove(fq_cor_path_2)

            if fq_path_3 is not None:

                linfo('Running Rcorrector in SE mode: ' + pe +
                      ' (Paired-read SRA run contains unpaired reads.)')

                run_rcorrector_se(rcorrector=rcorrector,
                                  in_file=fq_path_3,
                                  out_dir=dir_fq_cor_data_sample,
                                  threads=threads,
                                  dir_temp=dir_temp)

                fq_base_path_3 = opj(dir_fq_cor_data_sample,
                                     basename(fq_path_3))
                fq_cor_path_3 = splitext_gz(fq_base_path_3)[0] + '.cor.fq'
                log_f_3 = opj(dir_fq_cor_data_sample, pe + '_unpaired.txt')

                filter_unc_se(in_file=fq_cor_path_3, out_file=out_f_3,
                              log_file=log_f_3)

                osremove(fq_cor_path_3)


def run_trimmomatic(se_fastq_files, pe_fastq_files, dir_fq_trim_data,
                    trimmomatic, adapters, fpatt, threads, linfo=print):  # noqa
    if len(se_fastq_files) > 0 or len(pe_fastq_files) > 0:
        if trimmomatic is None:
            linfo('trimmomatic is not available. ' +
                  'Cannot continue. Exiting.')
            exit(0)
    for se in se_fastq_files:
        dir_fq_trim_data_sample = opj(dir_fq_trim_data, se)
        fq_path = se_fastq_files[se]['cor_path_fq']
        r_mode, w_mode, a_mode, fqopen, ext = plain_or_gzip(fq_path)
        min_acc_len = se_fastq_files[se]['min_acc_len']
        stats_f = opj(dir_fq_trim_data_sample, se + '.txt')
        out_f = opj(dir_fq_trim_data_sample, se + '.fastq' + ext)
        se_fastq_files[se]['trim_path_fq'] = out_f

        if ope(dir_fq_trim_data_sample):
            linfo('Trimmed FASTQ file for sample ' + se + ' already exists')
        else:
            make_dir(dir_fq_trim_data_sample)
            linfo('Running Trimmomatic in SE mode: ' + se)
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
        fq_path_1 = pe_fastq_files[pe]['cor_path_fq'][0]
        fq_path_2 = pe_fastq_files[pe]['cor_path_fq'][1]
        fq_path_3 = None
        r_mode, w_mode, a_mode, fqopen, ext = plain_or_gzip(fq_path_1)
        if len(pe_fastq_files[pe]['cor_path_fq']) == 3:
            fq_path_3 = pe_fastq_files[pe]['cor_path_fq'][2]
        min_acc_len = pe_fastq_files[pe]['min_acc_len']
        stats_f = opj(dir_fq_trim_data_sample, pe + '.txt')
        out_fs = [x.replace('@D@', dir_fq_trim_data_sample) for x in fpatt]
        out_fs = [x.replace('@N@', pe) for x in out_fs]
        out_fs = [x + ext for x in out_fs]
        pe_fastq_files[pe]['trim_path_fq'] = out_fs

        if ope(dir_fq_trim_data_sample):
            linfo('Trimmed FASTQ files for sample ' + pe + ' already exist')
        else:
            make_dir(dir_fq_trim_data_sample)
            linfo('Running Trimmomatic in PE mode: ' + pe)
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

                out_f = opj(dir_fq_trim_data_sample, 'unpaired.fastq' + ext)
                stats_f = opj(dir_fq_trim_data_sample, pe + '_unpaired.txt')

                linfo('Running Trimmomatic in SE mode: ' + pe +
                      ' (Paired-read SRA run contains unpaired reads.)')

                trimmomatic_se(
                    trimmomatic=trimmomatic,
                    adapters=adapters,
                    in_file=fq_path_3,
                    out_file=out_f,
                    stats_file=stats_f,
                    threads=threads,
                    minlen=min_acc_len)

                _ = opj(dir_fq_trim_data_sample, 'temp.fastq' + ext)
                f_temp = fqopen(_, a_mode)
                fi = fileinput.FileInput(openhook=fileinput.hook_compressed)
                with fi.input(files=[out_fs[2], out_f]) as f:
                    for line in f:
                        f_temp.write(line)
                f_temp.close()

                osremove(out_fs[2])
                osremove(out_f)
                copyfile(_, out_fs[2])
                osremove(_)


def run_kraken2(order, dbs, se_fastq_files, pe_fastq_files, dir_fq_filter_data,
                confidence, kraken2, threads, dir_temp, fpatt, linfo=print):  # noqa
    if (len(se_fastq_files) > 0 or len(pe_fastq_files) > 0) and len(order) > 0:
        if kraken2 is None:
            linfo('kraken2 is not available. ' +
                  'Cannot continue. Exiting.')
            exit(0)
    nuclear = ''
    for nuc in order:
        if nuc[1] == 'nuclear':
            nuclear = nuc[0]
            break

    for se in se_fastq_files:
        fq_path = se_fastq_files[se]['trim_path_fq']

        if len(order) == 0:
            se_fastq_files[se]['filter_path_fq'] = fq_path
            continue

        dir_fq_filter_data_sample = opj(dir_fq_filter_data, se)
        out_f = opj(dir_fq_filter_data_sample, nuclear, se + '.fastq')
        se_fastq_files[se]['filter_path_fq'] = out_f
        if ope(dir_fq_filter_data_sample):
            linfo('Kraken2 filtered FASTQ file for sample ' + se +
                  ' already exists')
        else:
            make_dir(dir_fq_filter_data_sample)
            linfo('Running Kraken2 in SE mode: ' + se)
            run_kraken_filters(
                order=order,
                dbs=dbs,
                base_name=se,
                in_files=fq_path,
                dir_out=dir_fq_filter_data_sample,
                confidence=confidence,
                kraken2=kraken2,
                threads=threads,
                dir_temp=dir_temp,
                linfo=linfo)

    for pe in pe_fastq_files:
        fq_path = pe_fastq_files[pe]['trim_path_fq']

        if len(order) == 0:
            pe_fastq_files[pe]['filter_path_fq'] = fq_path
            continue

        dir_fq_filter_data_sample = opj(dir_fq_filter_data, pe)
        dir_name_nuclear = dir_fq_filter_data_sample + ops + nuclear
        out_fs = [x.replace('@D@', dir_name_nuclear) for x in fpatt]
        out_fs = [x.replace('@N@', pe) for x in out_fs]
        pe_fastq_files[pe]['filter_path_fq'] = out_fs
        if ope(dir_fq_filter_data_sample):
            linfo('Kraken2 filtered FASTQ files for sample ' + pe +
                  ' already exist')
        else:
            make_dir(dir_fq_filter_data_sample)
            linfo('Running Kraken2 in PE mode: ' + pe)
            run_kraken_filters(
                order=order,
                dbs=dbs,
                base_name=pe,
                in_files=fq_path,
                dir_out=dir_fq_filter_data_sample,
                confidence=confidence,
                kraken2=kraken2,
                threads=threads,
                dir_temp=dir_temp,
                linfo=linfo)


def dnld_refseqs_for_taxid(taxid, filter_term, taxonomy,
                           dir_cache_refseqs, db='nuccore', linfo=print):  # noqa
    tax_terms = tuple(reversed(taxonomy.lineage_for_taxid(taxid)['names']))
    for tax_term in tax_terms:
        if tax_term is None:
            tax_term = taxonomy.scientific_name_for_taxid(taxid)
        term = '"RefSeq"[Keyword] AND "{}"[Primary Organism] AND "{}"[filter]'.format(tax_term, filter_term)
        accs = set(accessions_ncbi(esearch(term=term, db=db)))
        if len(accs) > 0:
            plural = 'sequences'
            if len(accs) == 1:
                plural = 'sequence'
            linfo('Found {} RefSeq {} {} for {}.'.format(len(accs), filter_term, plural, tax_term))
            break
        else:
            linfo('No RefSeq {} sequences were found for {}.'.format(filter_term, tax_term))
    cache_path = opj(dir_cache_refseqs, filter_term + '_' + tax_term + '.fasta')
    parsed_fasta_cache = {}
    if ope(cache_path):
        parsed_fasta_cache = read_fasta(cache_path, def_to_first_space=True)
        for acc in parsed_fasta_cache:
            if acc in accs:
                accs.remove(acc)
    if len(accs) > 0:
        parsed_fasta = dnld_seqs_fasta_format(list(accs), db)
        parsed_fasta.update(parsed_fasta_cache)
    else:
        parsed_fasta = parsed_fasta_cache
    write_fasta(parsed_fasta, cache_path)
    return cache_path


def run_bt2_fq(se_fastq_files, pe_fastq_files, dir_fq_filter_data,
               bowtie2, bowtie2_build, threads, dir_temp, filter_dir, dbs,
               fpatt, taxonomy, dir_cache_refseqs, linfo=print):  # noqa

    if len(dbs) > 0:

        if len(se_fastq_files) > 0 or len(pe_fastq_files) > 0:

            if bowtie2 is None:
                linfo('bowtie2 is not available. ' +
                      'Cannot continue. Exiting.')
                exit(0)

            if bowtie2_build is None:
                linfo('bowtie2-build is not available. ' +
                      'Cannot continue. Exiting.')
                exit(0)

    new_se_fastq_files = dict()
    new_pe_fastq_files = dict()

    for db in dbs:

        for se in se_fastq_files:
            dir_fq_bt_data_sample = opj(dir_fq_filter_data, se, filter_dir, db)
            dir_fq_filter_data_sample = opj(dir_fq_filter_data, se, filter_dir)
            in_f = opj(dir_fq_filter_data_sample, se + '.fastq')
            new_se = se + '_' + db
            out_f = opj(dir_fq_bt_data_sample, new_se + '.fastq')
            sam_f = opj(dir_fq_bt_data_sample, new_se + '.sam')
            new_se_fastq_files[new_se] = deepcopy(se_fastq_files[se])
            new_se_fastq_files[new_se]['path'] = None
            new_se_fastq_files[new_se]['cor_path_fq'] = None
            new_se_fastq_files[new_se]['trim_path_fq'] = None
            taxid = new_se_fastq_files[new_se]['tax_id']
            gc = new_se_fastq_files[new_se]['gc_id']
            if db == 'mitochondrion':
                gc = taxonomy.mito_genetic_code_for_taxid(taxid)
                new_se_fastq_files[new_se]['gc_id'] = gc
            elif db == 'chloroplast':
                gc = taxonomy.plastid_genetic_code()
                new_se_fastq_files[new_se]['gc_id'] = gc
            new_se_fastq_files[new_se]['gc_tt'] = TranslationTable(gc)
            new_se_fastq_files[new_se]['filter_path_fq'] = out_f
            if ope(dir_fq_bt_data_sample):
                linfo('Bowtie2 filtered FASTQ file for sample ' + new_se +
                      ' already exists')
            else:
                linfo('Running Bowtie2 in SE mode: ' + new_se)
                make_dir(dir_fq_bt_data_sample)
                db_fasta_path = dnld_refseqs_for_taxid(
                    taxid, db, taxonomy, dir_cache_refseqs,
                    db='nuccore', linfo=linfo)
                bt2_idx_path = db_fasta_path.replace('.fasta', '')
                build_bt2_index(bowtie2_build, [db_fasta_path], bt2_idx_path,
                                threads)

                run_bowtie2_se(bowtie2=bowtie2,
                               input_file=in_f,
                               output_file=out_f,
                               sam_output_file=sam_f,
                               index=bt2_idx_path,
                               threads=threads,
                               dir_temp=dir_temp)

        for pe in pe_fastq_files:
            dir_fq_bt_data_sample = opj(dir_fq_filter_data, pe, filter_dir, db)
            dir_fq_filter_data_sample = opj(dir_fq_filter_data, pe, filter_dir)
            in_fs = [x.replace('@D@', dir_fq_filter_data_sample) for x in fpatt]
            in_fs = [x.replace('@N@', pe) for x in in_fs]
            new_pe = pe + '_' + db
            out_fs = [x.replace('@D@', dir_fq_bt_data_sample) for x in fpatt]
            out_fs = [x.replace('@N@', new_pe) for x in out_fs]
            sam_f = opj(dir_fq_bt_data_sample, new_pe + '.sam')
            new_pe_fastq_files[new_pe] = deepcopy(pe_fastq_files[pe])
            new_pe_fastq_files[new_pe]['path'] = None
            new_pe_fastq_files[new_pe]['cor_path_fq'] = None
            new_pe_fastq_files[new_pe]['trim_path_fq'] = None
            taxid = new_pe_fastq_files[new_pe]['tax_id']
            gc = new_pe_fastq_files[new_pe]['gc_id']
            if db == 'mitochondrion':
                gc = taxonomy.mito_genetic_code_for_taxid(taxid)
                new_pe_fastq_files[new_pe]['gc_id'] = gc
            elif db == 'chloroplast':
                gc = taxonomy.plastid_genetic_code()
                new_pe_fastq_files[new_pe]['gc_id'] = gc
            new_pe_fastq_files[new_pe]['gc_tt'] = TranslationTable(gc)
            new_pe_fastq_files[new_pe]['filter_path_fq'] = out_fs
            if ope(dir_fq_bt_data_sample):
                linfo('Bowtie2 filtered FASTQ files for sample ' + new_pe +
                      ' already exist')
            else:
                linfo('Running Bowtie2 in PE mode: ' + new_pe)
                make_dir(dir_fq_bt_data_sample)
                db_fasta_path = dnld_refseqs_for_taxid(
                    taxid, db, taxonomy, dir_cache_refseqs,
                    db='nuccore', linfo=linfo)
                bt2_idx_path = db_fasta_path.replace('.fasta', '')
                build_bt2_index(bowtie2_build, [db_fasta_path], bt2_idx_path,
                                threads)

                paired_out_pattern = out_fs[0].replace(
                    '_paired_1.fastq', '_paired_%.fastq')

                run_bowtie2_pe(bowtie2=bowtie2,
                               input_files=in_fs,
                               paired_out_pattern=paired_out_pattern,
                               unpaired_out_1=out_fs[2],
                               unpaired_out_2=out_fs[3],
                               sam_output_file=sam_f,
                               index=bt2_idx_path,
                               threads=threads,
                               dir_temp=dir_temp)

    se_fastq_files.update(new_se_fastq_files)
    pe_fastq_files.update(new_pe_fastq_files)


def filtered_fq_to_fa(se_fastq_files, pe_fastq_files, dir_fa_trim_data, seqtk,
                      fpatt, linfo=print): # noqa
    if len(se_fastq_files) > 0 or len(pe_fastq_files) > 0:
        if seqtk is None:
            linfo('seqtk is not available. ' +
                  'Cannot continue. Exiting.')
            exit(0)
    for se in se_fastq_files:
        dir_fa_trim_data_sample = opj(dir_fa_trim_data, se)
        fq_path = se_fastq_files[se]['filter_path_fq']
        out_f = opj(dir_fa_trim_data_sample, se + '.fasta')
        se_fastq_files[se]['filter_path_fa'] = out_f

        if ope(dir_fa_trim_data_sample):
            linfo('Filtered FASTA files for sample ' + se + ' already exist')
        else:
            make_dir(dir_fa_trim_data_sample)
            linfo('Converting FASTQ to FASTA using Seqtk: ' + fq_path)
            seqtk_fq_to_fa(seqtk, fq_path, out_f)

    for pe in pe_fastq_files:
        dir_fa_trim_data_sample = opj(dir_fa_trim_data, pe)
        fq_paths = pe_fastq_files[pe]['filter_path_fq']
        out_fs = [x.replace('@D@', dir_fa_trim_data_sample) for x in fpatt]
        out_fs = [x.replace('@N@', pe) for x in out_fs]
        pe_fastq_files[pe]['filter_path_fa'] = out_fs

        if ope(dir_fa_trim_data_sample):
            linfo('Filtered FASTA files for sample ' + pe + ' already exist')
        else:
            make_dir(dir_fa_trim_data_sample)
            pe_trim_files = zip(fq_paths, out_fs)
            for x in pe_trim_files:
                linfo('Converting FASTQ to FASTA using Seqtk: ' + x[0])
                seqtk_fq_to_fa(seqtk, x[0], x[1])


def makeblastdb_fq(se_fastq_files, pe_fastq_files, dir_blast_fa_trim,
                   makeblastdb, fpatt, linfo=print): # noqa
    if len(se_fastq_files) > 0 or len(pe_fastq_files) > 0:
        if makeblastdb is None:
            linfo('makeblastdb is not available. ' +
                  'Cannot continue. Exiting.')
            exit(0)
    for se in se_fastq_files:
        dir_blast_fa_trim_sample = opj(dir_blast_fa_trim, se)
        fa_path = se_fastq_files[se]['filter_path_fa']
        out_f = opj(dir_blast_fa_trim_sample, se)
        se_fastq_files[se]['blast_db_path'] = out_f

        if ope(dir_blast_fa_trim_sample):
            linfo('BLAST database for sample ' + se + ' already exists')
        else:
            make_dir(dir_blast_fa_trim_sample)
            linfo('Building BLAST database for: ' + fa_path)
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
            linfo('BLAST database for sample ' + pe + ' already exists')
        else:
            make_dir(dir_blast_fa_trim_sample)
            pe_trim_files = zip(fa_paths, out_fs)
            for x in pe_trim_files:
                linfo('Building BLAST database for: ' + x[0])
                make_blast_db(
                    exec_file=makeblastdb,
                    in_file=x[0],
                    out_file=x[1],
                    title=basename(x[1]),
                    dbtype='nucl')


def run_tblastn_on_reads(se_fastq_files, pe_fastq_files, aa_queries_file,
                         tblastn, blast_1_evalue, blast_1_max_hsps,
                         blast_1_qcov_hsp_perc, blast_1_best_hit_overhang,
                         blast_1_best_hit_score_edge, blast_1_max_target_seqs,
                         dir_blast_results_fa_trim, fpatt, ss, threads,
                         seqtk, vsearch, linfo=print): # noqa

    if len(se_fastq_files) > 0 or len(pe_fastq_files) > 0:

        if tblastn is None:
            linfo('tblastn is not available. ' +
                  'Cannot continue. Exiting.')
            exit(0)

        if vsearch is None:
            linfo('vsearch is not available. ' +
                  'Cannot continue. Exiting.')
            exit(0)

        if seqtk is None:
            linfo('seqtk is not available. ' +
                  'Cannot continue. Exiting.')
            exit(0)

    ident = 0.85

    for se in se_fastq_files:
        dir_results = opj(dir_blast_results_fa_trim, se)
        blast_db_path = se_fastq_files[se]['blast_db_path']
        fq_path = se_fastq_files[se]['filter_path_fq']
        out_f = opj(dir_results, se + '__' + ss + '.txt')
        out_f_fastq = out_f.replace('.txt', '.fastq')
        out_f_fasta = out_f.replace('.txt', '.fasta')
        se_fastq_files[se]['blast_results_path' + '__' + ss] = out_f_fasta
        genetic_code = se_fastq_files[se]['gc_id']

        if ope(out_f_fasta):
            linfo('BLAST results for sample ' + se + ' already exists [' + ss + ']')
        else:
            make_dir(dir_results)
            linfo('Running tblastn on: ' + blast_db_path + ' [' + ss + ']')
            run_blast(exec_file=tblastn,
                      task='tblastn',
                      threads=threads,
                      db_path=blast_db_path,
                      queries_file=aa_queries_file,
                      out_file=out_f,
                      evalue=blast_1_evalue,
                      max_hsps=blast_1_max_hsps,
                      qcov_hsp_perc=blast_1_qcov_hsp_perc,
                      best_hit_overhang=blast_1_best_hit_overhang,
                      best_hit_score_edge=blast_1_best_hit_score_edge,
                      max_target_seqs=blast_1_max_target_seqs,
                      db_genetic_code=genetic_code,
                      out_cols=BLST_RES_COLS_1)

            linfo('Extracting unique BLAST hits using Seqtk [' + ss + ']')

            keep_unique_lines_in_file(out_f)

            seqtk_extract_reads(seqtk, fq_path, out_f_fastq, out_f)
            seqtk_fq_to_fa(seqtk, out_f_fastq, out_f_fasta)

            osremove(out_f)
            osremove(out_f_fastq)

            out_f_fasta_temp = out_f_fasta + '_temp'
            copyfile(out_f_fasta, out_f_fasta_temp)
            run_cluster_fast(vsearch, ident, out_f_fasta_temp, out_f_fasta)
            osremove(out_f_fasta_temp)

    for pe in pe_fastq_files:
        dir_results = opj(dir_blast_results_fa_trim, pe)
        blast_db_paths = pe_fastq_files[pe]['blast_db_path']
        fq_paths = pe_fastq_files[pe]['filter_path_fq']
        out_fs = [x.replace('@D@', dir_results) for x in fpatt]
        out_fs = [x.replace('@N@', pe) for x in out_fs]
        out_fs = [x.replace('@Q@', ss) for x in out_fs]
        out_fs_fastq = [x.replace('.txt', '.fastq') for x in out_fs]
        out_fs_fasta = [x.replace('.txt', '.fasta') for x in out_fs]
        out_f_fasta = opj(dir_results, pe + '__' + ss + '.fasta')
        pe_fastq_files[pe]['blast_results_path' + '__' + ss] = out_f_fasta
        genetic_code = pe_fastq_files[pe]['gc_id']

        if ope(out_f_fasta):
            linfo('BLAST results for sample ' + pe + ' already exist [' + ss + ']')
        else:
            make_dir(dir_results)
            pe_trim_files = zip(blast_db_paths, out_fs, fq_paths, out_fs_fastq,
                                out_fs_fasta)
            for x in pe_trim_files:
                linfo('Running tblastn on: ' + x[0] + ' [' + ss + ']')
                run_blast(exec_file=tblastn,
                          task='tblastn',
                          threads=threads,
                          db_path=x[0],
                          queries_file=aa_queries_file,
                          out_file=x[1],
                          evalue=blast_1_evalue,
                          max_hsps=blast_1_max_hsps,
                          qcov_hsp_perc=blast_1_qcov_hsp_perc,
                          best_hit_overhang=blast_1_best_hit_overhang,
                          best_hit_score_edge=blast_1_best_hit_score_edge,
                          max_target_seqs=blast_1_max_target_seqs,
                          db_genetic_code=genetic_code,
                          out_cols=BLST_RES_COLS_1)

                linfo('Extracting unique BLAST hits using Seqtk [' + ss + ']')

                keep_unique_lines_in_file(x[1])

                seqtk_extract_reads(seqtk, x[2], x[3], x[1])
                seqtk_fq_to_fa(seqtk, x[3], x[4])

                osremove(x[1])
                osremove(x[3])

            combine_text_files(out_fs_fasta, out_f_fasta)

            out_f_fasta_temp = out_f_fasta + '_temp'
            copyfile(out_f_fasta, out_f_fasta_temp)
            run_cluster_fast(vsearch, ident, out_f_fasta_temp, out_f_fasta)
            osremove(out_f_fasta_temp)

            for x in out_fs_fasta:
                osremove(x)


def run_vsearch_on_reads(se_fastq_files, pe_fastq_files, vsearch,
                         dir_vsearch_results_fa_trim, fpatt, ss, seqtk,
                         linfo=print): # noqa

    if len(se_fastq_files) > 0 or len(pe_fastq_files) > 0:
        if vsearch is None:
            linfo('vsearch is not available. ' +
                  'Cannot continue. Exiting.')
            exit(0)
        if seqtk is None:
            linfo('seqtk is not available. ' +
                  'Cannot continue. Exiting.')
            exit(0)

    ident = 0.85

    for se in se_fastq_files:
        dir_results = opj(dir_vsearch_results_fa_trim, se)
        min_acc_len = se_fastq_files[se]['min_acc_len']
        blast_results_fa_path = se_fastq_files[se]['blast_results_path' + '__' + ss]
        fq_path = se_fastq_files[se]['filter_path_fq']
        out_f = opj(dir_results, se + '__' + ss + '.txt')
        out_f_fastq = out_f.replace('.txt', '.fastq')
        se_fastq_files[se]['vsearch_results_path' + '__' + ss] = out_f_fastq

        if ope(out_f_fastq):
            linfo('Vsearch results for sample ' + se + ' already exists [' + ss + ']')
        else:
            make_dir(dir_results)
            linfo('Running vsearch on: ' + fq_path + ' [' + ss + ']')
            run_vsearch(vsearch,
                        ident=ident,
                        q_file=blast_results_fa_path,
                        db_file=fq_path,
                        out_file=out_f,
                        minlen=min_acc_len)

            linfo('Extracting unique vsearch hits using Seqtk [' + ss + ']')
            keep_unique_lines_in_file(out_f)
            seqtk_extract_reads(seqtk, fq_path, out_f_fastq, out_f)
            osremove(out_f)

    for pe in pe_fastq_files:
        dir_results = opj(dir_vsearch_results_fa_trim, pe)
        min_acc_len = pe_fastq_files[pe]['min_acc_len']
        blast_results_fa_path = pe_fastq_files[pe]['blast_results_path' + '__' + ss]
        fq_paths = pe_fastq_files[pe]['filter_path_fq']
        out_fs = [x.replace('@D@', dir_results) for x in fpatt]
        out_fs = [x.replace('@N@', pe) for x in out_fs]
        out_fs = [x.replace('@Q@', ss) for x in out_fs]
        out_fs_fastq = [x.replace('.txt', '.fastq') for x in out_fs]
        pe_fastq_files[pe]['vsearch_results_path' + '__' + ss] = out_fs_fastq

        if ope(out_fs_fastq[0]) and ope(out_fs_fastq[1]) and \
           ope(out_fs_fastq[2]) and ope(out_fs_fastq[3]):
            linfo('Vsearch results for sample ' + pe + ' already exist [' + ss + ']')
        else:
            make_dir(dir_results)
            pe_trim_files = zip(fq_paths, out_fs, out_fs_fastq)
            for x in pe_trim_files:
                linfo('Running vsearch on: ' + x[0] + ' [' + ss + ']')
                run_vsearch(vsearch,
                            ident=ident,
                            q_file=blast_results_fa_path,
                            db_file=x[0],
                            out_file=x[1],
                            minlen=min_acc_len)

            linfo('Extracting unique vsearch hits from paired files '
                  'using Seqtk [' + ss + ']')

            p1txt = out_fs[0]
            p2txt = out_fs[1]

            p1fq = fq_paths[0]
            p2fq = fq_paths[1]

            p1fq_out = out_fs_fastq[0]
            p2fq_out = out_fs_fastq[1]

            p12txt_temp = opj(dir_results, pe + '__' + ss + '_paired.txt')

            combine_text_files([p1txt, p2txt], p12txt_temp)
            keep_unique_lines_in_file(p12txt_temp)

            seqtk_extract_reads(seqtk, p1fq, p1fq_out, p12txt_temp)
            seqtk_extract_reads(seqtk, p2fq, p2fq_out, p12txt_temp)

            osremove(p1txt)
            osremove(p2txt)
            osremove(p12txt_temp)

            linfo('Extracting unique vsearch hits from unpaired files '
                  'using Seqtk [' + ss + ']')

            u1txt = out_fs[2]
            u2txt = out_fs[3]

            u1fq = fq_paths[2]
            u2fq = fq_paths[3]

            u1fq_out = out_fs_fastq[2]
            u2fq_out = out_fs_fastq[3]

            keep_unique_lines_in_file(u1txt)
            keep_unique_lines_in_file(u2txt)

            seqtk_extract_reads(seqtk, u1fq, u1fq_out, u1txt)
            seqtk_extract_reads(seqtk, u2fq, u2fq_out, u2txt)

            osremove(u1txt)
            osremove(u2txt)
