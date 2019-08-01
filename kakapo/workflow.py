#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""kakapo workflow"""

import json
import pickle
import re

from os import remove as osremove
from os.path import basename
from os.path import commonprefix
from os.path import exists as ope
from os.path import join as opj
from os.path import splitext
from shutil import copyfile
from sys import exit
from time import sleep

from kakapo.bioio import filter_fasta_text_by_length
from kakapo.bioio import parse_fasta_text
from kakapo.bioio import read_fasta_file, read_fasta_file_dict
from kakapo.bioio import standardize_fasta_text
from kakapo.bioio import trim_desc_to_first_space_in_fasta_text
from kakapo.bioio import write_fasta_file
from kakapo.blast import BLST_RES_COLS_1, BLST_RES_COLS_2
from kakapo.blast import collate_blast_results
from kakapo.blast import make_blast_db, run_blast
from kakapo.blast import parse_blast_results_file
from kakapo.config import PICKLE_PROTOCOL, CONSRED, CONYELL, CONSDFL
from kakapo.ebi_domain_search import pfam_entry
from kakapo.ebi_domain_search import pfam_seqs
from kakapo.ebi_domain_search import prot_ids_for_tax_ids
from kakapo.ebi_iprscan5 import job_runner
from kakapo.ebi_iprscan5 import result_json
from kakapo.ebi_proteins import fasta_by_accession_list
from kakapo.entrez import dnld_seqs as dnld_ncbi_seqs
from kakapo.entrez import sra_run_info
from kakapo.entrez import summary as entrez_summary
from kakapo.gff3 import gff_from_kakapo_ips5_json_file
from kakapo.helpers import combine_text_files
from kakapo.helpers import keep_unique_lines_in_file
from kakapo.helpers import make_dir
from kakapo.orf import find_orf_for_blast_hit
from kakapo.seq import reverse_complement, translate
from kakapo.seqtk import seqtk_fq_to_fa, seqtk_extract_reads
from kakapo.shell import call
from kakapo.spades import run_spades_se, run_spades_pe
from kakapo.trimmomatic import trimmomatic_se, trimmomatic_pe
from kakapo.vsearch import run_cluster_fast, run_vsearch

def prepare_output_directories(dir_out, prj_name):  # noqa

    # ToDo: Lock cache files in case of parallel execution -------------------

    dir_temp = opj(dir_out, '00-temp')
    make_dir(dir_temp)

    dir_cache = opj(dir_out, '00-cache')
    make_dir(dir_cache)

    dir_cache_pfam_acc = opj(dir_cache, 'pfam-uniprot-accessions')
    make_dir(dir_cache_pfam_acc)

    dir_cache_fq_minlen = opj(dir_cache, 'min-acceptable-read-lengths')
    make_dir(dir_cache_fq_minlen)

    dir_cache_prj = opj(dir_cache, 'projects', prj_name)
    make_dir(dir_cache_prj)

    dir_prj = opj(dir_out, '01-project-specific', prj_name)
    make_dir(dir_prj)

    dir_prj_queries = opj(dir_prj, '01-queries')
    make_dir(dir_prj_queries)

    dir_prj_blast_results_fa_trim = opj(dir_prj, '02-trimmed-fa-blast-results')
    make_dir(dir_prj_blast_results_fa_trim)

    dir_prj_vsearch_results_fa_trim = opj(dir_prj,
                                          '03-trimmed-fa-vsearch-results')
    make_dir(dir_prj_vsearch_results_fa_trim)

    dir_prj_spades_assemblies = opj(dir_prj, '04-spades-assemblies')
    make_dir(dir_prj_spades_assemblies)

    dir_prj_blast_assmbl = opj(dir_prj, '05-assemblies-blast-db-data')
    make_dir(dir_prj_blast_assmbl)

    dir_prj_assmbl_blast_results = opj(dir_prj, '06-assemblies-blast-results')
    make_dir(dir_prj_assmbl_blast_results)

    dir_prj_transcripts = opj(dir_prj, '07-transcripts')
    make_dir(dir_prj_transcripts)

    dir_prj_ips = dir_prj_transcripts

    dir_prj_transcripts_combined = opj(dir_prj, '08-transcripts-combined')
    make_dir(dir_prj_transcripts_combined)

    dir_global = opj(dir_out, '02-global')
    make_dir(dir_global)

    dir_fq_data = opj(dir_global, '01-sra-fq-data')
    make_dir(dir_fq_data)

    dir_fq_trim_data = opj(dir_global, '02-trimmed-fq-data')
    make_dir(dir_fq_trim_data)

    dir_fa_trim_data = opj(dir_global, '03-trimmed-fa-data')
    make_dir(dir_fa_trim_data)

    dir_blast_fa_trim = opj(dir_global, '04-trimmed-fa-blast-db-data')
    make_dir(dir_blast_fa_trim)

    ret_dict = {'dir_temp': dir_temp,
                'dir_cache': dir_cache,
                'dir_cache_pfam_acc': dir_cache_pfam_acc,
                'dir_cache_fq_minlen': dir_cache_fq_minlen,
                'dir_cache_prj': dir_cache_prj,
                'dir_prj': dir_prj,
                'dir_prj_queries': dir_prj_queries,
                'dir_fq_data': dir_fq_data,
                'dir_fq_trim_data': dir_fq_trim_data,
                'dir_fa_trim_data': dir_fa_trim_data,
                'dir_blast_fa_trim': dir_blast_fa_trim,
                'dir_prj_blast_results_fa_trim': dir_prj_blast_results_fa_trim,
                'dir_prj_vsearch_results_fa_trim':
                    dir_prj_vsearch_results_fa_trim,
                'dir_prj_spades_assemblies': dir_prj_spades_assemblies,
                'dir_prj_blast_assmbl': dir_prj_blast_assmbl,
                'dir_prj_assmbl_blast_results': dir_prj_assmbl_blast_results,
                'dir_prj_transcripts': dir_prj_transcripts,
                'dir_prj_ips': dir_prj_ips,
                'dir_prj_transcripts_combined': dir_prj_transcripts_combined}

    return ret_dict

def descending_tax_ids(tax_ids_user, taxonomy):  # noqa
    if len(tax_ids_user) > 0:
        print('Resolving descending nodes for:\n')
    tax_ids = []
    for tx in tax_ids_user:
        tx_name = taxonomy.scientific_name_for_taxid(taxid=tx)
        print('\t' + str(tx) + ':\t' + tx_name)
        tax_ids_for_tx = taxonomy.all_descending_taxids(taxid=tx)
        tax_ids = tax_ids + tax_ids_for_tx
    tax_ids = list(set(tax_ids))
    # ToDo: Fix Taxonomy class, so it returns integers by default
    tax_ids = [int(x) for x in tax_ids]

    return tax_ids


def pfam_uniprot_accessions(pfam_acc, tax_ids, dir_cache_pfam_acc):  # noqa
    if len(pfam_acc) > 0:
        print('\nDownloading UniProt accessions for Pfam families:\n')
    pfam_seqs_list = []
    for pa in pfam_acc:
        pfam_id = pfam_entry(pa)[0]['id']
        print('\t' + pa + ': ' + pfam_id)
        __ = opj(dir_cache_pfam_acc, pa)
        if ope(__):
            with open(__, 'rb') as f:
                acc = pickle.load(f)
            pfam_seqs_list = pfam_seqs_list + acc
        else:
            acc = pfam_seqs(query=pa)
            pfam_seqs_list = pfam_seqs_list + acc
            with open(__, 'wb') as f:
                pickle.dump(acc, f, protocol=PICKLE_PROTOCOL)

    pfam_uniprot_acc = prot_ids_for_tax_ids(pfam_seqs_list, tax_ids)
    print()
    return pfam_uniprot_acc


def dnld_pfam_uniprot_seqs(uniprot_acc, aa_uniprot_file, dir_cache_prj):  # noqa
    if len(uniprot_acc) != 0:
        __ = opj(dir_cache_prj, 'aa_uniprot_acc_cache')
        prev_uniprot_acc = []
        if ope(__):
            with open(__, 'rb') as f:
                prev_uniprot_acc = pickle.load(f)

        with open(__, 'wb') as f:
            pickle.dump(uniprot_acc, f, protocol=PICKLE_PROTOCOL)

        if (set(uniprot_acc) != set(prev_uniprot_acc)) or \
           (not ope(aa_uniprot_file)):

            print('Downloading Pfam protein sequences from UniProt.\n')
            __ = fasta_by_accession_list(uniprot_acc)
            __ = standardize_fasta_text(__)

            with open(aa_uniprot_file, 'w') as f:
                f.write(__)
    else:
        if ope(aa_uniprot_file):
            osremove(aa_uniprot_file)


def user_protein_accessions(prot_acc_user):  # noqa
    if len(prot_acc_user) > 0:
        print('Reading user provided protein accessions:\n')
        pa_info = entrez_summary(prot_acc_user, 'protein')
        prot_acc = []
        for pa in pa_info:
            title = pa['Title']
            title = title[0].upper() + title[1:]
            acc = pa['AccessionVersion']
            prot_acc.append(acc)

            if len(title) > 60:
                title = title[0:57] + '...'

            print('\t' + acc + ': ' + title)
        print()

        return prot_acc

    else:
        return prot_acc_user


def dnld_prot_seqs(prot_acc_user, aa_prot_ncbi_file, dir_cache_prj):  # noqa
    if len(prot_acc_user) != 0:
        __ = opj(dir_cache_prj, 'aa_prot_ncbi_acc_cache')
        prev_prot_acc_user = []
        if ope(__):
            with open(__, 'rb') as f:
                prev_prot_acc_user = pickle.load(f)

        with open(__, 'wb') as f:
            pickle.dump(prot_acc_user, f, protocol=PICKLE_PROTOCOL)

        if (set(prot_acc_user) != set(prev_prot_acc_user)) or \
           (not ope(aa_prot_ncbi_file)):

            print('Downloading protein sequences from NCBI.\n')
            __ = dnld_ncbi_seqs(prot_acc_user, 'protein')

            write_fasta_file(__, aa_prot_ncbi_file)

    else:
        if ope(aa_prot_ncbi_file):
            osremove(aa_prot_ncbi_file)


def user_aa_fasta(user_queries, aa_prot_user_file):  # noqa
    __ = ''
    if len(user_queries) > 0:
        print('Reading user provided AA sequences:')
        for ap in user_queries:
            print('\t\t' + ap)
            with open(ap, 'r') as f:
                __ = __ + f.read()
        print()
    if __ != '':
        with open(aa_prot_user_file, 'w') as f:
            f.write(standardize_fasta_text(__))


def combine_aa_fasta(fasta_files, aa_queries_file):  # noqa
    print('Combining all AA query sequences.\n')
    __ = ''
    for fasta_file in fasta_files:
        if ope(fasta_file):
            with open(fasta_file, 'r') as f:
                __ = __ + f.read()

    if __ != '':
        with open(aa_queries_file, 'w') as f:
            f.write(__)
    else:
        print('No queries were provided. Exiting.')
        exit(0)


def filter_queries(aa_queries_file, min_query_length, max_query_length): # noqa
    __ = ''
    with open(aa_queries_file, 'r') as f:
        __ = f.read()

    print('Filtering AA query sequences:\n')
    print('\tmin_query_length: ' + str(min_query_length))
    print('\tmax_query_length: ' + str(max_query_length))

    __ = filter_fasta_text_by_length(__, min_query_length, max_query_length)

    with open(aa_queries_file, 'w') as f:
        f.write(__)


def dnld_sra_info(sras, dir_cache_prj):  # noqa

    sra_runs_info = {}
    sras_acceptable = []

    if len(sras) > 0:
        print('\nDownloading SRA run information:\n')
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
            sra_spots_with_mates = int(info['spots_with_mates'])

            sample_base_name = (sra_species.replace(' ', '_') + '_' +
                                sra_taxid + '_' + sra)

            sra_runs_info[sra]['KakapoSampleBaseName'] = sample_base_name

            if sra_lib_source != 'transcriptomic':
                sra_info_str = (
                    '{sra}: the SRA library source type "{ltype}" '
                    'is not supported.').format(
                    sra=sra, ltype=sra_lib_source)

            elif sra_seq_platform != 'Illumina':
                sra_info_str = (
                    '{sra}: the SRA library sequencing platform "{plat}" '
                    'is not supported.').format(
                    sra=sra, plat=sra_seq_platform)

            else:
                sra_info_str = ('\tSRA run {sra}\n\t{source} '
                                '{strategy} {layout}-end library.\n'
                                '\tSourced from {species} '
                                '(Tax ID: {txid}).\n'
                                '\tSequenced using {platform} platform on '
                                '{model}.\n').format(
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

                if sra_lib_layout == 'paired':
                    if sra_spots_with_mates == 0:
                        sra_runs_info[sra]['KakapoLibraryLayout'] = 'SINGLE'
                        sra_info_str = (
                            sra_info_str + CONSRED + '\t>>> ' + CONYELL +
                            sra + CONSDFL + ' is listed as containing '
                            'paired-end reads, but only a single set of reads '
                            'is available. Treating as single-ended.\n')

                sras_acceptable.append(sra)

            print(sra_info_str)

    with open(__, 'wb') as f:
        pickle.dump(sra_runs_info, f, protocol=PICKLE_PROTOCOL)

    return sra_runs_info, sras_acceptable


def dnld_sra_fastq_files(sras, sra_runs_info, dir_fq_data, fasterq_dump,
                         threads, dir_temp): # noqa

    se_fastq_files = {}
    pe_fastq_files = {}

    for sra in sras:
        sra_run_info = sra_runs_info[sra]
        sra_lib_layout = sra_run_info['LibraryLayout'].lower()
        sra_lib_layout_k = sra_run_info['KakapoLibraryLayout'].lower()
        sample_base_name = sra_run_info['KakapoSampleBaseName']
        avg_len = int(sra_run_info['avgLength'])

        sra_dnld_needed = False

        if sra_lib_layout == 'single' or sra_lib_layout_k == 'single':
            se_file = opj(dir_fq_data, sra + '.fastq')
            se_fastq_files[sample_base_name] = {'path': se_file}
            se_fastq_files[sample_base_name]['src'] = 'sra'
            se_fastq_files[sample_base_name]['avg_len'] = avg_len
            if not ope(se_file):
                sra_dnld_needed = True

        elif sra_lib_layout == 'paired':
            pe_file_1 = opj(dir_fq_data, sra + '_1.fastq')
            pe_file_2 = opj(dir_fq_data, sra + '_2.fastq')
            pe_fastq_files[sample_base_name] = {'path': [pe_file_1, pe_file_2]}
            pe_fastq_files[sample_base_name]['src'] = 'sra'
            pe_fastq_files[sample_base_name]['avg_len'] = avg_len // 2
            if not ope(pe_file_1) or not ope(pe_file_2):
                sra_dnld_needed = True

        if sra_dnld_needed:
            print('\nDownloading FASTQ reads for ' + sample_base_name + '\n')
            cmd = [fasterq_dump,
                   '--threads', str(threads),
                   '--split-files',
                   '--outdir', dir_fq_data,
                   '--temp', dir_temp, sra]
            call(cmd)

        else:
            print('\tFASTQ reads for the SRA run ' + sample_base_name +
                  ' are available locally.')
    print()

    return se_fastq_files, pe_fastq_files, sra_runs_info


def user_fastq_files(fq_se, fq_pe): # noqa

    if len(fq_se) > 0 or len(fq_pe) > 0:
        print('\nPreparing user provided FASTQ files:\n')

    se_fastq_files = {}
    pe_fastq_files = {}

    for se in fq_se:
        sample_base_name = splitext(basename(se))[0]
        se_fastq_files[sample_base_name] = {'path': se}
        se_fastq_files[sample_base_name]['src'] = 'usr'
        se_fastq_files[sample_base_name]['avg_len'] = None
        print('\t' + sample_base_name + ':\n\t\t' + se)
        print()

    for pe in fq_pe:
        sample_base_name = basename(commonprefix(pe))
        sample_base_name = sample_base_name.rstrip('_- R')
        pe_fastq_files[sample_base_name] = {'path': pe}
        pe_fastq_files[sample_base_name]['src'] = 'usr'
        pe_fastq_files[sample_base_name]['avg_len'] = None
        print('\t' + sample_base_name + ':\n\t\t' + pe[0] + '\n\t\t' + pe[1])
        print()

    return se_fastq_files, pe_fastq_files


def min_accept_read_len(se_fastq_files, pe_fastq_files, dir_temp,
                        dir_cache_fq_minlen, vsearch): # noqa

    # lowest allowable
    low = 30

    if len(se_fastq_files) > 0 or len(pe_fastq_files) > 0:
        print('\nCalculating minimum acceptable read length:\n')
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
            ml = max(avg_len // 2, low)
            se_fastq_files[se]['min_acc_len'] = ml
            print('\t' + str(ml) + ' nt: ' + se)
            continue

        fq_path = se_fastq_files[se]['path']
        stats_file = opj(dir_temp, se + '_stats.txt')
        queue.append([se, fq_path, stats_file, 'se'])

    for pe in pe_fastq_files:
        src = pe_fastq_files[pe]['src']
        avg_len = pe_fastq_files[pe]['avg_len']
        if src == 'sra':
            ml = max(avg_len // 2, low)
            pe_fastq_files[pe]['min_acc_len'] = ml
            print('\t' + str(ml) + ' nt: ' + pe)
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
                ml = max(int(ml[0]) // 2, low)
            else:
                ml = None

            pickled[x[0]] = ml

        if ml is not None:
            print('\t' + str(ml) + ' nt: ' + x[0])
        else:
            print('\t' + ' ?' + ' nt: ' + x[0])
            ml = low

        if x[3] == 'se':
            se_fastq_files[x[0]]['min_acc_len'] = ml

        elif x[3] == 'pe':
            pe_fastq_files[x[0]]['min_acc_len'] = ml

        with open(__, 'wb') as f:
            pickle.dump(pickled, f, protocol=PICKLE_PROTOCOL)

    print()


def run_trimmomatic(se_fastq_files, pe_fastq_files, dir_fq_trim_data,
                    trimmomatic, adapters, fpatt, threads): # noqa

    for se in se_fastq_files:
        dir_fq_trim_data_sample = opj(dir_fq_trim_data, se)
        fq_path = se_fastq_files[se]['path']
        min_acc_len = se_fastq_files[se]['min_acc_len']
        stats_f = opj(dir_fq_trim_data_sample, se + '.txt')
        out_f = opj(dir_fq_trim_data_sample, se + '.fastq')
        se_fastq_files[se]['trim_path_fq'] = out_f

        if ope(dir_fq_trim_data_sample):
            print('Trimmed FASTQ files for sample ' + se + ' already exist.')
        else:
            make_dir(dir_fq_trim_data_sample)
            print('Running Trimmomatic in SE mode: ' + se)
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
        fq_path_1 = pe_fastq_files[pe]['path'][0]
        fq_path_2 = pe_fastq_files[pe]['path'][1]
        min_acc_len = pe_fastq_files[pe]['min_acc_len']
        stats_f = opj(dir_fq_trim_data_sample, pe + '.txt')
        out_fs = [x.replace('@D@', dir_fq_trim_data_sample) for x in fpatt]
        out_fs = [x.replace('@N@', pe) for x in out_fs]
        pe_fastq_files[pe]['trim_path_fq'] = out_fs

        if ope(dir_fq_trim_data_sample):
            print('Trimmed FASTQ files for sample ' + pe + ' already exist.')
        else:
            make_dir(dir_fq_trim_data_sample)
            print('Running Trimmomatic in PE mode: ' + pe)
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

    if len(se_fastq_files) > 0 or len(pe_fastq_files) > 0:
        print()


def trimmed_fq_to_fa(se_fastq_files, pe_fastq_files, dir_fa_trim_data, seqtk,
                     fpatt): # noqa

    for se in se_fastq_files:
        dir_fa_trim_data_sample = opj(dir_fa_trim_data, se)
        fq_path = se_fastq_files[se]['trim_path_fq']
        out_f = opj(dir_fa_trim_data_sample, se + '.fasta')
        se_fastq_files[se]['trim_path_fa'] = out_f

        if ope(dir_fa_trim_data_sample):
            print('Trimmed FASTA files for sample ' + se + ' already exist.')
        else:
            make_dir(dir_fa_trim_data_sample)
            print('Converting FASTQ to FASTA using Seqtk: ' + fq_path)
            seqtk_fq_to_fa(seqtk, fq_path, out_f)

    for pe in pe_fastq_files:
        dir_fa_trim_data_sample = opj(dir_fa_trim_data, pe)
        fq_paths = pe_fastq_files[pe]['trim_path_fq']
        out_fs = [x.replace('@D@', dir_fa_trim_data_sample) for x in fpatt]
        out_fs = [x.replace('@N@', pe) for x in out_fs]
        pe_fastq_files[pe]['trim_path_fa'] = out_fs

        if ope(dir_fa_trim_data_sample):
            print('Trimmed FASTA files for sample ' + pe + ' already exist.')
        else:
            make_dir(dir_fa_trim_data_sample)
            pe_trim_files = zip(fq_paths, out_fs)
            for x in pe_trim_files:
                print('Converting FASTQ to FASTA using Seqtk: ' + x[0])
                seqtk_fq_to_fa(seqtk, x[0], x[1])


def makeblastdb_fq(se_fastq_files, pe_fastq_files, dir_blast_fa_trim,
                   makeblastdb, fpatt): # noqa

    if len(se_fastq_files) > 0 or len(pe_fastq_files) > 0:
        print()

    for se in se_fastq_files:
        dir_blast_fa_trim_sample = opj(dir_blast_fa_trim, se)
        fa_path = se_fastq_files[se]['trim_path_fa']
        out_f = opj(dir_blast_fa_trim_sample, se)
        se_fastq_files[se]['blast_db_path'] = out_f

        if ope(dir_blast_fa_trim_sample):
            print('BLAST database for sample ' + se + ' already exists.')
        else:
            make_dir(dir_blast_fa_trim_sample)
            print('Building BLAST database for: ' + fa_path)
            make_blast_db(
                exec_file=makeblastdb,
                in_file=fa_path,
                out_file=out_f,
                title=se,
                dbtype='nucl')

    for pe in pe_fastq_files:
        dir_blast_fa_trim_sample = opj(dir_blast_fa_trim, pe)
        fa_paths = pe_fastq_files[pe]['trim_path_fa']
        out_fs = [x.replace('@D@', dir_blast_fa_trim_sample) for x in fpatt]
        out_fs = [x.replace('@N@', pe) for x in out_fs]
        pe_fastq_files[pe]['blast_db_path'] = out_fs

        if ope(dir_blast_fa_trim_sample):
            print('BLAST database for sample ' + pe + ' already exists.')
        else:
            make_dir(dir_blast_fa_trim_sample)
            pe_trim_files = zip(fa_paths, out_fs)
            for x in pe_trim_files:
                print('Building BLAST database for: ' + x[0])
                make_blast_db(
                    exec_file=makeblastdb,
                    in_file=x[0],
                    out_file=x[1],
                    title=basename(x[1]),
                    dbtype='nucl')


def run_tblastn_on_reads(se_fastq_files, pe_fastq_files, aa_queries_file,
                         tblastn, blast_1_evalue, blast_1_max_target_seqs,
                         blast_1_culling_limit, blast_1_qcov_hsp_perc,
                         dir_blast_results_fa_trim, fpatt, threads,
                         genetic_code, seqtk, vsearch): # noqa

    if len(se_fastq_files) > 0 or len(pe_fastq_files) > 0:
        print()

    ident = 0.85

    for se in se_fastq_files:
        dir_results = opj(dir_blast_results_fa_trim, se)
        blast_db_path = se_fastq_files[se]['blast_db_path']
        fq_path = se_fastq_files[se]['trim_path_fq']
        out_f = opj(dir_results, se + '.txt')
        out_f_fastq = out_f.replace('.txt', '.fastq')
        out_f_fasta = out_f.replace('.txt', '.fasta')
        se_fastq_files[se]['blast_results_path'] = out_f_fasta

        if ope(dir_results):
            print('BLAST results for sample ' + se + ' already exists.')
        else:
            make_dir(dir_results)
            print('Running tblastn on: ' + blast_db_path)
            run_blast(exec_file=tblastn,
                      task='tblastn',
                      threads=threads,
                      db_path=blast_db_path,
                      queries_file=aa_queries_file,
                      out_file=out_f,
                      evalue=blast_1_evalue,
                      qcov_hsp_perc=blast_1_qcov_hsp_perc,
                      culling_limit=blast_1_culling_limit,
                      max_target_seqs=blast_1_max_target_seqs,
                      db_genetic_code=genetic_code,
                      out_cols=BLST_RES_COLS_1)

            print('\tExtracting unique BLAST hits using Seqtk.\n')

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
        fq_paths = pe_fastq_files[pe]['trim_path_fq']
        out_fs = [x.replace('@D@', dir_results) for x in fpatt]
        out_fs = [x.replace('@N@', pe) for x in out_fs]
        out_fs_fastq = [x.replace('.txt', '.fastq') for x in out_fs]
        out_fs_fasta = [x.replace('.txt', '.fasta') for x in out_fs]
        out_f_fasta = opj(dir_results, pe + '.fasta')
        pe_fastq_files[pe]['blast_results_path'] = out_f_fasta

        if ope(dir_results):
            print('BLAST results for sample ' + pe + ' already exist.')
        else:
            make_dir(dir_results)
            pe_trim_files = zip(blast_db_paths, out_fs, fq_paths, out_fs_fastq,
                                out_fs_fasta)
            for x in pe_trim_files:
                print('Running tblastn on: ' + x[0])
                run_blast(exec_file=tblastn,
                          task='tblastn',
                          threads=threads,
                          db_path=x[0],
                          queries_file=aa_queries_file,
                          out_file=x[1],
                          evalue=blast_1_evalue,
                          qcov_hsp_perc=blast_1_qcov_hsp_perc,
                          culling_limit=blast_1_culling_limit,
                          max_target_seqs=blast_1_max_target_seqs,
                          db_genetic_code=genetic_code,
                          out_cols=BLST_RES_COLS_1)

                print('\tExtracting unique BLAST hits using Seqtk.\n')

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
                         dir_vsearch_results_fa_trim, fpatt, seqtk): # noqa

    if len(se_fastq_files) > 0 or len(pe_fastq_files) > 0:
        print()

    ident = 0.85

    for se in se_fastq_files:
        dir_results = opj(dir_vsearch_results_fa_trim, se)
        min_acc_len = se_fastq_files[se]['min_acc_len']
        blast_results_fa_path = se_fastq_files[se]['blast_results_path']
        fq_path = se_fastq_files[se]['trim_path_fq']
        out_f = opj(dir_results, se + '.txt')
        out_f_fastq = out_f.replace('.txt', '.fastq')
        se_fastq_files[se]['vsearch_results_path'] = out_f_fastq

        if ope(dir_results):
            print('Vsearch results for sample ' + se + ' already exists.')
        else:
            make_dir(dir_results)
            print('Running vsearch on: ' + fq_path)
            run_vsearch(vsearch,
                        ident=ident,
                        q_file=blast_results_fa_path,
                        db_file=fq_path,
                        out_file=out_f,
                        minlen=min_acc_len)

            print('\tExtracting unique vsearch hits using Seqtk.\n')
            keep_unique_lines_in_file(out_f)
            seqtk_extract_reads(seqtk, fq_path, out_f_fastq, out_f)
            osremove(out_f)

    for pe in pe_fastq_files:
        dir_results = opj(dir_vsearch_results_fa_trim, pe)
        min_acc_len = pe_fastq_files[pe]['min_acc_len']
        blast_results_fa_path = pe_fastq_files[pe]['blast_results_path']
        fq_paths = pe_fastq_files[pe]['trim_path_fq']
        out_fs = [x.replace('@D@', dir_results) for x in fpatt]
        out_fs = [x.replace('@N@', pe) for x in out_fs]
        out_fs_fastq = [x.replace('.txt', '.fastq') for x in out_fs]
        pe_fastq_files[pe]['vsearch_results_path'] = out_fs_fastq

        if ope(dir_results):
            print('Vsearch results for sample ' + pe + ' already exist.')
        else:
            make_dir(dir_results)
            pe_trim_files = zip(fq_paths, out_fs, out_fs_fastq)
            for x in pe_trim_files:
                print('Running vsearch on: ' + x[0])
                run_vsearch(vsearch,
                            ident=ident,
                            q_file=blast_results_fa_path,
                            db_file=x[0],
                            out_file=x[1],
                            minlen=min_acc_len)

            print('\tExtracting unique vsearch hits from paired files '
                  'using Seqtk.')

            p1txt = out_fs[0]
            p2txt = out_fs[1]

            p1fq = fq_paths[0]
            p2fq = fq_paths[1]

            p1fq_out = out_fs_fastq[0]
            p2fq_out = out_fs_fastq[1]

            p12txt_temp = opj(dir_results, pe + '_paired.txt')

            combine_text_files([p1txt, p2txt], p12txt_temp)
            keep_unique_lines_in_file(p12txt_temp)

            seqtk_extract_reads(seqtk, p1fq, p1fq_out, p12txt_temp)
            seqtk_extract_reads(seqtk, p2fq, p2fq_out, p12txt_temp)

            osremove(p1txt)
            osremove(p2txt)
            osremove(p12txt_temp)

            print('\tExtracting unique vsearch hits from unpaired files '
                  'using Seqtk.\n')

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


def run_spades(se_fastq_files, pe_fastq_files, dir_spades_assemblies,
               spades, threads, ram):  # noqa

    if len(se_fastq_files) > 0 or len(pe_fastq_files) > 0:
        print()

    for se in se_fastq_files:
        dir_results = opj(dir_spades_assemblies, se)
        fq_path = se_fastq_files[se]['vsearch_results_path']
        se_fastq_files[se]['spades_assembly'] = None

        if ope(dir_results):
            print('SPAdes results for sample ' + se + ' already exist.')
        else:
            make_dir(dir_results)
            print('Running SPAdes on: ' + se)
            run_spades_se(spades,
                          out_dir=dir_results,
                          input_file=fq_path,
                          threads=threads,
                          memory=ram,
                          rna=True)

        assmbl_path = opj(dir_results, 'transcripts.fasta')
        if ope(assmbl_path):
            count = len(read_fasta_file(assmbl_path))
            tr_str = ' transcripts.'
            if count == 1:
                tr_str = ' transcript.'
            print('\n\tSPAdes produced ' + str(count) + tr_str)
            se_fastq_files[se]['spades_assembly'] = assmbl_path
        else:
            print('\n\tSPAdes produced no transcripts.')
        print()

    for pe in pe_fastq_files:
        dir_results = opj(dir_spades_assemblies, pe)
        fq_paths = pe_fastq_files[pe]['vsearch_results_path']
        pe_fastq_files[pe]['spades_assembly'] = None

        if ope(dir_results):
            print('SPAdes results for sample ' + pe + ' already exist.')
        else:
            make_dir(dir_results)
            print('Running SPAdes on: ' + pe)
            run_spades_pe(spades,
                          out_dir=dir_results,
                          input_files=fq_paths,
                          threads=threads,
                          memory=ram,
                          rna=True)

        assmbl_path = opj(dir_results, 'transcripts.fasta')
        if ope(assmbl_path):
            count = len(read_fasta_file(assmbl_path))
            tr_str = ' transcripts.'
            if count == 1:
                tr_str = ' transcript.'
            print('\n\tSPAdes produced ' + str(count) + tr_str)
            pe_fastq_files[pe]['spades_assembly'] = assmbl_path
        else:
            print('\n\tSPAdes produced no transcripts.')
        print()


def makeblastdb_assemblies(assemblies, dir_prj_blast_assmbl, makeblastdb):  # noqa

    if len(assemblies) > 0:
        print('Building BLAST databases for assemblies:\n')
    for a in assemblies:
        assmbl_name = a['name']

        assmbl_blast_db_dir = opj(dir_prj_blast_assmbl, assmbl_name)
        assmbl_blast_db_file = opj(assmbl_blast_db_dir, assmbl_name)

        a['blast_db_path'] = assmbl_blast_db_file

        if ope(assmbl_blast_db_dir):
            print('\tBLAST database for ' + assmbl_name + ' already exists.')
        else:
            print('\t' + assmbl_name)
            make_dir(assmbl_blast_db_dir)
            make_blast_db(exec_file=makeblastdb,
                          in_file=a['path'],
                          out_file=assmbl_blast_db_file,
                          title=assmbl_name)

    if len(assemblies) > 0:
        print()


def run_tblastn_on_assemblies(assemblies, aa_queries_file, tblastn,
                              dir_prj_assmbl_blast_results, blast_2_evalue,
                              blast_2_max_target_seqs, blast_2_culling_limit,
                              blast_2_qcov_hsp_perc, threads, genetic_code):  # noqa
    if len(assemblies) > 0:
        print('Running BLAST on assemblies:\n')

    for a in assemblies:

        assmbl_name = a['name']
        assmbl_blast_db_path = a['blast_db_path']

        __ = opj(dir_prj_assmbl_blast_results, assmbl_name + '.tsv')

        if ope(__):
            print('\tBLAST results for assembly ' + assmbl_name +
                  ' already exist.')
        else:

            print('\tRunning tblastn on: ' + assmbl_name)

            run_blast(exec_file=tblastn,
                      task='tblastn',
                      threads=threads,
                      db_path=assmbl_blast_db_path,
                      queries_file=aa_queries_file,
                      out_file=__,
                      evalue=blast_2_evalue,
                      qcov_hsp_perc=blast_2_qcov_hsp_perc,
                      culling_limit=blast_2_culling_limit,
                      max_target_seqs=blast_2_max_target_seqs,
                      db_genetic_code=genetic_code,
                      out_cols=BLST_RES_COLS_2)

        a['blast_hits_aa'] = parse_blast_results_file(__, BLST_RES_COLS_2)

    if len(assemblies) > 0:
        print()


def find_orfs_translate(assemblies, dir_prj_transcripts, gc_tt, seqtk,
                        dir_temp, prepend_assmbl, min_target_orf_len,
                        max_target_orf_len, allow_non_aug, allow_no_strt_cod,
                        allow_no_stop_cod):  # noqa

    if len(assemblies) > 0:
        print('Analyzing BLAST hits for assemblies:\n')

    for a in assemblies:

        assmbl_name = a['name']

        parsed_hits = a['blast_hits_aa']
        a_path = a['path']

        transcripts_nt_fasta_file = opj(
            dir_prj_transcripts, assmbl_name + '_transcripts_nt.fasta')

        transcripts_nt_orf_fasta_file = opj(
            dir_prj_transcripts, assmbl_name + '_transcripts_nt_orf.fasta')

        transcripts_aa_orf_fasta_file = opj(
            dir_prj_transcripts, assmbl_name + '_transcripts_aa_orf.fasta')

        transcripts_nt = []
        transcripts_nt_orf = []
        transcripts_aa_orf = []

        a['annotations'] = {}

        collated = collate_blast_results(parsed_hits)

        ######################################################################
        # Use seqtk to sample the assembly FASTA file for sequences with
        # BLAST hits. This increases the speed substantially when the assembly
        # file is large.
        temp_a_file = opj(dir_temp, 'temp.fasta')
        temp_s_file = opj(dir_temp, 'temp.txt')
        sseqids_subsample = []
        for hit in collated:
            target_name = hit['sseqid']
            sseqids_subsample.append(target_name)
        sseqids_subsample_text = '\n'.join(sseqids_subsample)
        with open(temp_s_file, 'w') as f:
            f.write(sseqids_subsample_text)
        seqtk_extract_reads(seqtk,
                            in_file=a_path,
                            out_file=temp_a_file,
                            ids_file=temp_s_file)

        with open(temp_a_file, 'r') as f:
            __ = f.read()

        if __.strip() == '':
            continue

        print('\t' + assmbl_name)
        print('\t' + '-' * 80)

        __ = trim_desc_to_first_space_in_fasta_text(__)

        parsed_fasta = parse_fasta_text(__)
        ######################################################################

        all_kakapo_results = {}
        json_dump_file_path = opj(dir_prj_transcripts, assmbl_name +
                                  '_ann_kakapo.json')

        for hit in collated:

            target_name = hit['sseqid']
            target_seq = parsed_fasta[target_name]

            # Prepend assembly name to the sequence name:
            if prepend_assmbl is True:
                target_name = assmbl_name + '__' + target_name

            hit_start = hit['start']
            hit_end = hit['end']
            hit_frame = hit['frame']

            if allow_non_aug is True:
                start_codons = gc_tt['start_codons']
            else:
                start_codons = ['ATG']

            stop_codons = gc_tt['stop_codons']

            orf = find_orf_for_blast_hit(
                seq=target_seq,
                frame=hit_frame,
                hit_start=hit_start,
                hit_end=hit_end,
                hit_length_adjust=0.75,
                stop_codons=stop_codons,
                start_codons=start_codons,
                include_terminal_codon=True)

            if hit_frame > 0:
                ann_hit_b = hit_start
                ann_hit_e = hit_end
            else:
                target_seq = reverse_complement(target_seq)
                ann_hit_b = len(target_seq) - hit_start
                ann_hit_e = len(target_seq) - hit_end
                target_name = target_name + '__revcomp'

            a['annotations'][target_name] = {}

            if orf is not None:

                if hit_frame > 0:
                    orf_seq = target_seq[orf[0]:orf[1]]
                    ann_orf_b = orf[0]
                    ann_orf_e = orf[1]
                else:
                    ann_orf_b = len(target_seq) - orf[1]
                    ann_orf_e = len(target_seq) - orf[0]
                    orf_seq = target_seq[ann_orf_b:ann_orf_e]

                ##############################################################
                good_orf = True

                if allow_non_aug is False and \
                        orf_seq[0:3] != 'ATG':

                    good_orf = False
                    # print('allow_non_aug', allow_non_aug, orf_seq[0:3])

                elif allow_no_strt_cod is False and \
                        orf_seq[0:3] not in start_codons:

                    good_orf = False
                    # print('allow_no_strt_cod', allow_no_strt_cod,
                    #       orf_seq[0:3])

                elif allow_no_stop_cod is False and \
                        orf_seq[-3:] not in stop_codons:

                    good_orf = False
                    # print('allow_no_stop_cod', allow_no_stop_cod,
                    #       orf_seq[-3:0])

                elif len(orf_seq) < min_target_orf_len:
                    good_orf = False
                    # print('min_target_orf_len', min_target_orf_len,
                    #       len(orf_seq))

                elif len(orf_seq) > max_target_orf_len:
                    good_orf = False
                    # print('max_target_orf_len', max_target_orf_len,
                    #       len(orf_seq))
                ##############################################################

                if good_orf is True:

                    a['annotations'][target_name]['orf_begin'] = ann_orf_b
                    a['annotations'][target_name]['orf_end'] = ann_orf_e

                    transcripts_nt_orf.append({'name': target_name,
                                               'seq': orf_seq})

                    transl_seq = translate(orf_seq, gc_tt)
                    transcripts_aa_orf.append({'name': target_name,
                                               'seq': transl_seq})

            transcripts_nt.append({'name': target_name, 'seq': target_seq})

            a['annotations'][target_name]['frame'] = abs(hit_frame)
            a['annotations'][target_name]['blast_hit_begin'] = ann_hit_b
            a['annotations'][target_name]['blast_hit_end'] = ann_hit_e

            ##################################################################
            # Collect ORF and BLAST hit annotations for downstream use.
            kakapo_json = [{}]
            kakapo_json[0]['kakapo_annotations'] = (
                a['annotations'][target_name])
            all_kakapo_results[target_name] = kakapo_json
            ##################################################################

        # --------------------------------------------------------------------

        print('\t                     Transcripts: ' +
              str(len(transcripts_nt)))

        if len(transcripts_nt) > 0:
            write_fasta_file(transcripts_nt, transcripts_nt_fasta_file)
            a['transcripts_nt_fasta_file'] = transcripts_nt_fasta_file
        else:
            a['transcripts_nt_fasta_file'] = None

        print('\tTranscripts with acceptable ORFs: ' +
              str(len(transcripts_nt_orf)))

        if len(transcripts_nt_orf) > 0:
            write_fasta_file(transcripts_nt_orf, transcripts_nt_orf_fasta_file)
            a['transcripts_nt_orf_fasta_file'] = transcripts_nt_orf_fasta_file
        else:
            a['transcripts_nt_orf_fasta_file'] = None

        if len(transcripts_aa_orf) > 0:
            write_fasta_file(transcripts_aa_orf, transcripts_aa_orf_fasta_file)
            a['transcripts_aa_orf_fasta_file'] = transcripts_aa_orf_fasta_file
        else:
            a['transcripts_aa_orf_fasta_file'] = None

        # --------------------------------------------------------------------
        # Save ORF and BLAST hit annotations for downstream use.
        with open(json_dump_file_path, 'w') as f:
            json.dump(all_kakapo_results, f, sort_keys=True, indent=4)
        # --------------------------------------------------------------------

        print('\t' + '-' * 80 + '\n\n')

        # --------------------------------------------------------------------


def run_inter_pro_scan(assemblies, email, dir_prj_ips, dir_cache_prj):  # noqa

    delay = 1

    if len(assemblies) > 0:
        print('\nRunning InterProScan on translated transcripts:\n')

    for a in assemblies:

        if 'transcripts_aa_orf_fasta_file' not in a:
            continue

        aa_file = a['transcripts_aa_orf_fasta_file']

        if aa_file is None:
            continue

        assmbl_name = a['name']

        json_dump_file_path = opj(dir_prj_ips, assmbl_name + '_ann_ips.json')

        if ope(json_dump_file_path):
            print('\tInterProScan results for assembly ' + assmbl_name +
                  ' have already been downloaded.')
            continue

        seqs = read_fasta_file_dict(aa_file)

        __ = opj(dir_cache_prj, assmbl_name + '_ips_jobs')

        if ope(__):
            with open(__, 'rb') as f:
                jobs = pickle.load(f)

        else:
            print()
            jobs = job_runner(email=email, dir_cache=dir_cache_prj, seqs=seqs)

            with open(__, 'wb') as f:
                pickle.dump(jobs, f, protocol=PICKLE_PROTOCOL)

        print('\n\tDownloading InterProScan results for transcripts in ' +
              assmbl_name + '\n')

        all_ips_results = {}

        for i, job in enumerate(jobs['finished']):
            sleep(delay)
            job_id = jobs['finished'][job]
            ips_json = result_json(job_id)
            # ips_version = ips_json['interproscan-version']
            ips_json = ips_json['results']

            # These fields are set to 'EMBOSS_001' by default
            # Delete them
            del ips_json[0]['xref']

            # kakapo annotations
            ips_json[0]['kakapo_annotations'] = a['annotations'][job]

            all_ips_results[job] = ips_json

        with open(json_dump_file_path, 'w') as f:
            json.dump(all_ips_results, f, sort_keys=True, indent=4)

        # Removes cached jobs file.
        # ToDo: Check if there are no failed jobs, before deleting
        osremove(__)


def gff_from_json(assemblies, dir_prj_ips, dir_prj_transcripts_combined,
                  prj_name):  # noqa

    if len(assemblies) > 0:
        print('\nProducing GFF3 files:\n')

    all_fas_paths = []
    all_gff_paths = []

    combined_fas_path = opj(dir_prj_transcripts_combined, prj_name + '.fasta')
    combined_gff_path = opj(dir_prj_transcripts_combined, prj_name + '.gff')

    for a in assemblies:

        if 'transcripts_nt_fasta_file' not in a:
            continue

        assmbl_name = a['name']
        transcripts_nt_path = a['transcripts_nt_fasta_file']

        kakapo_json_path = opj(dir_prj_ips, assmbl_name + '_ann_kakapo.json')
        ips_json_path = opj(dir_prj_ips, assmbl_name + '_ann_ips.json')
        json_path = opj(dir_prj_ips, assmbl_name + '_ann.json')

        gff_path = transcripts_nt_path.replace('.fasta', '.gff')

        ips_json_dict = {}
        kakapo_json_dict = {}

        if ope(ips_json_path):
            with open(ips_json_path, 'r') as f:
                ips_json_dict = json.load(f)

        if ope(kakapo_json_path):
            with open(kakapo_json_path, 'r') as f:
                kakapo_json_dict = json.load(f)

        json_dict = kakapo_json_dict.copy()
        json_dict.update(ips_json_dict)

        with open(json_path, 'w') as f:
            json.dump(json_dict, f, sort_keys=True, indent=4)

        if ope(json_path):
            print('\t' + assmbl_name)
            gff_from_kakapo_ips5_json_file(json_path, gff_path)

            osremove(json_path)

            all_gff_paths.append(gff_path)
            all_fas_paths.append(transcripts_nt_path)

    combine_text_files(all_fas_paths, combined_fas_path)
    combine_text_files(all_gff_paths, combined_gff_path)
