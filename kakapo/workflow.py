#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""kakapo workflow"""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import with_statement

import fileinput
import json
import pickle
import re

from functools import partial
from os import remove as osremove
from os import stat as osstat
from os.path import basename
from os.path import commonprefix
from os.path import exists as ope
from os.path import join as opj
from os.path import splitext
from shutil import copyfile
from shutil import rmtree
from sys import exit
from time import sleep

from kakapo.bioio import filter_fasta_text_by_length
from kakapo.bioio import read_fasta
from kakapo.bioio import standardize_fasta_text
from kakapo.bioio import trim_desc_to_first_space_in_fasta_text
from kakapo.bioio import write_fasta
from kakapo.blast import BLST_RES_COLS_1, BLST_RES_COLS_2
from kakapo.blast import collate_blast_results
from kakapo.blast import make_blast_db, run_blast
from kakapo.blast import parse_blast_results_file
from kakapo.config import PICKLE_PROTOCOL
from kakapo.data.start_codon_context import contexts as atg_contexts
from kakapo.ebi_domain_search import pfam_entry
from kakapo.ebi_domain_search import pfam_seqs
from kakapo.ebi_domain_search import prot_ids_for_tax_ids
from kakapo.ebi_iprscan5 import job_runner
from kakapo.ebi_iprscan5 import result_json
from kakapo.ebi_proteins import fasta_by_accession_list
from kakapo.entrez import cds_acc_for_prot_acc
from kakapo.entrez import dnld_cds_nt_fasta as dnld_ncbi_cds_nt_fasta
from kakapo.entrez import dnld_seqs as dnld_ncbi_seqs
from kakapo.entrez import sra_run_info
from kakapo.entrez import summary as entrez_summary
from kakapo.entrez import taxids_for_acc
from kakapo.gff3 import gff_from_kakapo_ips5_json_file
from kakapo.helpers import combine_text_files
from kakapo.helpers import keep_unique_lines_in_file
from kakapo.helpers import make_dir
from kakapo.orf import find_orf_for_blast_hit
from kakapo.py_v_diffs import StringIO
from kakapo.rcorrector import filter_unc_se, filter_unc_pe
from kakapo.rcorrector import run_rcorrector_se, run_rcorrector_pe
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

    dir_prj = opj(dir_out, '02-project-specific', prj_name)
    make_dir(dir_prj)

    dir_prj_logs = opj(dir_prj, '00-logs')
    make_dir(dir_prj_logs)

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

    dir_global = opj(dir_out, '01-global')
    make_dir(dir_global)

    dir_fq_data = opj(dir_global, '01-sra-fq-data')
    make_dir(dir_fq_data)

    dir_fq_cor_data = opj(dir_global, '02-corrected-fq-data')
    make_dir(dir_fq_cor_data)

    dir_fq_trim_data = opj(dir_global, '03-trimmed-fq-data')
    make_dir(dir_fq_trim_data)

    dir_fq_filter_data = opj(dir_global, '04-filtered-fq-data')
    make_dir(dir_fq_filter_data)

    dir_fa_trim_data = opj(dir_global, '05-fa-data')
    make_dir(dir_fa_trim_data)

    dir_blast_fa_trim = opj(dir_global, '06-fa-blast-db-data')
    make_dir(dir_blast_fa_trim)

    ret_dict = {'dir_blast_fa_trim': dir_blast_fa_trim,
                'dir_cache': dir_cache,
                'dir_cache_fq_minlen': dir_cache_fq_minlen,
                'dir_cache_pfam_acc': dir_cache_pfam_acc,
                'dir_cache_prj': dir_cache_prj,
                'dir_fa_trim_data': dir_fa_trim_data,
                'dir_fq_cor_data': dir_fq_cor_data,
                'dir_fq_data': dir_fq_data,
                'dir_fq_trim_data': dir_fq_trim_data,
                'dir_fq_filter_data': dir_fq_filter_data,
                'dir_prj': dir_prj,
                'dir_prj_logs': dir_prj_logs,
                'dir_prj_assmbl_blast_results': dir_prj_assmbl_blast_results,
                'dir_prj_blast_assmbl': dir_prj_blast_assmbl,
                'dir_prj_blast_results_fa_trim': dir_prj_blast_results_fa_trim,
                'dir_prj_ips': dir_prj_ips,
                'dir_prj_queries': dir_prj_queries,
                'dir_prj_spades_assemblies': dir_prj_spades_assemblies,
                'dir_prj_transcripts': dir_prj_transcripts,
                'dir_prj_transcripts_combined': dir_prj_transcripts_combined,
                'dir_prj_vsearch_results_fa_trim':
                    dir_prj_vsearch_results_fa_trim,
                'dir_temp': dir_temp}

    return ret_dict

def descending_tax_ids(tax_ids_user, taxonomy, linfo=print):  # noqa
    # if len(tax_ids_user) > 0:
    #     linfo('Resolving descending nodes for TaxIds')
    shared = taxonomy.shared_taxid_for_taxids(tax_ids_user)
    if shared is None:
        return None
    tax_ids = taxonomy.all_descending_taxids(taxid=shared)
    if tax_ids is None:
        tax_ids = [shared]
    tax_ids = [int(x) for x in tax_ids]
    return tax_ids


def pfam_uniprot_accessions(pfam_acc, tax_ids, dir_cache_pfam_acc,
                            linfo=print):  # noqa
    if len(pfam_acc) > 0:
        linfo('Downloading UniProt accessions for Pfam families')
    pfam_seqs_list = []
    for pa in pfam_acc:
        pfam_id = pfam_entry(pa)[0]['id']
        linfo(pa + ': ' + pfam_id)
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
    return pfam_uniprot_acc


def dnld_pfam_uniprot_seqs(uniprot_acc, aa_uniprot_file, dir_cache_prj,
                           linfo=print):  # noqa
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

            linfo('Downloading Pfam protein sequences from UniProt')
            __ = fasta_by_accession_list(uniprot_acc)
            __ = standardize_fasta_text(__)

            with open(aa_uniprot_file, 'w') as f:
                f.write(__)
    else:
        if ope(aa_uniprot_file):
            osremove(aa_uniprot_file)


def user_protein_accessions(prot_acc_user, linfo=print):  # noqa
    if len(prot_acc_user) > 0:
        linfo('Reading user provided protein accessions')
        pa_info = entrez_summary(prot_acc_user, 'protein')
        prot_acc = []
        for pa in pa_info:
            title = pa['title']
            title = title[0].upper() + title[1:]
            acc = pa['accessionversion']
            prot_acc.append(acc)

            if len(title) > 60:
                title = title[0:57] + '...'

            linfo(acc + ': ' + title)

        return prot_acc

    else:
        return prot_acc_user


def dnld_prot_seqs(prot_acc_user, aa_prot_ncbi_file, dir_cache_prj,
                   linfo=print):  # noqa
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

            linfo('Downloading protein sequences from NCBI')
            __ = dnld_ncbi_seqs(prot_acc_user, 'protein')

            write_fasta(__, aa_prot_ncbi_file)

    else:
        if ope(aa_prot_ncbi_file):
            osremove(aa_prot_ncbi_file)


def user_aa_fasta(user_queries, aa_prot_user_file, linfo=print):  # noqa
    __ = ''
    if len(user_queries) > 0:
        linfo('Reading user provided AA sequences')
        for ap in user_queries:
            linfo(ap)
            with open(ap, 'r') as f:
                __ = __ + f.read()
    if __ != '':
        with open(aa_prot_user_file, 'w') as f:
            f.write(standardize_fasta_text(__))


def combine_aa_fasta(fasta_files, aa_queries_file, linfo=print):  # noqa
    linfo('Combining all AA query sequences')
    __ = ''
    for fasta_file in fasta_files:
        if ope(fasta_file):
            with open(fasta_file, 'r') as f:
                __ = __ + f.read()

    if __ != '':
        with open(aa_queries_file, 'w') as f:
            f.write(__)
    else:
        linfo('No queries were provided. Exiting.')
        exit(0)


def filter_queries(aa_queries_file, min_query_length, max_query_length,
                   linfo=print): # noqa
    __ = ''
    with open(aa_queries_file, 'r') as f:
        __ = f.read()

    linfo('Filtering AA query sequences')
    linfo('min_query_length: ' + str(min_query_length))
    linfo('max_query_length: ' + str(max_query_length))

    __ = filter_fasta_text_by_length(__, min_query_length, max_query_length)

    with open(aa_queries_file, 'w') as f:
        f.write(__)


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

    for se in fq_se:
        tax_id = se[0]
        path = se[1]
        sample_base_name = splitext(basename(path))[0].split('_R')[0]
        se_fastq_files[sample_base_name] = {'path': path}
        se_fastq_files[sample_base_name]['src'] = 'usr'
        se_fastq_files[sample_base_name]['avg_len'] = None
        se_fastq_files[sample_base_name]['tax_id'] = tax_id
        linfo(sample_base_name + ': ' + path)

    for pe in fq_pe:
        tax_id = pe[0]
        path = pe[1]
        sample_base_name = basename(commonprefix(path)).rstrip('_- R')
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
    for se in se_fastq_files:
        dir_fq_cor_data_sample = opj(dir_fq_cor_data, se)
        fq_path = se_fastq_files[se]['path']
        log_f = opj(dir_fq_cor_data_sample, se + '.txt')
        out_f = opj(dir_fq_cor_data_sample, se + '.fastq')
        se_fastq_files[se]['cor_path_fq'] = out_f

        if ope(dir_fq_cor_data_sample):
            linfo('Corrected FASTQ file for sample ' + se + ' already exist')
        else:
            make_dir(dir_fq_cor_data_sample)
            linfo('Running Rcorrector in SE mode: ' + se)
            run_rcorrector_se(rcorrector=rcorrector,
                              in_file=fq_path,
                              out_dir=dir_fq_cor_data_sample,
                              threads=threads,
                              dir_temp=dir_temp)

            fq_base_path = opj(dir_fq_cor_data_sample, basename(fq_path))
            fq_cor_path = splitext(fq_base_path)[0] + '.cor.fq'

            filter_unc_se(in_file=fq_cor_path, out_file=out_f, log_file=log_f)

            osremove(fq_cor_path)

    for pe in pe_fastq_files:
        dir_fq_cor_data_sample = opj(dir_fq_cor_data, pe)
        log_f = opj(dir_fq_cor_data_sample, pe + '.txt')
        out_f_1 = opj(dir_fq_cor_data_sample, pe + '_R1.fastq')
        out_f_2 = opj(dir_fq_cor_data_sample, pe + '_R2.fastq')
        pe_fastq_files[pe]['cor_path_fq'] = [out_f_1, out_f_2]

        fq_path_1 = pe_fastq_files[pe]['path'][0]
        fq_path_2 = pe_fastq_files[pe]['path'][1]
        fq_path_3 = None
        out_f_3 = None
        if len(pe_fastq_files[pe]['path']) == 3:
            fq_path_3 = pe_fastq_files[pe]['path'][2]
            out_f_3 = opj(dir_fq_cor_data_sample, pe + '_R3.fastq')
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
            fq_cor_path_1 = splitext(fq_base_path_1)[0] + '.cor.fq'
            fq_base_path_2 = opj(dir_fq_cor_data_sample, basename(fq_path_2))
            fq_cor_path_2 = splitext(fq_base_path_2)[0] + '.cor.fq'

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
                fq_cor_path_3 = splitext(fq_base_path_3)[0] + '.cor.fq'
                log_f_3 = opj(dir_fq_cor_data_sample, pe + '_unpaired.txt')

                filter_unc_se(in_file=fq_cor_path_3, out_file=out_f_3,
                              log_file=log_f_3)

                osremove(fq_cor_path_3)


def run_trimmomatic(se_fastq_files, pe_fastq_files, dir_fq_trim_data,
                    trimmomatic, adapters, fpatt, threads, linfo=print):  # noqa
    for se in se_fastq_files:
        dir_fq_trim_data_sample = opj(dir_fq_trim_data, se)
        fq_path = se_fastq_files[se]['cor_path_fq']
        min_acc_len = se_fastq_files[se]['min_acc_len']
        stats_f = opj(dir_fq_trim_data_sample, se + '.txt')
        out_f = opj(dir_fq_trim_data_sample, se + '.fastq')
        se_fastq_files[se]['trim_path_fq'] = out_f

        if ope(dir_fq_trim_data_sample):
            linfo('Trimmed FASTQ files for sample ' + se + ' already exist')
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
        if len(pe_fastq_files[pe]['cor_path_fq']) == 3:
            fq_path_3 = pe_fastq_files[pe]['cor_path_fq'][2]
        min_acc_len = pe_fastq_files[pe]['min_acc_len']
        stats_f = opj(dir_fq_trim_data_sample, pe + '.txt')
        out_fs = [x.replace('@D@', dir_fq_trim_data_sample) for x in fpatt]
        out_fs = [x.replace('@N@', pe) for x in out_fs]
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

                out_f = opj(dir_fq_trim_data_sample, 'unpaired.fastq')
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

                _ = opj(dir_fq_trim_data_sample, 'temp.fastq')
                f_temp = open(_, 'a')
                with fileinput.input(files=[out_fs[2], out_f]) as f:
                    for line in f:
                        f_temp.write(line)
                f_temp.close()

                osremove(out_fs[2])
                osremove(out_f)
                copyfile(_, out_fs[2])
                osremove(_)


def trimmed_fq_to_fa(se_fastq_files, pe_fastq_files, dir_fa_trim_data, seqtk,
                     fpatt, linfo=print): # noqa
    for se in se_fastq_files:
        dir_fa_trim_data_sample = opj(dir_fa_trim_data, se)
        fq_path = se_fastq_files[se]['trim_path_fq']
        out_f = opj(dir_fa_trim_data_sample, se + '.fasta')
        se_fastq_files[se]['trim_path_fa'] = out_f

        if ope(dir_fa_trim_data_sample):
            linfo('Trimmed FASTA files for sample ' + se + ' already exist')
        else:
            make_dir(dir_fa_trim_data_sample)
            linfo('Converting FASTQ to FASTA using Seqtk: ' + fq_path)
            seqtk_fq_to_fa(seqtk, fq_path, out_f)

    for pe in pe_fastq_files:
        dir_fa_trim_data_sample = opj(dir_fa_trim_data, pe)
        fq_paths = pe_fastq_files[pe]['trim_path_fq']
        out_fs = [x.replace('@D@', dir_fa_trim_data_sample) for x in fpatt]
        out_fs = [x.replace('@N@', pe) for x in out_fs]
        pe_fastq_files[pe]['trim_path_fa'] = out_fs

        if ope(dir_fa_trim_data_sample):
            linfo('Trimmed FASTA files for sample ' + pe + ' already exist')
        else:
            make_dir(dir_fa_trim_data_sample)
            pe_trim_files = zip(fq_paths, out_fs)
            for x in pe_trim_files:
                linfo('Converting FASTQ to FASTA using Seqtk: ' + x[0])
                seqtk_fq_to_fa(seqtk, x[0], x[1])


def makeblastdb_fq(se_fastq_files, pe_fastq_files, dir_blast_fa_trim,
                   makeblastdb, fpatt, linfo=print): # noqa
    for se in se_fastq_files:
        dir_blast_fa_trim_sample = opj(dir_blast_fa_trim, se)
        fa_path = se_fastq_files[se]['trim_path_fa']
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
        fa_paths = pe_fastq_files[pe]['trim_path_fa']
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
                         tblastn, blast_1_evalue, blast_1_max_target_seqs,
                         blast_1_culling_limit, blast_1_qcov_hsp_perc,
                         dir_blast_results_fa_trim, fpatt, threads,
                         genetic_code, seqtk, vsearch, linfo=print): # noqa
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
            linfo('BLAST results for sample ' + se + ' already exists')
        else:
            make_dir(dir_results)
            linfo('Running tblastn on: ' + blast_db_path)
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

            linfo('Extracting unique BLAST hits using Seqtk')

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
            linfo('BLAST results for sample ' + pe + ' already exist')
        else:
            make_dir(dir_results)
            pe_trim_files = zip(blast_db_paths, out_fs, fq_paths, out_fs_fastq,
                                out_fs_fasta)
            for x in pe_trim_files:
                linfo('Running tblastn on: ' + x[0])
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

                linfo('\tExtracting unique BLAST hits using Seqtk')

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
                         dir_vsearch_results_fa_trim, fpatt, seqtk,
                         linfo=print): # noqa
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
            linfo('Vsearch results for sample ' + se + ' already exists')
        else:
            make_dir(dir_results)
            linfo('Running vsearch on: ' + fq_path)
            run_vsearch(vsearch,
                        ident=ident,
                        q_file=blast_results_fa_path,
                        db_file=fq_path,
                        out_file=out_f,
                        minlen=min_acc_len)

            linfo('Extracting unique vsearch hits using Seqtk')
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
            linfo('Vsearch results for sample ' + pe + ' already exist')
        else:
            make_dir(dir_results)
            pe_trim_files = zip(fq_paths, out_fs, out_fs_fastq)
            for x in pe_trim_files:
                linfo('Running vsearch on: ' + x[0])
                run_vsearch(vsearch,
                            ident=ident,
                            q_file=blast_results_fa_path,
                            db_file=x[0],
                            out_file=x[1],
                            minlen=min_acc_len)

            linfo('Extracting unique vsearch hits from paired files '
                  'using Seqtk')

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

            linfo('Extracting unique vsearch hits from unpaired files '
                  'using Seqtk')

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
               spades, dir_temp, threads, ram, linfo=print):  # noqa
    for se in se_fastq_files:
        dir_results = opj(dir_spades_assemblies, se)
        fq_path = se_fastq_files[se]['vsearch_results_path']
        se_fastq_files[se]['spades_assembly'] = None

        if ope(dir_results):
            linfo('SPAdes results for sample ' + se + ' already exist')
        else:
            make_dir(dir_results)
            linfo('Running SPAdes on: ' + se)
            run_spades_se(spades,
                          out_dir=dir_results,
                          input_file=fq_path,
                          threads=threads,
                          memory=ram,
                          rna=True)

        assmbl_path = opj(dir_results, 'transcripts.fasta')
        if ope(assmbl_path):
            count = len(read_fasta(assmbl_path))
            tr_str = ' transcripts'
            if count == 1:
                tr_str = ' transcript'
            linfo('SPAdes produced ' + str(count) + tr_str)
            se_fastq_files[se]['spades_assembly'] = assmbl_path
        else:
            linfo('SPAdes produced no transcripts')

    for pe in pe_fastq_files:
        dir_results = opj(dir_spades_assemblies, pe)
        fq_paths = pe_fastq_files[pe]['vsearch_results_path']
        pe_fastq_files[pe]['spades_assembly'] = None

        if ope(dir_results):
            linfo('SPAdes results for sample ' + pe + ' already exist')
        else:
            make_dir(dir_results)
            linfo('Running SPAdes on: ' + pe)

            if osstat(fq_paths[0]).st_size > 0 and \
               osstat(fq_paths[1]).st_size > 0:

                run_spades_pe(spades,
                              out_dir=dir_results,
                              input_files=fq_paths,
                              threads=threads,
                              memory=ram,
                              rna=True)

            else:
                _ = opj(dir_temp, 'temp.fasta')
                combine_text_files(fq_paths, _)
                run_spades_se(spades,
                              out_dir=dir_results,
                              input_file=_,
                              threads=threads,
                              memory=ram,
                              rna=True)
                osremove(_)

        assmbl_path = opj(dir_results, 'transcripts.fasta')
        if ope(assmbl_path):
            count = len(read_fasta(assmbl_path))
            tr_str = ' transcripts'
            if count == 1:
                tr_str = ' transcript'
            linfo('SPAdes produced ' + str(count) + tr_str)
            pe_fastq_files[pe]['spades_assembly'] = assmbl_path
        else:
            linfo('SPAdes produced no transcripts')


def makeblastdb_assemblies(assemblies, dir_prj_blast_assmbl, makeblastdb,
                           linfo=print):  # noqa
    if len(assemblies) > 0:
        linfo('Building BLAST databases for assemblies')
    for a in assemblies:
        assmbl_name = a['name']

        assmbl_blast_db_dir = opj(dir_prj_blast_assmbl, assmbl_name)
        assmbl_blast_db_file = opj(assmbl_blast_db_dir, assmbl_name)

        a['blast_db_path'] = assmbl_blast_db_file

        if ope(assmbl_blast_db_dir):
            linfo('BLAST database for ' + assmbl_name + ' already exists')
        else:
            linfo(assmbl_name)
            make_dir(assmbl_blast_db_dir)
            make_blast_db(exec_file=makeblastdb,
                          in_file=a['path'],
                          out_file=assmbl_blast_db_file,
                          title=assmbl_name)


def run_tblastn_on_assemblies(assemblies, aa_queries_file, tblastn,
                              dir_prj_assmbl_blast_results, blast_2_evalue,
                              blast_2_max_target_seqs, blast_2_culling_limit,
                              blast_2_qcov_hsp_perc, threads, genetic_code,
                              linfo=print):  # noqa
    if len(assemblies) > 0:
        linfo('Running BLAST on assemblies')

    for a in assemblies:

        assmbl_name = a['name']
        assmbl_blast_db_path = a['blast_db_path']

        __ = opj(dir_prj_assmbl_blast_results, assmbl_name + '.tsv')

        if ope(__):
            linfo('BLAST results for assembly ' + assmbl_name +
                  ' already exist')
        else:
            linfo('Running tblastn on: ' + assmbl_name)

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


def find_orfs_translate(assemblies, dir_prj_transcripts, gc_tt, seqtk,
                        dir_temp, prepend_assmbl, min_target_orf_len,
                        max_target_orf_len, allow_non_aug, allow_no_strt_cod,
                        allow_no_stop_cod, tax, tax_group, tax_ids_user,
                        linfo=print):  # noqa
    if len(assemblies) > 0:
        linfo('Analyzing BLAST hits for assemblies')

    for a in assemblies:

        assmbl_name = a['name']
        tax_id = a['tax_id']

        parsed_hits = a['blast_hits_aa']
        a_path = a['path']

        transcripts_nt_fasta_file = opj(
            dir_prj_transcripts, assmbl_name + '_transcripts_nt.fasta')

        transcripts_nt_orf_fasta_file = opj(
            dir_prj_transcripts, assmbl_name + '_transcripts_nt_orf.fasta')

        transcripts_aa_orf_fasta_file = opj(
            dir_prj_transcripts, assmbl_name + '_transcripts_aa_orf.fasta')

        transcripts_nt = {}
        transcripts_nt_orf = {}
        transcripts_aa_orf = {}

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

        linfo(assmbl_name)

        __ = trim_desc_to_first_space_in_fasta_text(__)

        parsed_fasta = read_fasta(StringIO(__))
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
                # Also prepend taxonomic info to the sequence name:
                if tax_id is not None:
                    fm = tax.higher_rank_for_taxid(tax_id, rank='family')
                    cn = tax.genbank_common_name_for_taxid(tax_id)

                    if fm is not None:
                        target_name = fm + '__' + target_name

                    if cn is not None:
                        cn = cn.lower()
                        cn = cn.replace(' ', '_')
                        target_name = cn + '__' + target_name

            hit_start = hit['start']
            hit_end = hit['end']
            hit_frame = hit['frame']

            if allow_non_aug is True:
                start_codons = gc_tt.start_codons_ambiguous
            else:
                start_codons = ['ATG']

            stop_codons = gc_tt.stop_codons_ambiguous

            ##################################################################
            if tax_id is not None:
                tax_ids_for_orf = (tax_id, )
            else:
                tax_ids_for_orf = tax_ids_user

            cntx_txids_avail = tuple(
                sorted(set(map(lambda x: int(x.split('_')[0]),
                               atg_contexts.keys()))))

            cntx_taxid = set()
            for txid in tax_ids_for_orf:
                tax_path = partial(tax.path_between_taxids, txid)
                path_len = tuple(map(len,
                                     tuple(map(tax_path, cntx_txids_avail))))
                cntx_taxid.add(cntx_txids_avail[path_len.index(min(path_len))])
            cntx_taxid = tuple(cntx_taxid)[0]

            cntx_l_key = str(cntx_taxid) + '_L'
            cntx_r_key = str(cntx_taxid) + '_R'

            cntx_l = atg_contexts[cntx_l_key]
            cntx_r = atg_contexts[cntx_r_key]
            ##################################################################

            orf = find_orf_for_blast_hit(
                seq=target_seq,
                frame=hit_frame,
                hit_start=hit_start,
                hit_end=hit_end,
                stop_codons=stop_codons,
                start_codons=start_codons,
                context_l=cntx_l,
                context_r=cntx_r)

            if hit_frame > 0:
                ann_hit_b = hit_start
                ann_hit_e = hit_end
            else:
                target_seq = reverse_complement(target_seq)
                ann_hit_b = len(target_seq) - hit_start
                ann_hit_e = len(target_seq) - hit_end
                target_name = target_name + '__revcomp'

            a['annotations'][target_name] = {}

            good_orf = orf[0]
            if good_orf is not None:

                print(' :: ' + target_name)

                if hit_frame > 0:
                    ann_orf_b = good_orf[0]
                    ann_orf_e = good_orf[1] + 3
                    orf_seq = target_seq[ann_orf_b:ann_orf_e]
                else:
                    ann_orf_b = len(target_seq) - good_orf[1]
                    ann_orf_e = len(target_seq) - good_orf[0] + 3
                    orf_seq = target_seq[ann_orf_b:ann_orf_e]

                ##############################################################
                valid_orf = True

                if allow_non_aug is False and \
                        orf_seq[0:3] != 'ATG':

                    valid_orf = False

                elif allow_no_strt_cod is False and \
                        orf_seq[0:3] not in start_codons:

                    valid_orf = False

                elif allow_no_stop_cod is False and \
                        orf_seq[-3:] not in stop_codons:

                    valid_orf = False

                elif len(orf_seq) < min_target_orf_len:
                    valid_orf = False

                elif len(orf_seq) > max_target_orf_len:
                    valid_orf = False
                ##############################################################

                if valid_orf is True:

                    a['annotations'][target_name]['orf_begin'] = ann_orf_b
                    a['annotations'][target_name]['orf_end'] = ann_orf_e

                    transcripts_nt_orf[target_name] = orf_seq

                    transl_seq = translate(orf_seq,
                                           gc_tt.table_ambiguous,
                                           start_codons)

                    transcripts_aa_orf[target_name] = transl_seq

            bad_orfs = orf[1]
            if len(bad_orfs) > 0:
                a['annotations'][target_name]['orfs_bad'] = list()
                orfs_bad_list = a['annotations'][target_name]['orfs_bad']

            for bad_orf in bad_orfs:

                if hit_frame > 0:
                    ann_orf_b = bad_orf[0]
                    ann_orf_e = bad_orf[1] + 3
                    orf_seq = target_seq[ann_orf_b:ann_orf_e]
                else:
                    ann_orf_b = len(target_seq) - bad_orf[1]
                    ann_orf_e = len(target_seq) - bad_orf[0] + 3
                    orf_seq = target_seq[ann_orf_b:ann_orf_e]

                orf_bad_dict = dict()
                orf_bad_dict['orf_begin'] = ann_orf_b
                orf_bad_dict['orf_end'] = ann_orf_e
                orfs_bad_list.append(orf_bad_dict)

            transcripts_nt[target_name] = target_seq

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

        linfo('Transcripts: ' + str(len(transcripts_nt)))

        if len(transcripts_nt) > 0:
            write_fasta(transcripts_nt, transcripts_nt_fasta_file)
            a['transcripts_nt_fasta_file'] = transcripts_nt_fasta_file
        else:
            a['transcripts_nt_fasta_file'] = None

        linfo('Transcripts with acceptable ORFs: ' +
              str(len(transcripts_nt_orf)))

        if len(transcripts_nt_orf) > 0:
            write_fasta(transcripts_nt_orf, transcripts_nt_orf_fasta_file)
            a['transcripts_nt_orf_fasta_file'] = transcripts_nt_orf_fasta_file
        else:
            a['transcripts_nt_orf_fasta_file'] = None

        if len(transcripts_aa_orf) > 0:
            write_fasta(transcripts_aa_orf, transcripts_aa_orf_fasta_file)
            a['transcripts_aa_orf_fasta_file'] = transcripts_aa_orf_fasta_file
        else:
            a['transcripts_aa_orf_fasta_file'] = None

        # --------------------------------------------------------------------
        # Save ORF and BLAST hit annotations for downstream use.
        with open(json_dump_file_path, 'w') as f:
            json.dump(all_kakapo_results, f, sort_keys=True, indent=4)
        # --------------------------------------------------------------------


def run_inter_pro_scan(assemblies, email, dir_prj_ips, dir_cache_prj,
                       linfo=print):  # noqa
    delay = 1

    if len(assemblies) > 0:
        linfo('Running InterProScan on translated transcripts')

    for a in assemblies:

        if 'transcripts_aa_orf_fasta_file' not in a:
            continue

        aa_file = a['transcripts_aa_orf_fasta_file']

        if aa_file is None:
            continue

        assmbl_name = a['name']

        json_dump_file_path = opj(dir_prj_ips, assmbl_name + '_ann_ips.json')

        if ope(json_dump_file_path):
            linfo('InterProScan results for assembly ' + assmbl_name +
                  ' have already been downloaded')
            continue

        seqs = read_fasta(aa_file)

        __ = opj(dir_cache_prj, assmbl_name + '_ips_jobs')

        if ope(__):
            with open(__, 'rb') as f:
                jobs = pickle.load(f)

        else:
            jobs = job_runner(email=email, dir_cache=dir_cache_prj,
                              seqs=seqs, logger=linfo)

            with open(__, 'wb') as f:
                pickle.dump(jobs, f, protocol=PICKLE_PROTOCOL)

        linfo('Downloading InterProScan results for transcripts in ' +
              assmbl_name)

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
                  prj_name, linfo=print):  # noqa
    if len(assemblies) > 0:
        linfo('Producing GFF3 files')

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
            linfo(assmbl_name)
            gff_from_kakapo_ips5_json_file(json_path, gff_path)

            osremove(json_path)

            all_gff_paths.append(gff_path)
            all_fas_paths.append(transcripts_nt_path)

    combine_text_files(all_fas_paths, combined_fas_path)
    combine_text_files(all_gff_paths, combined_gff_path)


def dnld_cds_for_ncbi_prot_acc(prot_acc_user, nt_prot_ncbi_file, tax,
                               linfo=print):  # noqa
    linfo('Downloading CDS for user provided NCBI protein accessions')
    cds_acc_dict = cds_acc_for_prot_acc(prot_acc_user)

    cds_accessions = []
    for prot_acc in cds_acc_dict:
        cds_acc = cds_acc_dict[prot_acc]
        cds_accessions.append(cds_acc)

    temp = dnld_ncbi_cds_nt_fasta(cds_accessions)

    taxids = taxids_for_acc(prot_acc_user, 'protein')

    cds_seqs_fasta_list = []

    for x in temp:
        description = x.split('\n')[0]
        description = description.split('|')[1]
        prot_id = re.findall(r'\[protein_id=(.*?)\]', x)
        if len(prot_id) == 1:
            prot_id = prot_id[0]
            if prot_id in prot_acc_user:
                tax_id = taxids[prot_id]
                taxon = tax.scientific_name_for_taxid(tax_id)
                cds_acc = re.findall(r'^(.*?)\s',
                                     description)[0].split('_cds_')[0]
                prot_name = re.findall(r'\[protein=(.*?)\]', x)[0]
                gene_name = re.findall(r'\[gene=(.*?)\]', x)
                if len(gene_name) == 0:
                    gene_name = ''
                else:
                    gene_name = '__' + gene_name[0]
                x_seq = ''.join(x.split('\n')[1:])
                x_desc = ('>' + taxon + gene_name + '__' + prot_name +
                          '__QUERY__' + cds_acc + '__' + prot_id)
                x_desc = x_desc.replace(',', '')
                x_desc = x_desc.replace(';', '')
                x_desc = x_desc.replace(':', '')
                x_desc = x_desc.replace(' ', '_')
                x = x_desc + '\n' + x_seq
                cds_seqs_fasta_list.append(x)

    cds_seqs_fasta_text = '\n'.join(cds_seqs_fasta_list)

    with open(nt_prot_ncbi_file, 'w') as f:
        f.write(cds_seqs_fasta_text)

    return cds_seqs_fasta_text
