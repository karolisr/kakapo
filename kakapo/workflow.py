#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""kakapo workflow"""

import pickle

from os.path import exists as ope
from os.path import join as opj
from os.path import splitext
from os.path import basename
from os.path import commonprefix
from os import remove as osremove

import re

from kakapo.bioio import dnld_ncbi_seqs
from kakapo.bioio import entrez_summary
from kakapo.bioio import filter_fasta_text_by_length
from kakapo.bioio import sra_info
from kakapo.bioio import standardize_fasta_text
from kakapo.bioio import write_fasta_file
from kakapo.config import PICKLE_PROTOCOL
from kakapo.ebi_domain_search import pfam_entry
from kakapo.ebi_domain_search import pfam_seqs
from kakapo.ebi_domain_search import prot_ids_for_tax_ids
from kakapo.ebi_proteins import fasta_by_accession_list
from kakapo.helpers import make_dir
from kakapo.shell import call
from kakapo.trimmomatic import trimmomatic_se, trimmomatic_pe


def prepare_output_directories(dir_out, prj_name):  # noqa

    # -- ToDo: Lock cache files in case of parallel execution ------------

    dir_temp = opj(dir_out, '00-b-temp')
    make_dir(dir_temp)

    dir_cache = opj(dir_out, '00-a-cache')
    make_dir(dir_cache)

    dir_cache_pfam_acc = opj(dir_cache, 'pfam-uniprot-accessions')
    make_dir(dir_cache_pfam_acc)

    dir_cache_prj = opj(dir_cache, 'projects', prj_name)
    make_dir(dir_cache_prj)

    dir_prj = opj(dir_out, '01-projects', prj_name)
    make_dir(dir_prj)

    dir_prj_queries = opj(dir_prj, '01-queries')
    make_dir(dir_prj_queries)

    dir_fq_data = opj(dir_out, '11-sra-fq-data')
    make_dir(dir_fq_data)

    dir_fq_trim_data = opj(dir_out, '12-trimmed-fq-data')
    make_dir(dir_fq_trim_data)

    dir_fa_trim_data = opj(dir_out, '13-trimmed-fa-data')
    make_dir(dir_fa_trim_data)

    ret_dict = {'dir_temp': dir_temp,
                'dir_cache': dir_cache,
                'dir_cache_pfam_acc': dir_cache_pfam_acc,
                'dir_cache_prj': dir_cache_prj,
                'dir_prj': dir_prj,
                'dir_prj_queries': dir_prj_queries,
                'dir_fq_data': dir_fq_data,
                'dir_fq_trim_data': dir_fq_trim_data,
                'dir_fa_trim_data': dir_fa_trim_data}

    return ret_dict

def descending_tax_ids(tax_ids_user, taxonomy):  # noqa
    if len(tax_ids_user) > 0:
        print('Resolving descending nodes for:')
    tax_ids = []
    for tx in tax_ids_user:
        tx_name = taxonomy.scientific_name_for_taxid(taxid=tx)
        print('\t\t' + str(tx) + ': ' + tx_name)
        tax_ids_for_tx = taxonomy.all_descending_taxids(taxid=tx)
        tax_ids = tax_ids + tax_ids_for_tx
    tax_ids = list(set(tax_ids))
    # ToDo: Fix Taxonomy class, so it returns integers by default
    tax_ids = [int(x) for x in tax_ids]

    return tax_ids


def pfam_uniprot_accessions(pfam_acc, tax_ids, dir_cache_pfam_acc):  # noqa
    if len(pfam_acc) > 0:
        print('\nDownloading UniProt accessions for Pfam families:')
    pfam_seqs_list = []
    for pa in pfam_acc:
        pfam_id = pfam_entry(pa)[0]['id']
        print('\t\t' + pa + ': ' + pfam_id)
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
        print('Reading user provided protein accessions:')
        pa_info = entrez_summary(prot_acc_user, 'protein')
        prot_acc = []
        for pa in pa_info:
            title = pa['Title']
            title = title[0].upper() + title[1:]
            acc = pa['AccessionVersion']
            prot_acc.append(acc)
            print('\t\t' + acc + ': ' + title)
        print()

        return prot_acc


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
        with open(fasta_file, 'r') as f:
            __ = __ + f.read()

    if __ != '':
        with open(aa_queries_file, 'w') as f:
            f.write(__)


def filter_queries(aa_queries_file, min_query_length, max_query_length): # noqa
    __ = ''
    with open(aa_queries_file, 'r') as f:
        __ = f.read()

    print('Filtering AA query sequences:')
    print('\t\tmin_query_length: ' + str(min_query_length))
    print('\t\tmax_query_length: ' + str(max_query_length))

    __ = filter_fasta_text_by_length(__, min_query_length, max_query_length)

    with open(aa_queries_file, 'w') as f:
        f.write(__)


def dnld_sra_info(sras, dir_cache_prj):  # noqa
    if len(sras) > 0:
        print('\nDownloading SRA run information:\n')

    sra_runs_info = {}

    __ = opj(dir_cache_prj, 'sra_runs_info_cache')

    if ope(__):
        with open(__, 'rb') as f:
            sra_runs_info = pickle.load(f)

    for sra in sras:
        if sra not in sra_runs_info:
            parsed = sra_info(sra)[0]
            sra_lib_layout = parsed['LibraryLayout'].lower()
            sra_lib_source = parsed['LibrarySource'].lower()
            sra_lib_strategy = parsed['LibraryStrategy']
            sra_seq_platform = parsed['Platform'].lower().capitalize()
            sra_seq_platform_model = parsed['Model']
            sra_species = parsed['ScientificName']
            sra_taxid = parsed['TaxID']

            sra_run_info = {}

            sra_run_info['sra_lib_layout'] = sra_lib_layout
            sra_run_info['sra_lib_source'] = sra_lib_source
            sra_run_info['sra_lib_strategy'] = sra_lib_strategy
            sra_run_info['sra_seq_platform'] = sra_seq_platform
            sra_run_info['sra_seq_platform_model'] = sra_seq_platform_model
            sra_run_info['sra_species'] = sra_species
            sra_run_info['sra_taxid'] = sra_taxid

            sample_base_name = (sra_species.replace(' ', '_') + '_' +
                                sra_taxid + '_' + sra)

            sra_run_info['sample_base_name'] = sample_base_name

        else:
            sra_run_info = sra_runs_info[sra]

        if sra_run_info['sra_lib_source'] != 'transcriptomic':
            sra_info_str = (
                '{sra}: the SRA library source type "{ltype}" '
                'is not supported.').format(
                sra=sra, ltype=sra_run_info['sra_lib_source'])

        else:
            sra_info_str = ('\tSRA run {sra}\n\t{source} '
                            '{strategy} {layout}-end library.\n'
                            '\tSourced from {species} '
                            '(Tax ID: {txid}).\n'
                            '\tSequenced using {platform} platform on '
                            '{model}.\n').format(
                                sra=sra,
                                source=sra_run_info['sra_lib_source'].title(),
                                strategy=sra_run_info['sra_lib_strategy'],
                                layout=sra_run_info['sra_lib_layout'],
                                platform=sra_run_info['sra_seq_platform'],
                                model=sra_run_info['sra_seq_platform_model'],
                                species=sra_run_info['sra_species'],
                                txid=sra_run_info['sra_taxid'])

            sra_runs_info[sra] = sra_run_info

        print(sra_info_str)

    with open(__, 'wb') as f:
        pickle.dump(sra_runs_info, f, protocol=PICKLE_PROTOCOL)

    return sra_runs_info


def dnld_sra_fastq_files(sras, sra_runs_info, dir_fq_data, fasterq_dump,
                         threads, dir_temp): # noqa

    se_fastq_files = {}
    pe_fastq_files = {}

    for sra in sras:
        sra_run_info = sra_runs_info[sra]
        sra_lib_layout = sra_run_info['sra_lib_layout']
        sample_base_name = sra_run_info['sample_base_name']

        sra_dnld_needed = False

        if sra_lib_layout == 'single':
            se_file = opj(dir_fq_data, sra + '.fastq')
            se_fastq_files[sample_base_name] = {'path': se_file}
            if not ope(se_file):
                sra_dnld_needed = True

        elif sra_lib_layout == 'paired':
            pe_file_1 = opj(dir_fq_data, sra + '_1.fastq')
            pe_file_2 = opj(dir_fq_data, sra + '_2.fastq')
            pe_fastq_files[sample_base_name] = {'path': [pe_file_1, pe_file_2]}
            if not ope(pe_file_1) or not ope(pe_file_2):
                sra_dnld_needed = True

        if sra_dnld_needed:
            print('Downloading FASTQ reads for the SRA accession ' + sra)
            cmd = [fasterq_dump,
                   '--threads', str(threads),
                   '--split-files',
                   '--outdir', dir_fq_data,
                   '--temp', dir_temp, sra]
            call(cmd)

        else:
            print('\tFASTQ reads for the SRA run ' + sra +
                  ' are available locally.')

    return se_fastq_files, pe_fastq_files


def user_fastq_files(fq_se, fq_pe): # noqa

    if len(fq_se) > 0 or len(fq_pe) > 0:
        print('\nPreparing user provided FASTQ files:\n')

    se_fastq_files = {}
    pe_fastq_files = {}

    for se in fq_se:
        sample_base_name = splitext(basename(se))[0]
        se_fastq_files[sample_base_name] = {'path': se}
        print('\t' + sample_base_name + ':\n\t\t' + se)
        print()

    for pe in fq_pe:
        sample_base_name = basename(commonprefix(pe))
        sample_base_name = sample_base_name.rstrip('_- R')
        pe_fastq_files[sample_base_name] = {'path': pe}
        print('\t' + sample_base_name + ':\n\t\t' + pe[0] + '\n\t\t' + pe[1])
        print()

    return se_fastq_files, pe_fastq_files


def min_accept_read_len(se_fastq_files, pe_fastq_files, dir_temp, vsearch): # noqa

    print('Calculating minimum acceptable read length:\n')

    queue = []

    for se in se_fastq_files:
        fq_path = se_fastq_files[se]['path']
        stats_file = opj(dir_temp, se + '_stats.txt')
        queue.append([se, fq_path, stats_file, 'se'])

    for pe in pe_fastq_files:
        fq_path = pe_fastq_files[pe]['path'][0]
        stats_file = opj(dir_temp, pe + '_stats.txt')
        queue.append([pe, fq_path, stats_file, 'pe'])

    for x in queue:
        cmd = [vsearch, '--fastq_stats', x[1], '--log', x[2]]
        call(cmd)

        with open(x[2]) as f:
            stats = f.read()

        osremove(x[2])

        ml = re.findall(r'>=\s+(\d+)', stats)

        if len(ml) != 0:
            ml = int(ml[0]) // 2
            print('\t' + x[0] + ': ' + str(ml))
        else:
            ml = None
            print('\t' + x[0] +
                  ': could not be determined')

        if x[3] == 'se':
            se_fastq_files[x[0]]['min_acc_len'] = ml
        elif x[3] == 'pe':
            pe_fastq_files[x[0]]['min_acc_len'] = ml

    print()


def run_trimmomatic(se_fastq_files, pe_fastq_files, dir_fq_trim_data,
                    trimmomatic, adapters, fpatt, threads): # noqa

    for se in se_fastq_files:
        dir_fq_trim_data_sample = opj(dir_fq_trim_data, se)
        fq_path = se_fastq_files[se]['path']
        min_acc_len = se_fastq_files[se]['min_acc_len']
        stats_f = opj(dir_fq_trim_data_sample, se + '.txt')
        out_f = opj(dir_fq_trim_data_sample, se + '.fastq')
        se_fastq_files[se]['trim_path'] = out_f

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
        out_fs = [x.replace('xDIRx', dir_fq_trim_data_sample) for x in fpatt]
        out_fs = [x.replace('xBASENAMEx', pe) for x in out_fs]
        pe_fastq_files[pe]['trim_path'] = out_fs

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
