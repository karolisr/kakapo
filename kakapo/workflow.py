#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""kakapo workflow"""

import pickle

from os.path import exists as ope
from os.path import join as opj
from os import remove as osremove

from kakapo.helpers import make_dir
from kakapo.ebi_domain_search import pfam_seqs
from kakapo.ebi_domain_search import pfam_entry
from kakapo.ebi_domain_search import prot_ids_for_tax_ids
from kakapo.bioio import entrez_summary
from kakapo.bioio import dnld_ncbi_seqs
from kakapo.bioio import write_fasta_file
from kakapo.ebi_proteins import fasta_by_accession_list
from kakapo.bioio import standardize_fasta_text
from kakapo.bioio import filter_fasta_text_by_length
from kakapo.bioio import sra_info

from kakapo.config import PICKLE_PROTOCOL

def prepare_output_directories(dir_out, prj_name):  # noqa

    # -- ToDo: Lock cache files in case of parallel execution ------------
    dir_cache = opj(dir_out, '00-cache')
    make_dir(dir_cache)

    dir_cache_pfam_acc = opj(dir_cache, 'pfam-uniprot-accessions')
    make_dir(dir_cache_pfam_acc)

    dir_cache_prj = opj(dir_cache, 'projects', prj_name)
    make_dir(dir_cache_prj)

    dir_prj = opj(dir_out, '01-projects', prj_name)
    make_dir(dir_prj)

    dir_prj_queries = opj(dir_prj, '01-queries')
    make_dir(dir_prj_queries)

    ret_dict = {'dir_cache': dir_cache,
                'dir_cache_pfam_acc': dir_cache_pfam_acc,
                'dir_cache_prj': dir_cache_prj,
                'dir_prj': dir_prj,
                'dir_prj_queries': dir_prj_queries}

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

    print()


def dnld_sra_info(sras, dir_cache_prj):  # noqa
    if len(sras) > 0:
        print('Downloading SRA run information:\n')

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
