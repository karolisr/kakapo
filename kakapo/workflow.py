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

from collections import OrderedDict
from copy import deepcopy
from functools import partial
from os import remove as osremove
from os import stat as osstat
from os.path import basename
from os.path import commonprefix
from os.path import exists as ope
from os.path import join as opj
from os.path import sep as ops
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
from kakapo.bowtie2 import build_bt2_index, run_bowtie2_se, run_bowtie2_pe
from kakapo.config import PICKLE_PROTOCOL
from kakapo.data.start_codon_context import contexts as atg_contexts
from kakapo.ebi_domain_search import pfam_entry
from kakapo.ebi_domain_search import pfam_seqs
from kakapo.ebi_domain_search import prot_ids_for_tax_ids
from kakapo.ebi_iprscan5 import job_runner
from kakapo.ebi_iprscan5 import result_json
from kakapo.ebi_proteins import fasta_by_accession_list
from kakapo.entrez import accessions as accessions_ncbi
from kakapo.entrez import cds_acc_for_prot_acc
from kakapo.entrez import dnld_cds_nt_fasta as dnld_ncbi_cds_nt_fasta
from kakapo.entrez import dnld_seqs as dnld_ncbi_seqs
from kakapo.entrez import dnld_seqs_fasta_format
from kakapo.entrez import esearch
from kakapo.entrez import esummary as entrez_summary
from kakapo.entrez import sra_run_info
from kakapo.entrez import taxids_for_accs
from kakapo.gff3 import gff_from_kakapo_ips5_json_file
from kakapo.helpers import combine_text_files
from kakapo.helpers import keep_unique_lines_in_file
from kakapo.helpers import make_dir
from kakapo.helpers import plain_or_gzip
from kakapo.helpers import split_seq_defn_for_printing as split_seq_defn
from kakapo.helpers import splitext_gz
from kakapo.kraken import run_kraken_filters
from kakapo.orf import find_orf_for_blast_hit
from kakapo.py_v_diffs import StringIO
from kakapo.rcorrector import filter_unc_se, filter_unc_pe
from kakapo.rcorrector import run_rcorrector_se, run_rcorrector_pe
from kakapo.seq import reverse_complement, translate, untranslate
from kakapo.seqtk import seqtk_fq_to_fa, seqtk_extract_reads
from kakapo.shell import call
from kakapo.spades import run_spades_se, run_spades_pe
from kakapo.translation_tables import TranslationTable
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

    dir_cache_refseqs = opj(dir_cache, 'ref-seqs')
    make_dir(dir_cache_refseqs)

    dir_prj = opj(dir_out, '02-project-specific', prj_name)
    make_dir(dir_prj)

    dir_prj_logs = opj(dir_prj, '00-logs')
    make_dir(dir_prj_logs)

    dir_prj_queries = opj(dir_prj, '01-queries')
    make_dir(dir_prj_queries)

    dir_prj_blast_results_fa_trim = opj(dir_prj, '02-filtered-fa-blast-results')
    make_dir(dir_prj_blast_results_fa_trim)

    dir_prj_vsearch_results_fa_trim = opj(dir_prj,
                                          '03-filtered-fa-vsearch-results')
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
                'dir_cache_refseqs': dir_cache_refseqs,
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


def user_protein_accessions(prot_acc_user, dir_cache_prj, linfo=print):  # noqa
    if len(prot_acc_user) > 0:
        linfo('Reading user provided protein accessions')
        pickle_file = opj(dir_cache_prj, 'ncbi_prot_metadata_cache')
        acc_old = set()
        if ope(pickle_file):
            with open(pickle_file, 'rb') as f:
                pickled = pickle.load(f)
                acc_old = set([x['accessionversion'] for x in pickled])

        if acc_old == set(prot_acc_user):
            pa_info = pickled
        else:
            pa_info = entrez_summary(prot_acc_user, 'protein')

        prot_acc = []
        prot_info_to_print = []
        for pa in pa_info:
            acc = pa['accessionversion']
            prot_acc.append(acc)
            title = pa['title']
            title_split = title.split('[')
            organism = title_split[1].replace(']', '').strip().replace('_', ' ')
            title = title_split[0]
            title = title.lower().strip()
            title = title.replace('_', ' ').replace('-', ' ')
            title = title.replace(',', '')
            title = title[0].upper() + title[1:] + ' [' + organism + ']'
            prot_info_to_print.append((title, acc))

        prot_info_to_print = sorted(prot_info_to_print)
        for pi in prot_info_to_print:
            title = pi[0]
            acc = pi[1]
            if len(title) > 80:
                title = title[:77] + '...'
            linfo(acc + ': ' + title)

        with open(pickle_file, 'wb') as f:
            pickle.dump(pa_info, f, protocol=PICKLE_PROTOCOL)

        return prot_acc

    else:

        return prot_acc_user


def user_entrez_search(queries, dir_cache_prj, linfo=print):  # noqa
    accs = []
    if len(queries) != 0:
        linfo('Searching for protein sequences on NCBI')
        for q in queries:
            accs = accs + accessions_ncbi(esearch(term=q, db='protein'))

    return user_protein_accessions(accs, dir_cache_prj, linfo=linfo)


def dnld_prot_seqs(prot_acc_user, aa_prot_ncbi_file, linfo=print):  # noqa
    if len(prot_acc_user) != 0:
        acc_old = set()
        if ope(aa_prot_ncbi_file):
            _ = read_fasta(aa_prot_ncbi_file)
            acc_old = set([x.split('|')[0] for x in tuple(_.keys())])

        if acc_old == set(prot_acc_user):
            return prot_acc_user
        else:
            linfo('Downloading protein sequences from NCBI')
            _ = dnld_ncbi_seqs(prot_acc_user, 'protein')
            prot_acc_user_new = list()
            for rec in _:
                accession = rec['accession']
                version = rec['version']
                defn = rec['definition']
                organism = rec['organism']

                new_acc = accession + '.' + version
                prot_acc_user_new.append(new_acc)

                defn_new = defn.split('[' + organism + ']')[0]
                defn_new = defn_new.lower().strip()
                defn_new = defn_new.replace(' ', '_').replace('-', '_')
                defn_new = defn_new.replace(',', '')
                defn_new = defn_new[0].upper() + defn_new[1:]

                rec['definition'] = defn_new

            prot_acc_user = prot_acc_user_new
            write_fasta(_, aa_prot_ncbi_file)
    else:
        if ope(aa_prot_ncbi_file):
            osremove(aa_prot_ncbi_file)

    return prot_acc_user


def user_aa_fasta(user_queries, aa_prot_user_file, linfo=print):  # noqa
    _ = ''
    if len(user_queries) > 0:
        linfo('Reading user provided AA sequences')
        for ap in user_queries:
            linfo(ap)
            with open(ap, 'r') as f:
                _ = _ + f.read()
    if _ != '':
        with open(aa_prot_user_file, 'w') as f:
            f.write(standardize_fasta_text(_))


def combine_aa_fasta(fasta_files, aa_queries_file, linfo=print):  # noqa
    linfo('Combining all AA query sequences')
    _ = ''
    for fasta_file in fasta_files:
        if ope(fasta_file):
            with open(fasta_file, 'r') as f:
                _ = _ + f.read()

    if _ != '':
        with open(aa_queries_file, 'w') as f:
            f.write(_)
    else:
        linfo('No queries were provided. Exiting.')
        exit(0)


def filter_queries(aa_queries_file, min_query_length, max_query_length,
                   max_query_identity, vsearch, prot_acc_user, linfo=print): # noqa
    _ = ''
    with open(aa_queries_file, 'r') as f:
        _ = f.read()

    linfo('Filtering AA query sequences')
    linfo('min_query_length: ' + str(min_query_length))
    linfo('max_query_length: ' + str(max_query_length))

    _ = filter_fasta_text_by_length(_, min_query_length, max_query_length)

    tmp1 = aa_queries_file + '_temp1'
    tmp2 = aa_queries_file + '_temp2'
    tt = TranslationTable(1)
    parsed_fasta_1 = read_fasta(StringIO(_))
    for dfn in parsed_fasta_1:
        parsed_fasta_1[dfn] = untranslate(parsed_fasta_1[dfn], tt.table_inv)
    write_fasta(parsed_fasta_1, tmp1)
    run_cluster_fast(vsearch, max_query_identity, tmp1, tmp2)
    parsed_fasta_2 = read_fasta(tmp2)
    prot_acc_user_new = list()
    for dfn in parsed_fasta_2:
        parsed_fasta_2[dfn] = translate(parsed_fasta_2[dfn], tt.table,
                                        tt.start_codons)
        acc = dfn.split('|')[0]
        if acc in prot_acc_user:
            prot_acc_user_new.append(acc)

    write_fasta(parsed_fasta_2, aa_queries_file)

    osremove(tmp1)
    osremove(tmp2)

    return prot_acc_user_new


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

    nuclear = ''
    for nuc in order:
        if nuc[1] == 'nuclear':
            nuclear = nuc[0]
            break

    for se in se_fastq_files:
        fq_path = se_fastq_files[se]['trim_path_fq']
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


def run_bt2_fq(se_fastq_files, pe_fastq_files, dir_fq_filter_data,
               bowtie2, bowtie2_build, threads, dir_temp, filter_dir, dbs,
               fpatt, taxonomy, dir_cache_refseqs, linfo=print):  # noqa

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
                         dir_blast_results_fa_trim, fpatt, threads,
                         seqtk, vsearch, linfo=print): # noqa
    ident = 0.85

    for se in se_fastq_files:
        dir_results = opj(dir_blast_results_fa_trim, se)
        blast_db_path = se_fastq_files[se]['blast_db_path']
        fq_path = se_fastq_files[se]['filter_path_fq']
        out_f = opj(dir_results, se + '.txt')
        out_f_fastq = out_f.replace('.txt', '.fastq')
        out_f_fasta = out_f.replace('.txt', '.fasta')
        se_fastq_files[se]['blast_results_path'] = out_f_fasta
        genetic_code = se_fastq_files[se]['gc_id']

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
                      max_hsps=blast_1_max_hsps,
                      qcov_hsp_perc=blast_1_qcov_hsp_perc,
                      best_hit_overhang=blast_1_best_hit_overhang,
                      best_hit_score_edge=blast_1_best_hit_score_edge,
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
        fq_paths = pe_fastq_files[pe]['filter_path_fq']
        out_fs = [x.replace('@D@', dir_results) for x in fpatt]
        out_fs = [x.replace('@N@', pe) for x in out_fs]
        out_fs_fastq = [x.replace('.txt', '.fastq') for x in out_fs]
        out_fs_fasta = [x.replace('.txt', '.fasta') for x in out_fs]
        out_f_fasta = opj(dir_results, pe + '.fasta')
        pe_fastq_files[pe]['blast_results_path'] = out_f_fasta
        genetic_code = pe_fastq_files[pe]['gc_id']

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
                          max_hsps=blast_1_max_hsps,
                          qcov_hsp_perc=blast_1_qcov_hsp_perc,
                          best_hit_overhang=blast_1_best_hit_overhang,
                          best_hit_score_edge=blast_1_best_hit_score_edge,
                          max_target_seqs=blast_1_max_target_seqs,
                          db_genetic_code=genetic_code,
                          out_cols=BLST_RES_COLS_1)

                linfo('Extracting unique BLAST hits using Seqtk')

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
        fq_path = se_fastq_files[se]['filter_path_fq']
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
        fq_paths = pe_fastq_files[pe]['filter_path_fq']
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
                              blast_2_max_hsps, blast_2_qcov_hsp_perc,
                              blast_2_best_hit_overhang,
                              blast_2_best_hit_score_edge,
                              blast_2_max_target_seqs, threads, dir_cache_prj,
                              dir_prj_ips, linfo=print):  # noqa
    if len(assemblies) > 0:
        linfo('Running BLAST on assemblies')
    else:
        linfo('There are no assemblies. Nothing to do, stopping.')
        exit(0)

    cache_file = opj(dir_cache_prj, 'blast_2_settings_cache')

    pickled = dict()
    settings = {'blast_2_evalue': blast_2_evalue,
                'blast_2_max_hsps': blast_2_max_hsps,
                'blast_2_qcov_hsp_perc': blast_2_qcov_hsp_perc,
                'blast_2_best_hit_overhang': blast_2_best_hit_overhang,
                'blast_2_best_hit_score_edge': blast_2_best_hit_score_edge,
                'blast_2_max_target_seqs': blast_2_max_target_seqs,
                'queries': read_fasta(aa_queries_file)}

    linfo('evalue: ' + str(blast_2_evalue))
    linfo('max_hsps: ' + str(blast_2_max_hsps))
    linfo('qcov_hsp_perc: ' + str(blast_2_qcov_hsp_perc))
    linfo('best_hit_overhang: ' + str(blast_2_best_hit_overhang))
    linfo('best_hit_score_edge: ' + str(blast_2_best_hit_score_edge))
    linfo('max_target_seqs: ' + str(blast_2_max_target_seqs))

    for a in assemblies:
        assmbl_name = a['name']
        assmbl_blast_db_path = a['blast_db_path']
        assmbl_genetic_code = a['gc_id']

        ips_json_dump_path = opj(dir_prj_ips, assmbl_name + '_ann_ips.json')

        _ = opj(dir_prj_assmbl_blast_results, assmbl_name + '.tsv')

        if ope(_) and ope(cache_file):
            with open(cache_file, 'rb') as f:
                pickled = pickle.load(f)

        if ope(_) and pickled == settings:
            linfo('The provided BLAST settings and query sequences did not ' +
                  'change since the previous run. BLAST results for the ' +
                  'assembly "' + assmbl_name + '" already exist')

        else:
            linfo('Running tblastn on: ' + assmbl_name)

            if ope(ips_json_dump_path):
                osremove(ips_json_dump_path)

            run_blast(exec_file=tblastn,
                      task='tblastn',
                      threads=threads,
                      db_path=assmbl_blast_db_path,
                      queries_file=aa_queries_file,
                      out_file=_,
                      evalue=blast_2_evalue,
                      max_hsps=blast_2_max_hsps,
                      qcov_hsp_perc=blast_2_qcov_hsp_perc,
                      best_hit_overhang=blast_2_best_hit_overhang,
                      best_hit_score_edge=blast_2_best_hit_score_edge,
                      max_target_seqs=blast_2_max_target_seqs,
                      db_genetic_code=assmbl_genetic_code,
                      out_cols=BLST_RES_COLS_2)

        a['blast_hits_aa'] = parse_blast_results_file(_, BLST_RES_COLS_2)

    with open(cache_file, 'wb') as f:
        pickle.dump(settings, f, protocol=PICKLE_PROTOCOL)


def find_orfs_translate(assemblies, dir_prj_transcripts, seqtk,
                        dir_temp, prepend_assmbl, min_target_orf_len,
                        max_target_orf_len, allow_non_aug, allow_no_strt_cod,
                        allow_no_stop_cod, tax, tax_group, tax_ids_user,
                        min_overlap, linfo=print):  # noqa
    if len(assemblies) > 0:
        linfo('Analyzing BLAST hits for assemblies')

    for a in assemblies:

        assmbl_name = a['name']
        tax_id = a['tax_id']

        parsed_hits = a['blast_hits_aa']
        a_path = a['path']
        gc_tt = a['gc_tt']

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
            _ = f.read()

        if _.strip() == '':
            continue

        linfo(assmbl_name)

        _ = trim_desc_to_first_space_in_fasta_text(_)

        parsed_fasta = read_fasta(StringIO(_))
        ######################################################################

        all_kakapo_results = {}
        json_dump_file_path = opj(dir_prj_transcripts, assmbl_name +
                                  '_ann_kakapo.json')

        if len(collated) > 0:
            print('\n' + '-' * 90 + '\n')

        for hit in collated:

            target_name = hit['sseqid']
            target_seq = parsed_fasta[target_name]
            query_name = hit['qseqid']
            hit_evalue = hit['evalue']

            # Prepend assembly name to the sequence name:
            if prepend_assmbl is True:
                target_name = assmbl_name + '__' + target_name
                # Also prepend taxonomic info to the sequence name:
                if tax_id is not None:
                    fm = tax.higher_rank_for_taxid(tax_id, rank='family')
                    if fm is not None:
                        target_name = fm + '__' + target_name
                    # cn = tax.genbank_common_name_for_taxid(tax_id)
                    # if cn is not None:
                    #     cn = cn.lower()
                    #     cn = cn.replace(' ', '_')
                    #     target_name = cn + '__' + target_name

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

            orf_log_str = target_name.center(90) + '\n'
            orf_log_str += query_name.center(90) + '\n\n'

            orf_log_str += ('grade'.rjust(6) + 'ovrlp'.rjust(7) +
                            'cntx'.rjust(6) + 'length'.center(9) +
                            'cntx_l'.rjust(7) + 'cntx_r'.rjust(15) + '\n')

            orf = find_orf_for_blast_hit(
                seq=target_seq,
                frame=hit_frame,
                hit_start=hit_start,
                hit_end=hit_end,
                stop_codons=stop_codons,
                start_codons=start_codons,
                context_l=cntx_l,
                context_r=cntx_r,
                min_overlap=min_overlap)

            orf_log_str += orf[2]

            rev_comp_def_str = ''
            if hit_frame > 0:
                ann_hit_b = hit_start
                ann_hit_e = hit_end
            else:
                target_seq = reverse_complement(target_seq)
                ann_hit_b = len(target_seq) - hit_start
                ann_hit_e = len(target_seq) - hit_end
                rev_comp_def_str = '; RevComp'

            target_def = target_name + ' ' + query_name + rev_comp_def_str

            a['annotations'][target_name] = {}

            good_orf = orf[0]
            bad_orfs = orf[1]

            if good_orf is not None:

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
                invalid_orf_reason = ''

                if allow_non_aug is False and \
                        orf_seq[0:3] != 'ATG':
                    valid_orf = False
                    invalid_orf_reason = 'Start codon is not ATG.'

                elif allow_no_strt_cod is False and \
                        orf_seq[0:3] not in start_codons:
                    valid_orf = False
                    invalid_orf_reason = 'No start codon.'

                elif allow_no_stop_cod is False and \
                        orf_seq[-3:] not in stop_codons:
                    valid_orf = False
                    invalid_orf_reason = 'No stop codon.'

                elif len(orf_seq) < min_target_orf_len:
                    valid_orf = False
                    invalid_orf_reason = 'ORF is not long enough.'

                elif len(orf_seq) > max_target_orf_len:
                    valid_orf = False
                    invalid_orf_reason = 'ORF is too long.'

                ##############################################################

                if valid_orf is True:

                    orf_log_str += '\n' + 'VALID'.center(90) + '\n'

                    a['annotations'][target_name]['orf_begin'] = ann_orf_b
                    a['annotations'][target_name]['orf_end'] = ann_orf_e
                    a['annotations'][target_name]['orf_grade'] = good_orf[3]

                    transcripts_nt_orf[target_def] = orf_seq

                    transl_seq = translate(orf_seq,
                                           gc_tt.table_ambiguous,
                                           start_codons)

                    transcripts_aa_orf[target_def] = transl_seq[:-1]

                else:
                    msg = 'INVALID: ' + invalid_orf_reason
                    orf_log_str += '\n' + msg.center(90) + '\n'
                    bad_orfs.append(good_orf)

            else:
                orf_log_str += '\n' + 'INVALID'.center(90) + '\n'

            orf_log_str += '\n' + '-' * 90 + '\n'
            print(orf_log_str)

            if len(bad_orfs) > 0:
                a['annotations'][target_name]['orfs_bad'] = list()
                orfs_bad_list = a['annotations'][target_name]['orfs_bad']

            for bad_orf in bad_orfs:

                bad_orf_frame = bad_orf[2]

                if bad_orf_frame > 0:
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
                orf_bad_dict['orf_frame'] = abs(bad_orf_frame)
                orf_bad_dict['orf_grade'] = bad_orf[3]

                orfs_bad_list.append(orf_bad_dict)

            transcripts_nt[target_def] = target_seq

            a['annotations'][target_name]['query_name'] = query_name
            a['annotations'][target_name]['evalue'] = hit_evalue
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
    delay = 0.25

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
        seqs = OrderedDict(sorted(seqs.items(),
                                  key=lambda x: x[0].split(' ')[1],
                                  reverse=True))

        _ = opj(dir_cache_prj, assmbl_name + '_ips_jobs')

        if ope(_):
            with open(_, 'rb') as f:
                jobs = pickle.load(f)

        else:
            jobs = job_runner(email=email, dir_cache=dir_cache_prj,
                              seqs=seqs, logger=linfo)

            with open(_, 'wb') as f:
                pickle.dump(jobs, f, protocol=PICKLE_PROTOCOL)

        print()
        linfo('Downloading InterProScan results for transcripts in ' +
              assmbl_name)
        print()

        all_ips_results = {}

        # Nicer printing
        max_title_a_len = 2 + max([len(split_seq_defn(x)[0]) for x in list(jobs['finished'].keys())])
        max_title_b_len = 2 + max([len(split_seq_defn(x)[1]) for x in list(jobs['finished'].keys())])

        for i, job in enumerate(jobs['finished']):

            job_id = jobs['finished'][job]

            titles_ab = split_seq_defn(job)
            title_a = titles_ab[0]
            title_b = titles_ab[1]

            progress = round(((i + 1) / len(jobs['finished'])) * 100)
            progress_str = '{:3d}'.format(progress) + '%'

            msg = (' ' * 12 +
                   title_a.ljust(max_title_a_len) +
                   title_b.ljust(max_title_b_len) +
                   progress_str + ' ' + job_id)

            linfo(msg)

            sleep(delay)

            ips_json = result_json(job_id)
            # ips_version = ips_json['interproscan-version']
            ips_json = ips_json['results']

            # These fields are set to 'EMBOSS_001' by default
            # Delete them
            del ips_json[0]['xref']

            job_no_def = job.split(' ')[0]

            # kakapo annotations
            ips_json[0]['kakapo_annotations'] = a['annotations'][job_no_def]

            all_ips_results[job_no_def] = ips_json

        print()

        with open(json_dump_file_path, 'w') as f:
            json.dump(all_ips_results, f, sort_keys=True, indent=4)

        # Removes cached jobs file.
        osremove(_)


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


def dnld_cds_for_ncbi_prot_acc(prot_acc_user, prot_cds_ncbi_file, tax,
                               dir_cache_prj, linfo=print):  # noqa

    pickle_file = opj(dir_cache_prj, 'ncbi_prot_cds_cache')
    acc_old = set()
    if ope(pickle_file):
        with open(pickle_file, 'rb') as f:
            pickled = pickle.load(f)
            acc_old = set(pickled[0].keys())

    if acc_old == set(prot_acc_user):
        cds_fasta = pickled[1]
        taxids = pickled[2]
    else:
        linfo('Downloading CDS for the dereplicated set of the user-provided '
              'NCBI protein accessions')
        cds_acc_dict = cds_acc_for_prot_acc(prot_acc_user)
        cds_accessions = []
        for prot_acc in cds_acc_dict:
            cds_acc = cds_acc_dict[prot_acc]
            cds_accessions.append(cds_acc)
        cds_accessions = sorted(set(cds_accessions))
        cds_fasta = dnld_ncbi_cds_nt_fasta(cds_accessions)
        taxids = taxids_for_accs(prot_acc_user, 'protein')
        with open(pickle_file, 'wb') as f:
            pickle.dump((cds_acc_dict, cds_fasta, taxids), f,
                        protocol=PICKLE_PROTOCOL)

    prot_ids_used = []
    cds_seqs_fasta_list = []
    for rec in cds_fasta:
        description = rec.split('|')[1]
        prot_id = re.findall(r'\[protein_id=(.*?)\]', rec)

        if len(prot_id) == 1:

            prot_id = prot_id[0]

            if prot_id in prot_acc_user:
                if prot_id in prot_ids_used:
                    continue

                prot_ids_used.append(prot_id)

                taxid = taxids[prot_id]
                taxon = tax.scientific_name_for_taxid(taxid)
                seq = cds_fasta[rec]
                cds_acc = re.findall(r'^(.*?)\s',
                                     description)[0].split('_cds_')[0]
                prot_name = re.findall(r'\[protein=(.*?)\]', rec)[0]

                prot_name = prot_name.lower().strip()
                prot_name = prot_name.replace(' ', '_').replace('-', '_')
                prot_name = prot_name.replace(',', '')
                prot_name = prot_name[0].upper() + prot_name[1:]

                defn = prot_name + '__' + prot_id + '__QUERY'
                defn = defn + ' ' + prot_name.replace('_', ' ')
                defn = defn + '; ' + taxon + ', ' + str(taxid)
                defn = defn + '; ' + prot_id + '; ' + cds_acc

                rec_new = '>' + defn + '\n' + seq
                cds_seqs_fasta_list.append(rec_new)

    cds_seqs_fasta_text = '\n'.join(cds_seqs_fasta_list)

    with open(prot_cds_ncbi_file, 'w') as f:
        f.write(cds_seqs_fasta_text)
