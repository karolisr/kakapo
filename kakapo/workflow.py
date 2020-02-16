#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Kakapo Workflow"""

import json
import pickle
import re

from collections import OrderedDict
from functools import partial
from io import StringIO
from os import remove as osremove
from os import stat as osstat
from os.path import exists as ope
from os.path import join as opj
from sys import exit
from time import sleep

from kakapo.bioio import read_fasta
from kakapo.bioio import trim_desc_to_first_space_in_fasta_text
from kakapo.bioio import write_fasta
from kakapo.blast import BLST_RES_COLS_2
from kakapo.blast import collate_blast_results
from kakapo.blast import make_blast_db, run_blast
from kakapo.blast import parse_blast_results_file
from kakapo.config import PICKLE_PROTOCOL
from kakapo.data.start_codon_context import contexts as atg_contexts
from kakapo.ebi_iprscan5 import job_runner
from kakapo.ebi_iprscan5 import result_json
from kakapo.entrez import cds_acc_for_prot_acc
from kakapo.entrez import dnld_cds_nt_fasta as dnld_ncbi_cds_nt_fasta
from kakapo.entrez import taxids_for_accs
from kakapo.gff3 import gff_from_kakapo_ips5_json_file
from kakapo.helpers import combine_text_files
from kakapo.helpers import make_dir
from kakapo.helpers import split_seq_defn_for_printing as split_seq_defn
from kakapo.orf import find_orf_for_blast_hit
from kakapo.seq import reverse_complement, translate
from kakapo.seqtk import seqtk_extract_reads
from kakapo.spades import run_spades_se, run_spades_pe


def prepare_output_directories(dir_out, prj_name):  # noqa
    # TODO: Lock cache files in case of parallel execution -------------------
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

    dir_prj_blast_results_fa_trim = opj(dir_prj,
                                        '02-filtered-fa-blast-results')
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


def descending_tax_ids(tax_ids_user, taxonomy, linfo=print):  # noqa
    shared = taxonomy.shared_taxid_for_taxids(tax_ids_user)
    if shared is None:
        return None
    tax_ids = taxonomy.all_descending_taxids(taxid=shared)
    if tax_ids is None:
        tax_ids = [shared]
    tax_ids = [int(x) for x in tax_ids]
    return tax_ids


def run_spades(se_fastq_files, pe_fastq_files, dir_spades_assemblies,
               spades, dir_temp, ss, threads, ram, linfo=print):  # noqa

    if len(se_fastq_files) > 0 or len(pe_fastq_files) > 0:
        if spades is None:
            linfo('spades is not available. ' +
                  'Cannot continue. Exiting.')
            exit(0)

    for se in se_fastq_files:
        dir_results = opj(dir_spades_assemblies, se + '__' + ss)
        fq_path = se_fastq_files[se]['vsearch_results_path' + '__' + ss]
        se_fastq_files[se]['spades_assembly' + '__' + ss] = None

        if ope(dir_results):
            linfo('SPAdes results for sample ' + se + ' already exist [' + ss + ']')
        else:
            make_dir(dir_results)
            linfo('Running SPAdes on: ' + se + ' [' + ss + ']')
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
            linfo('SPAdes produced ' + str(count) + tr_str + ' [' + ss + ']')
            se_fastq_files[se]['spades_assembly' + '__' + ss] = assmbl_path
        else:
            linfo('SPAdes produced no transcripts [' + ss + ']')

    for pe in pe_fastq_files:
        dir_results = opj(dir_spades_assemblies, pe + '__' + ss)
        fq_paths = pe_fastq_files[pe]['vsearch_results_path' + '__' + ss]
        pe_fastq_files[pe]['spades_assembly' + '__' + ss] = None

        if ope(dir_results):
            linfo('SPAdes results for sample ' + pe + ' already exist [' + ss + ']')
        else:
            make_dir(dir_results)
            linfo('Running SPAdes on: ' + pe + ' [' + ss + ']')

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
            linfo('SPAdes produced ' + str(count) + tr_str + ' [' + ss + ']')
            pe_fastq_files[pe]['spades_assembly' + '__' + ss] = assmbl_path
        else:
            linfo('SPAdes produced no transcripts [' + ss + ']')


def makeblastdb_assemblies(assemblies, dir_prj_blast_assmbl, makeblastdb,
                           linfo=print):  # noqa
    if len(assemblies) > 0:
        linfo('Building BLAST databases for assemblies')
        if makeblastdb is None:
            linfo('makeblastdb is not available. ' +
                  'Cannot continue. Exiting.')
            exit(0)
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


def run_tblastn_on_assemblies(ss, assemblies, aa_queries_file, tblastn,
                              dir_prj_assmbl_blast_results, blast_2_evalue,
                              blast_2_max_hsps, blast_2_qcov_hsp_perc,
                              blast_2_best_hit_overhang,
                              blast_2_best_hit_score_edge,
                              blast_2_max_target_seqs, threads, dir_cache_prj,
                              dir_prj_ips, linfo=print):  # noqa

    if len(assemblies) > 0:
        linfo('Running BLAST on assemblies [' + ss + ']')
        if tblastn is None:
            linfo('tblastn is not available. ' +
                  'Cannot continue. Exiting.')
            exit(0)
    else:
        linfo('There are no assemblies. Nothing to do, stopping.')
        exit(0)

    cache_file = opj(dir_cache_prj, 'blast_2_settings_cache__' + ss)

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
        assmbl_src = a['src']
        assmbl_name = a['name']

        if assmbl_src != 'user_fasta':
            if assmbl_name.endswith('__' + ss):
                assmbl_name = assmbl_name.replace('__' + ss, '')
            else:
                continue

        assmbl_blast_db_path = a['blast_db_path']
        assmbl_genetic_code = a['gc_id']

        ips_json_dump_path = opj(dir_prj_ips, assmbl_name + '_ann_ips__' + ss +
                                 '.json')

        _ = opj(dir_prj_assmbl_blast_results, assmbl_name + '__' + ss + '.tsv')

        if ope(_) and ope(cache_file):
            with open(cache_file, 'rb') as f:
                pickled = pickle.load(f)

        if ope(_) and pickled == settings:
            linfo('The provided BLAST settings and query sequences did not ' +
                  'change since the previous run. BLAST results for the ' +
                  'assembly "' + assmbl_name + '" already exist')

        else:
            linfo('Running tblastn on: ' + assmbl_name + ' [' + ss + ']')

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

        a['blast_hits_aa__' + ss] = parse_blast_results_file(_, BLST_RES_COLS_2)

    with open(cache_file, 'wb') as f:
        pickle.dump(settings, f, protocol=PICKLE_PROTOCOL)


def find_orfs_translate(ss, assemblies, dir_prj_transcripts, seqtk,
                        dir_temp, prepend_assmbl, min_target_orf_len,
                        max_target_orf_len, allow_non_aug, allow_no_strt_cod,
                        allow_no_stop_cod, tax, tax_group, tax_ids_user,
                        min_overlap, linfo=print):  # noqa

    if len(assemblies) > 0:
        linfo('Analyzing BLAST hits for assemblies [' + ss + ']')
        if seqtk is None:
            linfo('seqtk is not available. ' +
                  'Cannot continue. Exiting.')
            exit(0)

    for a in assemblies:

        if ('blast_hits_aa__' + ss) not in a:
            continue

        assmbl_name = a['name']
        tax_id = a['tax_id']

        parsed_hits = a['blast_hits_aa__' + ss]

        a_path = a['path']
        gc_tt = a['gc_tt']

        transcripts_nt_fasta_file = opj(
            dir_prj_transcripts, assmbl_name + '_transcripts_nt__' + ss + '.fasta')

        transcripts_nt_orf_fasta_file = opj(
            dir_prj_transcripts, assmbl_name + '_transcripts_nt_orf__' + ss + '.fasta')

        transcripts_aa_orf_fasta_file = opj(
            dir_prj_transcripts, assmbl_name + '_transcripts_aa_orf__' + ss + '.fasta')

        transcripts_nt = {}
        transcripts_nt_orf = {}
        transcripts_aa_orf = {}

        a['annotations__' + ss] = {}

        collated = collate_blast_results(parsed_hits)

        ######################################################################
        # Use seqtk to sample the assembly FASTA file for sequences with
        # BLAST hits. This increases the speed substantially when the assembly
        # file is large.
        temp_a_file = opj(dir_temp, 'temp__' + ss + '.fasta')
        temp_s_file = opj(dir_temp, 'temp__' + ss + '.txt')
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
                                  '_ann_kakapo__' + ss + '.json')

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

            a['annotations__' + ss][target_name] = {}

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

                    a['annotations__' + ss][target_name]['orf_begin'] = ann_orf_b
                    a['annotations__' + ss][target_name]['orf_end'] = ann_orf_e
                    a['annotations__' + ss][target_name]['orf_grade'] = good_orf[3]

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
                a['annotations__' + ss][target_name]['orfs_bad'] = list()
                orfs_bad_list = a['annotations__' + ss][target_name]['orfs_bad']

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

            a['annotations__' + ss][target_name]['query_name'] = query_name
            a['annotations__' + ss][target_name]['evalue'] = hit_evalue
            a['annotations__' + ss][target_name]['frame'] = abs(hit_frame)
            a['annotations__' + ss][target_name]['blast_hit_begin'] = ann_hit_b
            a['annotations__' + ss][target_name]['blast_hit_end'] = ann_hit_e

            ##################################################################
            # Collect ORF and BLAST hit annotations for downstream use.
            kakapo_json = [{}]
            kakapo_json[0]['kakapo_annotations__' + ss] = (
                a['annotations__' + ss][target_name])
            all_kakapo_results[target_name] = kakapo_json
            ##################################################################

        # --------------------------------------------------------------------

        linfo('Transcripts: ' + str(len(transcripts_nt)))

        if len(transcripts_nt) > 0:
            write_fasta(transcripts_nt, transcripts_nt_fasta_file)
            a['transcripts_nt_fasta_file__' + ss] = transcripts_nt_fasta_file
        else:
            a['transcripts_nt_fasta_file__' + ss] = None

        linfo('Transcripts with acceptable ORFs: ' +
              str(len(transcripts_nt_orf)))

        if len(transcripts_nt_orf) > 0:
            write_fasta(transcripts_nt_orf, transcripts_nt_orf_fasta_file)
            a['transcripts_nt_orf_fasta_file__' + ss] = transcripts_nt_orf_fasta_file
        else:
            a['transcripts_nt_orf_fasta_file__' + ss] = None

        if len(transcripts_aa_orf) > 0:
            write_fasta(transcripts_aa_orf, transcripts_aa_orf_fasta_file)
            a['transcripts_aa_orf_fasta_file__' + ss] = transcripts_aa_orf_fasta_file
        else:
            a['transcripts_aa_orf_fasta_file__' + ss] = None

        # --------------------------------------------------------------------
        # Save ORF and BLAST hit annotations for downstream use.
        with open(json_dump_file_path, 'w') as f:
            json.dump(all_kakapo_results, f, sort_keys=True, indent=4)
        # --------------------------------------------------------------------


def run_inter_pro_scan(ss, assemblies, email, dir_prj_ips, dir_cache_prj,
                       linfo=print):  # noqa
    delay = 0.25

    if len(assemblies) > 0:
        linfo('Running InterProScan on translated transcripts [' + ss + ']')

    for a in assemblies:

        if 'transcripts_aa_orf_fasta_file__' + ss not in a:
            continue

        aa_file = a['transcripts_aa_orf_fasta_file__' + ss]

        if aa_file is None:
            continue

        assmbl_name = a['name']

        json_dump_file_path = opj(dir_prj_ips, assmbl_name + '_ann_ips__' + ss + '.json')

        if ope(json_dump_file_path):
            linfo('InterProScan results for assembly ' + assmbl_name + ', ' +
                  ss + ' have already been downloaded [' + ss + ']')
            continue

        seqs = read_fasta(aa_file)
        seqs = OrderedDict(sorted(seqs.items(),
                                  key=lambda x: x[0].split(' ')[1],
                                  reverse=True))

        _ = opj(dir_cache_prj, assmbl_name + '_ips_jobs__' + ss)

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
              assmbl_name + ' [' + ss + ']')
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
            ips_json[0]['kakapo_annotations__' + ss] = a['annotations__' + ss][job_no_def]

            all_ips_results[job_no_def] = ips_json

        print()

        with open(json_dump_file_path, 'w') as f:
            json.dump(all_ips_results, f, sort_keys=True, indent=4)

        # Removes cached jobs file.
        osremove(_)


def gff_from_json(ss, assemblies, dir_prj_ips, dir_prj_transcripts_combined,
                  prj_name, linfo=print):  # noqa
    if len(assemblies) > 0:
        linfo('Producing GFF3 files [' + ss + ']')

    all_fas_paths = []
    all_gff_paths = []

    combined_fas_path = opj(dir_prj_transcripts_combined, prj_name + '__' + ss + '.fasta')
    combined_gff_path = opj(dir_prj_transcripts_combined, prj_name + '__' + ss + '.gff')

    for a in assemblies:

        if 'transcripts_nt_fasta_file__' + ss not in a:
            continue

        assmbl_name = a['name']
        transcripts_nt_path = a['transcripts_nt_fasta_file__' + ss]

        kakapo_json_path = opj(dir_prj_ips, assmbl_name + '_ann_kakapo__' + ss + '.json')
        ips_json_path = opj(dir_prj_ips, assmbl_name + '_ann_ips__' + ss + '.json')
        json_path = opj(dir_prj_ips, assmbl_name + '_ann__' + ss + '.json')

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
            gff_from_kakapo_ips5_json_file(ss, json_path, gff_path)

            osremove(json_path)

            all_gff_paths.append(gff_path)
            all_fas_paths.append(transcripts_nt_path)

    combine_text_files(all_fas_paths, combined_fas_path)
    combine_text_files(all_gff_paths, combined_gff_path)


def dnld_cds_for_ncbi_prot_acc(ss, prot_acc_user, prot_cds_ncbi_file, tax,
                               dir_cache_prj, linfo=print):  # noqa

    # TODO: The function downloads more CDS than are strictly required,
    #       then filters out the unneeded ones. Sometimes this causes
    #       large amounts of data to be downloaded unnecessarily.

    pickle_file = opj(dir_cache_prj, 'ncbi_prot_cds_cache__' + ss)
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
              'NCBI protein accessions [' + ss + ']')
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

                prot_name = re.findall(r'\[protein=(.*?)\]', rec)
                if len(prot_name) > 0:
                    prot_name = prot_name[0]
                else:
                    prot_name = re.findall(r'\[gene=(.*?)\]', rec)
                    if len(prot_name) > 0:
                        prot_name = prot_name[0]
                    else:
                        prot_name = ''

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
