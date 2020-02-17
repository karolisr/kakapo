# -*- coding: utf-8 -*-

"""Kakapo workflow: Search Assemblies."""

import pickle

from os import remove as osremove
from os.path import exists as ope
from os.path import join as opj

from kakapo.bioio import read_fasta
from kakapo.blast import BLST_RES_COLS_2
from kakapo.blast import parse_blast_results_file
from kakapo.blast import run_blast
from kakapo.config import PICKLE_PROTOCOL


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
