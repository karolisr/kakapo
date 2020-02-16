# -*- coding: utf-8 -*-

"""Kakapo workflow: Search Reads."""

from os import remove as osremove
from os.path import exists as ope
from os.path import join as opj
from shutil import copyfile

from kakapo.blast import BLST_RES_COLS_1
from kakapo.blast import run_blast
from kakapo.helpers import combine_text_files
from kakapo.helpers import keep_unique_lines_in_file
from kakapo.helpers import make_dir
from kakapo.seqtk import seqtk_fq_to_fa, seqtk_extract_reads
from kakapo.vsearch import run_cluster_fast, run_vsearch


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
