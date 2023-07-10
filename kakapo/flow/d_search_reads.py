"""Kakapo workflow: Search Reads."""

import pickle

from os import remove as osremove
from os.path import basename
from os.path import exists as ope
from os.path import join as opj
from shutil import copyfile

from kakapo.tools.bioio import read_fasta
from kakapo.tools.bioio import seq_records_to_dict
from kakapo.tools.blast import BLST_RES_COLS_1
from kakapo.tools.blast import run_blast
from kakapo.tools.config import PICKLE_PROTOCOL
from kakapo.tools.seq import SEQ_TYPE_AA
from kakapo.tools.seqtk import seqtk_fq_to_fa, seqtk_extract_reads
from kakapo.tools.vsearch import run_cluster_fast, run_vsearch
from kakapo.utils.logging import Log
from kakapo.utils.misc import combine_text_files
from kakapo.utils.misc import keep_unique_lines_in_file
from kakapo.utils.misc import make_dirs


MT = 'mitochondrion'
PT = 'plastid'
ORGANELLES = (MT, PT)


def run_tblastn_on_reads(se_fastq_files, pe_fastq_files, aa_queries_file,
                         tblastn, blast_1_evalue, blast_1_max_hsps,
                         blast_1_qcov_hsp_perc, blast_1_best_hit_overhang,
                         blast_1_best_hit_score_edge, blast_1_max_target_seqs,
                         dir_blast_results_fa_trim, fpatt, ss, threads,
                         seqtk, vsearch, dir_cache_prj, ss_organelle):

    changed_blast_1 = False

    if len(se_fastq_files) > 0 or len(pe_fastq_files) > 0:
        print()
        Log.inf('Running BLAST on reads:', ss)
        if tblastn is None:
            Log.err('tblastn is not available. Cannot continue. Exiting.')
            exit(0)

        if vsearch is None:
            Log.err('vsearch is not available. Cannot continue. Exiting.')
            exit(0)

        if seqtk is None:
            Log.err('seqtk is not available. Cannot continue. Exiting.')
            exit(0)

    cache_file = opj(dir_cache_prj, 'blast_1_settings_cache__' + ss)

    pickled = dict()
    settings = {'blast_1_evalue': blast_1_evalue,
                'blast_1_max_hsps': blast_1_max_hsps,
                'blast_1_qcov_hsp_perc': blast_1_qcov_hsp_perc,
                'blast_1_best_hit_overhang': blast_1_best_hit_overhang,
                'blast_1_best_hit_score_edge': blast_1_best_hit_score_edge,
                'blast_1_max_target_seqs': blast_1_max_target_seqs,
                'queries': seq_records_to_dict(
                    read_fasta(aa_queries_file, SEQ_TYPE_AA))}

    Log.msg('evalue:', str(blast_1_evalue))
    Log.msg('max_hsps:', str(blast_1_max_hsps))
    Log.msg('qcov_hsp_perc:', str(blast_1_qcov_hsp_perc))
    Log.msg('best_hit_overhang:', str(blast_1_best_hit_overhang))
    Log.msg('best_hit_score_edge:', str(blast_1_best_hit_score_edge))
    Log.msg('max_target_seqs:', str(blast_1_max_target_seqs))
    print()

    # FixMe: Expose in configuration files?
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
        fq_organelle = se_fastq_files[se]['organelle']

        if ope(out_f_fasta) and ope(cache_file):
            with open(cache_file, 'rb') as f:
                pickled = pickle.load(f)

        if ope(out_f_fasta) and pickled == settings:
            # Log.msg('The provided BLAST settings and query sequences did '
            #         'not change since the previous run.')
            Log.msg('BLAST results already exist:', se)

        else:
            if fq_organelle is not None:
                if fq_organelle != ss_organelle:
                    se_fastq_files[se]['blast_results_path' + '__' + ss] = None
                    continue
            else:
                if ss_organelle in ORGANELLES:
                    se_fastq_files[se]['blast_results_path' + '__' + ss] = None
                    continue
            changed_blast_1 = True
            make_dirs(dir_results)
            Log.msg('Running tblastn on:', basename(blast_db_path))
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

            # Log.inf('Extracting unique BLAST hits using Seqtk.')

            keep_unique_lines_in_file(out_f)

            seqtk_extract_reads(seqtk, fq_path, out_f_fastq, out_f)
            seqtk_fq_to_fa(seqtk, out_f_fastq, out_f_fasta, gz_out=False)

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
        fq_organelle = pe_fastq_files[pe]['organelle']

        if ope(out_f_fasta) and ope(cache_file):
            with open(cache_file, 'rb') as f:
                pickled = pickle.load(f)

        if ope(out_f_fasta) and pickled == settings:
            # Log.msg('The provided BLAST settings and query sequences did '
            #         'not change since the previous run.')
            Log.msg('BLAST results already exist:', pe)

        else:
            if fq_organelle is not None:
                if fq_organelle != ss_organelle:
                    pe_fastq_files[pe]['blast_results_path' + '__' + ss] = None
                    continue
            else:
                if ss_organelle in ORGANELLES:
                    pe_fastq_files[pe]['blast_results_path' + '__' + ss] = None
                    continue
            changed_blast_1 = True
            make_dirs(dir_results)
            pe_trim_files = zip(blast_db_paths, out_fs, fq_paths, out_fs_fastq,
                                out_fs_fasta)
            for x in pe_trim_files:
                Log.msg('Running tblastn on:', basename(x[0]))
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

                # Log.msg('Extracting unique BLAST hits using Seqtk.')

                keep_unique_lines_in_file(x[1])

                seqtk_extract_reads(seqtk, x[2], x[3], x[1])
                seqtk_fq_to_fa(seqtk, x[3], x[4], gz_out=False)

                osremove(x[1])
                osremove(x[3])

            combine_text_files(out_fs_fasta, out_f_fasta)

            out_f_fasta_temp = out_f_fasta + '_temp'
            copyfile(out_f_fasta, out_f_fasta_temp)
            run_cluster_fast(vsearch, ident, out_f_fasta_temp, out_f_fasta)
            osremove(out_f_fasta_temp)

            for x in out_fs_fasta:
                osremove(x)

    with open(cache_file, 'wb') as f:
        pickle.dump(settings, f, protocol=PICKLE_PROTOCOL)

    return changed_blast_1


def run_vsearch_on_reads(se_fastq_files, pe_fastq_files, vsearch,
                         dir_vsearch_results_fa_trim, fpatt, ss, seqtk):

    if len(se_fastq_files) > 0 or len(pe_fastq_files) > 0:
        if vsearch is None:
            Log.err('vsearch is not available. Cannot continue. Exiting.')
            exit(0)
        if seqtk is None:
            Log.err('seqtk is not available. Cannot continue. Exiting.')
            exit(0)

    # FixMe: Expose in configuration files?
    ident = 0.85

    for se in se_fastq_files:
        dir_results = opj(dir_vsearch_results_fa_trim, se)
        min_acc_len = se_fastq_files[se]['min_acc_len']
        blast_results_fa_path = se_fastq_files[se]['blast_results_path' + '__' + ss]
        fq_path = se_fastq_files[se]['filter_path_fq']
        out_f = opj(dir_results, se + '__' + ss + '.txt')
        out_f_fastq = out_f.replace('.txt', '.fastq')
        se_fastq_files[se]['vsearch_results_path' + '__' + ss] = out_f_fastq

        if blast_results_fa_path is None:
            se_fastq_files[se]['vsearch_results_path' + '__' + ss] = None
            continue

        if ope(out_f_fastq):
            Log.msg('Vsearch results already exist:', se)
        else:
            make_dirs(dir_results)
            Log.msg('Running vsearch on:', basename(fq_path))
            run_vsearch(vsearch,
                        ident=ident,
                        q_file=blast_results_fa_path,
                        db_file=fq_path,
                        out_file=out_f,
                        minlen=min_acc_len)

            # Log.msg('Extracting unique vsearch hits using Seqtk.')
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

        if blast_results_fa_path is None:
            pe_fastq_files[pe]['vsearch_results_path' + '__' + ss] = None
            continue

        if ope(out_fs_fastq[0]) and ope(out_fs_fastq[1]) and \
           ope(out_fs_fastq[2]) and ope(out_fs_fastq[3]):
            Log.msg('Vsearch results already exist:', pe)
        else:
            make_dirs(dir_results)
            pe_trim_files = zip(fq_paths, out_fs, out_fs_fastq)
            for x in pe_trim_files:
                Log.msg('Running vsearch on:', basename(x[0]))
                run_vsearch(vsearch,
                            ident=ident,
                            q_file=blast_results_fa_path,
                            db_file=x[0],
                            out_file=x[1],
                            minlen=min_acc_len)

            # Log.msg('Extracting unique vsearch hits from paired files '
            #         'using Seqtk.')

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

            # Log.msg('Extracting unique vsearch hits from unpaired files '
            #         'using Seqtk.')

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
