"""Kakapo workflow: Prepare."""

from os.path import join as opj
from kakapo.utils.misc import make_dirs


def prepare_output_directories(dir_out, prj_name):
    # ToDo: Lock cache files in case of parallel execution -------------------
    dir_temp = opj(dir_out, '00-temp')
    make_dirs(dir_temp)

    dir_cache = opj(dir_out, '00-cache')
    make_dirs(dir_cache)

    dir_cache_pfam_acc = opj(dir_cache, 'pfam-uniprot-accessions')
    make_dirs(dir_cache_pfam_acc)

    dir_cache_fq_minlen = opj(dir_cache, 'min-acceptable-read-lengths')
    make_dirs(dir_cache_fq_minlen)

    dir_cache_prj = opj(dir_cache, 'projects', prj_name)
    make_dirs(dir_cache_prj)

    dir_cache_refseqs = opj(dir_cache, 'ref-seqs')
    make_dirs(dir_cache_refseqs)

    dir_prj = opj(dir_out, '02-project-specific', prj_name)
    make_dirs(dir_prj)

    dir_prj_logs = opj(dir_prj, '00-logs')
    make_dirs(dir_prj_logs)

    dir_prj_queries = opj(dir_prj, '01-queries')
    make_dirs(dir_prj_queries)

    dir_prj_blast_results_fa_trim = opj(dir_prj,
                                        '02-filtered-fa-blast-results')
    make_dirs(dir_prj_blast_results_fa_trim)

    dir_prj_vsearch_results_fa_trim = opj(dir_prj,
                                          '03-filtered-fq-vsearch-results')
    make_dirs(dir_prj_vsearch_results_fa_trim)

    dir_prj_spades_assemblies = opj(dir_prj, '04-spades-assemblies')
    make_dirs(dir_prj_spades_assemblies)

    dir_prj_blast_assmbl = opj(dir_prj, '05-assemblies-blast-db-data')
    make_dirs(dir_prj_blast_assmbl)

    dir_prj_assmbl_blast_results = opj(dir_prj, '06-assemblies-blast-results')
    make_dirs(dir_prj_assmbl_blast_results)

    dir_prj_transcripts = opj(dir_prj, '07-transcripts')
    make_dirs(dir_prj_transcripts)

    dir_prj_ips = dir_prj_transcripts

    dir_prj_transcripts_combined = opj(dir_prj, '08-transcripts-combined')
    make_dirs(dir_prj_transcripts_combined)

    dir_global = opj(dir_out, '01-global')
    make_dirs(dir_global)

    dir_fq_data = opj(dir_global, '01-sra-fq-data')
    make_dirs(dir_fq_data)

    dir_fq_trim_data = opj(dir_global, '02-trimmed-fq-data')
    make_dirs(dir_fq_trim_data)

    dir_fq_cor_data = opj(dir_global, '03-corrected-fq-data')
    make_dirs(dir_fq_cor_data)

    dir_fq_filter_bt2_data = opj(dir_global, '04-bowtie2-filtered-fq-data')
    make_dirs(dir_fq_filter_bt2_data)

    dir_fq_filter_krkn2_data = opj(dir_global, '05-kraken2-filtered-fq-data')
    make_dirs(dir_fq_filter_krkn2_data)

    dir_fa_trim_data = opj(dir_global, '06-fa-data')
    make_dirs(dir_fa_trim_data)

    dir_blast_fa_trim = opj(dir_global, '07-fa-blast-db-data')
    make_dirs(dir_blast_fa_trim)

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
                'dir_fq_filter_bt2_data': dir_fq_filter_bt2_data,
                'dir_fq_filter_krkn2_data': dir_fq_filter_krkn2_data,
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
