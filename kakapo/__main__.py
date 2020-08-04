#!/usr/bin/env python3

"""Kakapo main file."""

import inspect
import os
import sys

##############################################################################
SCRIPT_FILE_PATH = inspect.getfile(inspect.currentframe())
SCRIPT_DIR_PATH = os.path.dirname(os.path.abspath(SCRIPT_FILE_PATH))
KAKAPO_DIR_PATH = os.path.sep.join(SCRIPT_DIR_PATH.split(os.path.sep)[0:-1])
sys.path.insert(0, KAKAPO_DIR_PATH)
##############################################################################

import argparse
import logging

from collections import OrderedDict
from io import StringIO
from os import stat
from os.path import exists as ope
from os.path import join as opj
from shutil import rmtree
from sys import exit

from joblib import Parallel, delayed

from ncbi_taxonomy_local import taxonomy

from kakapo import __version__, __script_name__

from kakapo.utils import dependencies as deps
from kakapo.utils.logging import prepare_logger
from kakapo.utils.misc import make_dirs
from kakapo.utils.misc import time_stamp

from kakapo.tools.config import CONYELL, CONSRED, CONGREE, CONBLUE, CONSDFL
from kakapo.tools.config import DIR_CFG, DIR_DEP, DIR_TAX, DIR_KRK
from kakapo.tools.config import SCRIPT_INFO
from kakapo.tools.config import OS_ID, OS_STR, DIST_ID, DEBIAN_DISTS, REDHAT_DISTS, SUPPORTED_DISTS
from kakapo.tools.config import THREADS, RAM
from kakapo.tools.config_file_parse import config_file_parse
from kakapo.tools.config_file_parse import ss_file_parse

from kakapo.tools.seq import SEQ_TYPE_AA

from kakapo.tools.bioio import read_fasta
from kakapo.tools.bioio import seq_records_to_dict
from kakapo.tools.transl_tables import TranslationTable

from kakapo.flow.a_prepare import prepare_output_directories

from kakapo.flow.b_process_queries import combine_aa_fasta
from kakapo.flow.b_process_queries import dnld_pfam_uniprot_seqs
from kakapo.flow.b_process_queries import dnld_prot_seqs
from kakapo.flow.b_process_queries import filter_queries
from kakapo.flow.b_process_queries import pfam_uniprot_accessions
from kakapo.flow.b_process_queries import user_aa_fasta
from kakapo.flow.b_process_queries import user_entrez_search
from kakapo.flow.b_process_queries import user_protein_accessions

from kakapo.flow.c_process_reads import dnld_sra_fastq_files
from kakapo.flow.c_process_reads import dnld_sra_info
from kakapo.flow.c_process_reads import file_name_patterns
from kakapo.flow.c_process_reads import filtered_fq_to_fa
from kakapo.flow.c_process_reads import makeblastdb_fq
from kakapo.flow.c_process_reads import min_accept_read_len
from kakapo.flow.c_process_reads import run_bt2_fq
from kakapo.flow.c_process_reads import run_kraken2
from kakapo.flow.c_process_reads import run_rcorrector
from kakapo.flow.c_process_reads import run_trimmomatic
from kakapo.flow.c_process_reads import user_fastq_files

from kakapo.flow.d_search_reads import run_tblastn_on_reads
from kakapo.flow.d_search_reads import run_vsearch_on_reads

from kakapo.flow.e_process_assmbl import combine_assemblies
from kakapo.flow.e_process_assmbl import makeblastdb_assemblies
from kakapo.flow.e_process_assmbl import run_spades

from kakapo.flow.f_search_assmbl import run_tblastn_on_assemblies
from kakapo.flow.g_find_orfs import find_orfs_translate
from kakapo.flow.h_prepare_gff import gff_from_json
from kakapo.flow.i_inter_pro_scan import run_inter_pro_scan
from kakapo.flow.j_dnld_aa_query_cds import dnld_cds_for_ncbi_prot_acc

# Silence urllib3 overly verbose logging.
logging.getLogger("urllib3").setLevel(logging.WARNING)

# Command line arguments -----------------------------------------------------
USAGE = '{} --cfg path/to/config_file ' \
        '--ss path/to/search_strategies_file'.format(__script_name__)

PARSER = argparse.ArgumentParser(
    prog=__script_name__,
    formatter_class=argparse.RawTextHelpFormatter,
    usage=USAGE,
    description=None,
    epilog=None,
    prefix_chars='-',
    add_help=False)

PARSER.add_argument(
    '--cfg',
    type=str,
    required=False,
    metavar='path',
    dest='CONFIG_FILE_PATH',
    help='Path to a {} project configuration file.'.format(__script_name__))

PARSER.add_argument(
    '--ss',
    type=str,
    required=False,
    metavar='path',
    dest='SS_FILE_PATH',
    help='Path to a {} search strategies file.'.format(__script_name__))

PARSER.add_argument(
    '--stop-after-filter',
    action='store_true',
    required=False,
    dest='STOP_AFTER_FILTER',
    help='Stop {} after Kraken2/Bowtie2 filtering '
         'step.'.format(__script_name__))

PARSER.add_argument(
    '--force-deps',
    action='store_true',
    default=False,
    required=False,
    dest='FORCE_DEPS',
    help='Force the use of {}-installed dependencies,\neven if they are '
         'already available on the system.'.format(__script_name__))

PARSER.add_argument(
    '--install-deps',
    action='store_true',
    default=False,
    required=False,
    dest='INSTALL_DEPS',
    help='Install {} dependencies and quit.'.format(__script_name__))

PARSER.add_argument(
    '--dnld-kraken-dbs',
    action='store_true',
    default=False,
    required=False,
    dest='DNLD_KRAKEN_DBS',
    help='Download Kraken2 databases and quit.')

PARSER.add_argument(
    '--clean-config-dir',
    action='store_true',
    required=False,
    dest='CLEAN_CONFIG_DIR',
    help='Remove all downloaded software dependencies and\n'
         'NCBI taxonomy data.')

PARSER.add_argument(
    '-v', '--version',
    action='store_true',
    required=False,
    dest='PRINT_VERSION',
    help='Print {} version.'.format(__script_name__))

PARSER.add_argument(
    '-h', '--help',
    action='store_true',
    required=False,
    dest='PRINT_HELP',
    help='Print {} help information.'.format(__script_name__))

ARGS = PARSER.parse_args()

CLEAN_CONFIG_DIR = ARGS.CLEAN_CONFIG_DIR
CONFIG_FILE_PATH = ARGS.CONFIG_FILE_PATH
SS_FILE_PATH = ARGS.SS_FILE_PATH
STOP_AFTER_FILTER = ARGS.STOP_AFTER_FILTER
FORCE_DEPS = ARGS.FORCE_DEPS
INSTALL_DEPS = ARGS.INSTALL_DEPS
DNLD_KRAKEN_DBS = ARGS.DNLD_KRAKEN_DBS
PRINT_VERSION = ARGS.PRINT_VERSION
PRINT_HELP = ARGS.PRINT_HELP

if PRINT_HELP is True:
    print(SCRIPT_INFO)
    PARSER.print_help()
    exit(0)

if PRINT_VERSION is True:
    print(__script_name__ + ' v' + __version__)
    exit(0)

if CLEAN_CONFIG_DIR is True and ope(DIR_CFG):
    print('Removing configuration directory: ' + DIR_CFG)
    rmtree(DIR_CFG)
    exit(0)
elif CLEAN_CONFIG_DIR is True:
    print('Configuration directory does not exist. Nothing to do.')
    exit(0)

if CLEAN_CONFIG_DIR is False and CONFIG_FILE_PATH is not None:
    if not ope(CONFIG_FILE_PATH):
        print('Configuration file ' + CONFIG_FILE_PATH + ' does not exist.')
        exit(0)
elif INSTALL_DEPS is True or DNLD_KRAKEN_DBS is True:
    pass
else:
    print(SCRIPT_INFO)
    print(CONYELL +
          'Configuration file was not provided. Nothing to do.' +
          CONSDFL)
    print()
    print('-' * 80)
    PARSER.print_help()
    print('-' * 80)
    print()
    exit(0)

if CLEAN_CONFIG_DIR is False and SS_FILE_PATH is not None:
    if not ope(SS_FILE_PATH):
        print('Search strategies file ' + SS_FILE_PATH + ' does not exist.')
        exit(0)
elif INSTALL_DEPS is True or DNLD_KRAKEN_DBS is True or SS_FILE_PATH is None:
    pass
else:
    pass

print(SCRIPT_INFO)

# ----------------------------------------------------------------------------


def main():
    """Run the script."""
    # Prepare initial logger (before we know the log file path) --------------
    prj_log_file_suffix = time_stamp() + '.log'
    log_stream = StringIO()
    log, _ = prepare_logger(console=True, stream=log_stream)
    linfo = log.info

    # Prepare configuration directory ----------------------------------------
    if ope(DIR_CFG):
        linfo(CONGREE + 'Found configuration directory: ' + DIR_CFG)
    else:
        linfo(CONBLUE + 'Creating configuration directory: ' + DIR_CFG)
        make_dirs(DIR_CFG)

    # Check for dependencies -------------------------------------------------
    linfo(CONBLUE + 'Checking for dependencies')
    make_dirs(DIR_DEP)
    make_dirs(DIR_KRK)
    seqtk = deps.dep_check_seqtk(DIR_DEP, FORCE_DEPS)
    trimmomatic, adapters = deps.dep_check_trimmomatic(DIR_DEP)
    fasterq_dump = deps.dep_check_sra_toolkit(DIR_DEP, OS_ID, DIST_ID,
                                              DEBIAN_DISTS, REDHAT_DISTS,
                                              FORCE_DEPS)
    makeblastdb, _, tblastn = deps.dep_check_blast(DIR_DEP, OS_ID, DIST_ID,
                                              DEBIAN_DISTS, REDHAT_DISTS,
                                              FORCE_DEPS)
    vsearch = deps.dep_check_vsearch(DIR_DEP, FORCE_DEPS)
    spades = deps.dep_check_spades(DIR_DEP, OS_ID, FORCE_DEPS)
    bowtie2, bowtie2_build = deps.dep_check_bowtie2(DIR_DEP, OS_ID, FORCE_DEPS)
    rcorrector = deps.dep_check_rcorrector(DIR_DEP, FORCE_DEPS)
    kraken2, kraken2_build = deps.dep_check_kraken2(DIR_DEP, FORCE_DEPS)

    linfo(CONBLUE + 'Checking for available Kraken2 databases')
    kraken2_dbs = deps.download_kraken2_dbs(DIR_KRK)
    for db in sorted(kraken2_dbs.keys()):
        linfo('Found Kraken2 database: ' + db)

    if INSTALL_DEPS is True or DNLD_KRAKEN_DBS is True:
        exit(0)

    # Initialize NCBI taxonomy database --------------------------------------
    tax = taxonomy(data_dir_path=DIR_TAX, linfo=linfo)

    # Parse configuration file -----------------------------------------------
    linfo(CONBLUE + 'Reading configuration file: ' + CONFIG_FILE_PATH)
    _ = config_file_parse(CONFIG_FILE_PATH, tax, linfo)

    allow_no_stop_cod = _['allow_no_stop_cod']
    allow_no_strt_cod = _['allow_no_strt_cod']
    allow_non_aug = _['allow_non_aug']

    blast_1_evalue = _['blast_1_evalue']
    blast_1_max_hsps = _['blast_1_max_hsps']
    blast_1_qcov_hsp_perc = _['blast_1_qcov_hsp_perc']
    blast_1_best_hit_overhang = _['blast_1_best_hit_overhang']
    blast_1_best_hit_score_edge = _['blast_1_best_hit_score_edge']
    blast_1_max_target_seqs = _['blast_1_max_target_seqs']

    blast_2_evalue = _['blast_2_evalue']
    blast_2_max_hsps = _['blast_2_max_hsps']
    blast_2_qcov_hsp_perc = _['blast_2_qcov_hsp_perc']
    blast_2_best_hit_overhang = _['blast_2_best_hit_overhang']
    blast_2_best_hit_score_edge = _['blast_2_best_hit_score_edge']
    blast_2_max_target_seqs = _['blast_2_max_target_seqs']

    dir_out = _['output_directory']
    email = _['email']
    ncbi_longevity = _['ncbi_longevity']
    fq_pe = _['fq_pe']
    fq_se = _['fq_se']
    should_run_rcorrector = _['should_run_rcorrector']
    should_run_ipr = _['should_run_ipr']
    bt2_order = _['bt2_order']
    kraken_confidence = _['kraken_confidence']
    krkn_order = _['krkn_order']
    prepend_assmbl = _['prepend_assmbl']
    prj_name = _['project_name']
    sras = _['sras']
    tax_group = _['tax_group']
    # tax_group_name = _['tax_group_name']
    tax_ids_user = _['tax_ids']
    user_assemblies = _['assmbl']

    # Parse search strategies file -------------------------------------------
    if SS_FILE_PATH is not None:
        linfo(CONBLUE + 'Reading search strategies file: ' + SS_FILE_PATH)
        sss = ss_file_parse(SS_FILE_PATH, linfo)
    else:
        linfo(CONSRED + 'Search strategies file was not provided. ' +
              'Will process reads, assemblies and then stop.' + CONSDFL)
        sss = dict()

    # Create output directory ------------------------------------------------
    if dir_out is not None:
        if ope(dir_out):
            linfo(CONGREE + 'Found output directory: ' + dir_out)
        else:
            linfo(CONBLUE + 'Creating output directory: ' + dir_out)
            make_dirs(dir_out)

    # Write Kakapo version information to the output directory ---------------
    version_file = opj(dir_out, 'kakapo_version.txt')
    if ope(version_file):
        with open(version_file, 'r') as f:
            version_prev = f.read().strip()
            if __version__ != version_prev:
                linfo(CONSRED + 'The output directory contains data produced by a ' +
                      'different version of Kakapo: ' + version_prev +
                      '. The currently running version is: ' + __version__ + '. ' +
                      'Delete "kakapo_version.txt" file located in the ' +
                      'output directory if you would like to continue.' +
                      CONSDFL)
                exit(0)
    with open(version_file, 'w') as f:
        f.write(__version__)

    # Create subdirectories in the output directory --------------------------
    _ = prepare_output_directories(dir_out, prj_name)

    dir_temp = _['dir_temp']
    dir_cache_pfam_acc = _['dir_cache_pfam_acc']
    dir_cache_fq_minlen = _['dir_cache_fq_minlen']
    dir_cache_prj = _['dir_cache_prj']
    dir_cache_refseqs = _['dir_cache_refseqs']
    dir_prj_logs = _['dir_prj_logs']
    dir_prj_queries = _['dir_prj_queries']
    dir_fq_data = _['dir_fq_data']
    dir_fq_cor_data = _['dir_fq_cor_data']
    dir_fq_trim_data = _['dir_fq_trim_data']
    dir_fq_filter_bt2_data = _['dir_fq_filter_bt2_data']
    dir_fq_filter_krkn2_data = _['dir_fq_filter_krkn2_data']
    dir_fa_trim_data = _['dir_fa_trim_data']
    dir_blast_fa_trim = _['dir_blast_fa_trim']
    dir_prj_blast_results_fa_trim = _['dir_prj_blast_results_fa_trim']
    dir_prj_vsearch_results_fa_trim = _['dir_prj_vsearch_results_fa_trim']
    dir_prj_spades_assemblies = _['dir_prj_spades_assemblies']
    dir_prj_blast_assmbl = _['dir_prj_blast_assmbl']
    dir_prj_assmbl_blast_results = _['dir_prj_assmbl_blast_results']
    dir_prj_transcripts = _['dir_prj_transcripts']
    dir_prj_ips = _['dir_prj_ips']
    dir_prj_transcripts_combined = _['dir_prj_transcripts_combined']

    # Prepare logger ---------------------------------------------------------
    prj_log_file = opj(dir_prj_logs, prj_name + '_' + prj_log_file_suffix)
    with open(prj_log_file, 'w') as f:
        f.write(SCRIPT_INFO.strip() + '\n\n' + log_stream.getvalue())
    log, _ = prepare_logger(console=True, stream=None, file=prj_log_file)
    linfo = log.info

    # Resolve descending taxonomy nodes --------------------------------------
    tax_ids = tax.all_descending_taxids_for_taxids([tax_group])

    # Pfam uniprot accessions ------------------------------------------------
    pfam_uniprot_acc = OrderedDict()
    for ss in sss:
        pfam_acc = sss[ss]['pfam_families']
        pfam_uniprot_acc[ss] = pfam_uniprot_accessions(ss, pfam_acc, tax_ids,
                                                       dir_cache_pfam_acc,
                                                       linfo)

    # Download Pfam uniprot sequences if needed ------------------------------
    aa_uniprot_files = OrderedDict()
    for ss in sss:
        aa_uniprot_files[ss] = opj(dir_prj_queries, 'aa_uniprot__' + ss +
                                   '.fasta')
        dnld_pfam_uniprot_seqs(ss, pfam_uniprot_acc[ss], aa_uniprot_files[ss],
                               dir_cache_prj, linfo)

    # User provided entrez query ---------------------------------------------
    prot_acc_user_from_query = OrderedDict()
    for ss in sss:
        entrez_queries = sss[ss]['entrez_search_queries']
        prot_acc_user_from_query[ss] = user_entrez_search(ss, entrez_queries,
                                                          dir_cache_prj,
                                                          ncbi_longevity,
                                                          linfo)

    # User provided protein accessions ---------------------------------------
    prot_acc_user = OrderedDict()
    for ss in sss:
        prot_acc_all = sorted(set(sss[ss]['ncbi_accessions_aa'] +
                                  prot_acc_user_from_query[ss]))
        prot_acc_user[ss] = user_protein_accessions(ss, prot_acc_all,
                                                    dir_cache_prj, tax, linfo)

    # Download from NCBI if needed -------------------------------------------
    aa_prot_ncbi_files = OrderedDict()
    for ss in sss:
        aa_prot_ncbi_files[ss] = opj(dir_prj_queries, 'aa_prot_ncbi__' + ss +
                                     '.fasta')
        prot_acc_user[ss] = dnld_prot_seqs(ss, prot_acc_user[ss],
                                           aa_prot_ncbi_files[ss], linfo)

    # User provided protein sequences ----------------------------------------
    aa_prot_user_files = OrderedDict()
    for ss in sss:
        user_queries = sss[ss]['fasta_files_aa']
        aa_prot_user_files[ss] = opj(dir_prj_queries, 'aa_prot_user__' + ss +
                                     '.fasta')
        user_aa_fasta(ss, user_queries, aa_prot_user_files[ss], linfo)

    # Combine all AA queries -------------------------------------------------
    aa_queries_files = OrderedDict()
    for ss in sss:
        aa_queries_files[ss] = opj(dir_prj_queries, 'aa_all__' + ss + '.fasta')
        combine_aa_fasta(ss, [aa_uniprot_files[ss], aa_prot_ncbi_files[ss],
                              aa_prot_user_files[ss]], aa_queries_files[ss],
                         linfo)

    # Filter AA queries ------------------------------------------------------
    prot_acc_user_filtered = OrderedDict()
    for ss in sss:
        min_query_length = sss[ss]['min_query_length']
        max_query_length = sss[ss]['max_query_length']
        max_query_identity = sss[ss]['max_query_identity']

        # Dereplicate all queries
        filter_queries(ss, aa_queries_files[ss], min_query_length,
                       max_query_length, max_query_identity,
                       vsearch, prot_acc_user[ss], overwrite=True, linfo=linfo)

        # Dereplicate only NCBI queries. CDS for these will be downloaded
        # later for reference.
        if ope(aa_prot_ncbi_files[ss]):
            prot_acc_user_filtered[ss] = filter_queries(
                ss, aa_prot_ncbi_files[ss], min_query_length, max_query_length,
                max_query_identity, vsearch, prot_acc_user[ss],
                overwrite=False, linfo=lambda x: x)

    # Download SRA run metadata if needed ------------------------------------
    sra_runs_info, sras_acceptable = dnld_sra_info(sras, dir_cache_prj, linfo)

    # Download SRA run FASTQ files if needed ---------------------------------
    x, y, z = dnld_sra_fastq_files(sras_acceptable, sra_runs_info, dir_fq_data,
                                   fasterq_dump, THREADS, dir_temp, linfo)

    se_fastq_files_sra = x
    pe_fastq_files_sra = y
    sra_runs_info = z

    # User provided FASTQ files ----------------------------------------------
    se_fastq_files_usr, pe_fastq_files_usr = user_fastq_files(fq_se, fq_pe,
                                                              linfo)

    # Collate FASTQ file info ------------------------------------------------
    se_fastq_files = se_fastq_files_sra.copy()
    se_fastq_files.update(se_fastq_files_usr)
    pe_fastq_files = pe_fastq_files_sra.copy()
    pe_fastq_files.update(pe_fastq_files_usr)

    def gc_tt(k, d, tax):
        taxid = d[k]['tax_id']

        gc = tax.genetic_code_for_taxid(taxid)

        d[k]['gc_id'] = gc
        d[k]['gc_tt'] = TranslationTable(gc)

        gc_mito = None
        tt_mito = None

        gc_plastid = None
        tt_plastid = None

        if tax.is_eukaryote(taxid) is True:
            gc_mito = tax.mito_genetic_code_for_taxid(taxid)
            if gc_mito != '0':
                tt_mito = TranslationTable(gc_mito)

            if tax.contains_plastid(taxid) is True:
                gc_plastid = tax.plastid_genetic_code_for_taxid(taxid)
                if gc_plastid != '0':
                    tt_plastid = TranslationTable(gc_plastid)

        d[k]['gc_id_mito'] = gc_mito
        d[k]['gc_tt_mito'] = tt_mito

        d[k]['gc_id_plastid'] = gc_plastid
        d[k]['gc_tt_plastid'] = tt_plastid

    for se in se_fastq_files:
        gc_tt(se, se_fastq_files, tax)

    for pe in pe_fastq_files:
        gc_tt(pe, pe_fastq_files, tax)

    # Minimum acceptable read length -----------------------------------------
    min_accept_read_len(se_fastq_files, pe_fastq_files, dir_temp,
                        dir_cache_fq_minlen, vsearch, linfo)

    # Run Rcorrector ---------------------------------------------------------
    run_rcorrector(se_fastq_files, pe_fastq_files, dir_fq_cor_data, rcorrector,
                   THREADS, dir_temp, should_run_rcorrector, linfo)

    # File name patterns -----------------------------------------------------
    a, b, c, d, e = file_name_patterns()

    pe_trim_fq_file_patterns = a
    pe_trim_fa_file_patterns = b
    pe_blast_db_file_patterns = c
    pe_blast_results_file_patterns = d
    pe_vsearch_results_file_patterns = e

    # Run Trimmomatic --------------------------------------------------------
    run_trimmomatic(se_fastq_files, pe_fastq_files, dir_fq_trim_data,
                    trimmomatic, adapters, pe_trim_fq_file_patterns, THREADS,
                    linfo)

    # Run Bowtie 2 -----------------------------------------------------------
    run_bt2_fq(se_fastq_files, pe_fastq_files, dir_fq_filter_bt2_data,
               bowtie2, bowtie2_build, THREADS, dir_temp, bt2_order,
               pe_trim_fq_file_patterns, tax, dir_cache_refseqs,
               linfo)

    # Run Kraken2 ------------------------------------------------------------
    run_kraken2(krkn_order, kraken2_dbs, se_fastq_files, pe_fastq_files,
                dir_fq_filter_krkn2_data, kraken_confidence, kraken2, THREADS,
                dir_temp, pe_trim_fq_file_patterns, linfo)

    se_fastq_files = OrderedDict(se_fastq_files)
    pe_fastq_files = OrderedDict(pe_fastq_files)

    se_fastq_files = OrderedDict(sorted(se_fastq_files.items(),
                                        key=lambda x: x[1]['filter_path_fq']))

    pe_fastq_files = OrderedDict(sorted(pe_fastq_files.items(),
                                        key=lambda x: x[1]['filter_path_fq']))

    # Stop After Filter ------------------------------------------------------
    if STOP_AFTER_FILTER is True:
        linfo(CONSRED +
              'Stopping after Kraken2/Bowtie2 filtering step as requested.' +
              CONSDFL)
        exit(0)

    # Convert filtered FASTQ files to FASTA ----------------------------------
    filtered_fq_to_fa(se_fastq_files, pe_fastq_files, dir_fa_trim_data, seqtk,
                      pe_trim_fa_file_patterns, linfo)

    # Run makeblastdb on reads -----------------------------------------------
    makeblastdb_fq(se_fastq_files, pe_fastq_files, dir_blast_fa_trim,
                   makeblastdb, pe_blast_db_file_patterns, linfo)

    # Check if there are any query sequences.
    any_queries = False
    for ss in sss:
        if stat(aa_queries_files[ss]).st_size == 0:
            continue
        else:
            any_queries = True

    # Run tblastn on reads ---------------------------------------------------
    for ss in sss:
        if stat(aa_queries_files[ss]).st_size == 0:
            continue
        changed_blast_1 = run_tblastn_on_reads(
            se_fastq_files, pe_fastq_files, aa_queries_files[ss], tblastn,
            blast_1_evalue, blast_1_max_hsps, blast_1_qcov_hsp_perc,
            blast_1_best_hit_overhang, blast_1_best_hit_score_edge,
            blast_1_max_target_seqs, dir_prj_blast_results_fa_trim,
            pe_blast_results_file_patterns, ss, THREADS, seqtk, vsearch,
            dir_cache_prj, linfo)

        if changed_blast_1 is True:
            if ope(dir_prj_vsearch_results_fa_trim):
                rmtree(dir_prj_vsearch_results_fa_trim)
            if ope(dir_prj_spades_assemblies):
                rmtree(dir_prj_spades_assemblies)
            if ope(dir_prj_blast_assmbl):
                rmtree(dir_prj_blast_assmbl)
            if ope(dir_prj_assmbl_blast_results):
                rmtree(dir_prj_assmbl_blast_results)
            if ope(dir_prj_transcripts):
                rmtree(dir_prj_transcripts)
            if ope(dir_prj_transcripts_combined):
                rmtree(dir_prj_transcripts_combined)

    prepare_output_directories(dir_out, prj_name)

    # Run vsearch on reads ---------------------------------------------------
    for ss in sss:
        if stat(aa_queries_files[ss]).st_size == 0:
            continue
        run_vsearch_on_reads(se_fastq_files, pe_fastq_files, vsearch,
                             dir_prj_vsearch_results_fa_trim,
                             pe_vsearch_results_file_patterns, ss, seqtk,
                             linfo)

    # Run SPAdes -------------------------------------------------------------
    for ss in sss:
        if stat(aa_queries_files[ss]).st_size == 0:
            for se in se_fastq_files:
                se_fastq_files[se]['spades_assembly' + '__' + ss] = None
            for pe in pe_fastq_files:
                pe_fastq_files[pe]['spades_assembly' + '__' + ss] = None
            continue
        run_spades(se_fastq_files, pe_fastq_files, dir_prj_spades_assemblies,
                   spades, dir_temp, ss, THREADS, RAM, linfo)

    # Combine SPAdes and user provided assemblies ----------------------------
    assemblies = combine_assemblies(se_fastq_files, pe_fastq_files,
                                    user_assemblies, tax, sss)

    # Run makeblastdb on assemblies  -----------------------------------------
    makeblastdb_assemblies(assemblies, dir_prj_blast_assmbl, makeblastdb,
                           linfo)

    if any_queries is False:
        linfo(CONSRED + 'No query sequences were provided' + CONSDFL)

    # Run tblastn on assemblies ----------------------------------------------
    for ss in sss:

        if stat(aa_queries_files[ss]).st_size == 0:
            continue

        blast_2_evalue_ss = sss[ss]['blast_2_evalue']
        blast_2_max_hsps_ss = sss[ss]['blast_2_max_hsps']
        blast_2_qcov_hsp_perc_ss = sss[ss]['blast_2_qcov_hsp_perc']
        blast_2_best_hit_overhang_ss = sss[ss]['blast_2_best_hit_overhang']
        blast_2_best_hit_score_edge_ss = sss[ss]['blast_2_best_hit_score_edge']
        blast_2_max_target_seqs_ss = sss[ss]['blast_2_max_target_seqs']

        if blast_2_evalue_ss is None:
            blast_2_evalue_ss = blast_2_evalue
        if blast_2_max_hsps_ss is None:
            blast_2_max_hsps_ss = blast_2_max_hsps
        if blast_2_qcov_hsp_perc_ss is None:
            blast_2_qcov_hsp_perc_ss = blast_2_qcov_hsp_perc
        if blast_2_best_hit_overhang_ss is None:
            blast_2_best_hit_overhang_ss = blast_2_best_hit_overhang
        if blast_2_best_hit_score_edge_ss is None:
            blast_2_best_hit_score_edge_ss = blast_2_best_hit_score_edge
        if blast_2_max_target_seqs_ss is None:
            blast_2_max_target_seqs_ss = blast_2_max_target_seqs

        run_tblastn_on_assemblies(ss, assemblies, aa_queries_files[ss],
                                  tblastn, dir_prj_assmbl_blast_results,
                                  blast_2_evalue_ss, blast_2_max_hsps_ss,
                                  blast_2_qcov_hsp_perc_ss,
                                  blast_2_best_hit_overhang_ss,
                                  blast_2_best_hit_score_edge_ss,
                                  blast_2_max_target_seqs_ss, THREADS,
                                  dir_cache_prj, dir_prj_ips, linfo)

    # Prepare BLAST hits for analysis: find ORFs, translate ------------------
    for ss in sss:

        if stat(aa_queries_files[ss]).st_size == 0:
            continue

        min_target_orf_len_ss = sss[ss]['min_target_orf_length']
        max_target_orf_len_ss = sss[ss]['max_target_orf_length']
        organelle = sss[ss]['organelle']

        blast_2_qcov_hsp_perc_ss = sss[ss]['blast_2_qcov_hsp_perc']

        if blast_2_qcov_hsp_perc_ss is None:
            blast_2_qcov_hsp_perc_ss = blast_2_qcov_hsp_perc

        find_orfs_translate(ss, assemblies, dir_prj_transcripts, seqtk,
                            dir_temp, prepend_assmbl, min_target_orf_len_ss,
                            max_target_orf_len_ss, allow_non_aug,
                            allow_no_strt_cod,
                            allow_no_stop_cod, tax, tax_group, tax_ids_user,
                            blast_2_qcov_hsp_perc_ss, organelle, linfo)

    # GFF3 files from kakapo results JSON files ------------------------------
    for ss in sss:
        if stat(aa_queries_files[ss]).st_size == 0:
            continue
        gff_from_json(ss, assemblies, dir_prj_ips,
                      dir_prj_transcripts_combined, prj_name, linfo)

    # Run InterProScan 5 -----------------------------------------------------
    if should_run_ipr is True:

        ss_names = tuple(sss.keys())

        # Determine the length of printed strings, for better spacing --------
        max_title_a_len = 0
        max_run_id_len = 0
        for a in assemblies:
            for ss in ss_names:
                if 'transcripts_aa_orf_fasta_file__' + ss not in a:
                    continue

                aa_file = a['transcripts_aa_orf_fasta_file__' + ss]

                if aa_file is None:
                    continue

                assmbl_name = a['name']
                run_id = ss + '_' + assmbl_name
                max_run_id_len = max(len(run_id), max_run_id_len)

                seqs = seq_records_to_dict(read_fasta(aa_file, SEQ_TYPE_AA))

                # Filter all ORFs except the first one.
                for seq_def in tuple(seqs.keys()):
                    seq_def_prefix = seq_def.split(' ')[0]
                    if seq_def_prefix.endswith('ORF001'):
                        max_title_a_len = max(len(seq_def_prefix),
                                              max_title_a_len)

        max_title_a_len += 2
        max_run_id_len += 2
        # --------------------------------------------------------------------

        # Parallel version BEGIN ---------------------------------------------
        parallel_run_count = min(THREADS, len(ss_names))

        def run_inter_pro_scan_parallel(ss):

            if stat(aa_queries_files[ss]).st_size == 0:
                return

            run_inter_pro_scan(ss, assemblies, email, dir_prj_ips,
                               dir_cache_prj, parallel_run_count,
                               max_title_a_len, max_run_id_len, linfo)

            # GFF3 files from kakapo and InterProScan 5 results JSON files
            gff_from_json(ss, assemblies, dir_prj_ips,
                          dir_prj_transcripts_combined, prj_name, linfo)

        Parallel(n_jobs=parallel_run_count, verbose=0, require='sharedmem')(
            delayed(run_inter_pro_scan_parallel)(ss) for ss in ss_names)
        # Parallel version END -----------------------------------------------

        # Sequential version BEGIN -------------------------------------------
        # for ss in sss:
        #     if stat(aa_queries_files[ss]).st_size == 0:
        #         continue
        #     run_inter_pro_scan(ss, assemblies, email, dir_prj_ips,
        #                        dir_cache_prj, 1, max_title_a_len,
        #                        max_run_id_len, linfo)
        #     # GFF3 files from kakapo and InterProScan 5 results JSON files
        #     gff_from_json(ss, assemblies, dir_prj_ips,
        #                   dir_prj_transcripts_combined, prj_name, linfo)
        # Sequential version END ---------------------------------------------

    # Download CDS for NCBI protein queries ----------------------------------

    prot_cds_ncbi_files = OrderedDict()

    # Parallel version BEGIN -------------------------------------------------
    def dnld_cds_for_ncbi_prot_acc_parallel(ss):

        if stat(aa_queries_files[ss]).st_size == 0:
            return

        if ss not in prot_acc_user_filtered:
            return

        prot_cds_ncbi_files[ss] = opj(dir_prj_transcripts_combined, prj_name +
                                      '_ncbi_query_cds__' + ss + '.fasta')

        if len(prot_acc_user_filtered[ss]) > 0:
            dnld_cds_for_ncbi_prot_acc(ss, prot_acc_user_filtered[ss],
                                       prot_cds_ncbi_files[ss], tax,
                                       dir_cache_prj, linfo)

    ss_names = tuple(sss.keys())
    Parallel(n_jobs=2, verbose=0, require='sharedmem')(
        delayed(dnld_cds_for_ncbi_prot_acc_parallel)(ss) for ss in ss_names)
    # Parallel version END ---------------------------------------------------

    # Sequential version BEGIN -----------------------------------------------
    # for ss in sss:
    #     if stat(aa_queries_files[ss]).st_size == 0:
    #         continue
    #     if ss not in prot_acc_user_filtered:
    #         continue
    #     prot_cds_ncbi_files[ss] = opj(dir_prj_transcripts_combined, prj_name +
    #                                   '_ncbi_query_cds__' + ss + '.fasta')
    #     if len(prot_acc_user_filtered[ss]) > 0:
    #         dnld_cds_for_ncbi_prot_acc(ss, prot_acc_user_filtered[ss],
    #                                    prot_cds_ncbi_files[ss], tax,
    #                                    dir_cache_prj, linfo)
    # Sequential version END -------------------------------------------------

    # ------------------------------------------------------------------------

    rmtree(dir_temp)
    log, _ = prepare_logger(console=False)

    # ------------------------------------------------------------------------

    rerun = input('\nRepeat ([y]/n)? ').lower().strip()
    if rerun.startswith('y') or rerun == '':
        print()
        return False
    else:
        print('\nExiting...')
        return True

    # ------------------------------------------------------------------------


def run_kakapo():
    while True:
        stop = main()
        if stop is True:
            break


if __name__ == '__main__':
    run_kakapo()
