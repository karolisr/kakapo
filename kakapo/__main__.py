#!/usr/bin/env python3

"""Kakapo main file."""

import os

##############################################################################
# Uncomment to allow running this file without installing kakapo with pip ####
##############################################################################
# import inspect
# import sys
# _ = inspect.currentframe()
# assert _ is not None
# SCRIPT_FILE_PATH = inspect.getfile(_)
# SCRIPT_DIR_PATH = os.path.dirname(os.path.abspath(SCRIPT_FILE_PATH))
# KAKAPO_DIR_PATH = os.path.sep.join(SCRIPT_DIR_PATH.split(os.path.sep)[0:-1])
# sys.path.insert(0, KAKAPO_DIR_PATH)
##############################################################################

import argparse
import logging

from collections import OrderedDict
from io import StringIO
from os import stat
from os.path import exists as ope
from os.path import join as opj
from shutil import copyfile
from shutil import rmtree
from sys import exit

from joblib import Parallel, delayed

from ncbi_taxonomy_local import Taxonomy, TaxonomyRAM, TaxonomySQL

from . import __version__ as KKPV
from . import __script_name__ as KKP

from .utils import dependencies as deps
from .utils.logging import Log
from .utils.misc import make_dirs
from .utils.misc import time_stamp

from .tools.config import CONSRED, CONSDFL, CONBLUE
from .tools.config import DIR_DAT, DIR_DEP, DIR_TAX, DIR_KRK
from .tools.config import MACHINE_TYPE
from .tools.config import OS_ID, DIST_ID, DEBIAN_DISTS, REDHAT_DISTS
from .tools.config import RELEASE_NAME
from .tools.config import SCRIPT_INFO
from .tools.config import NCPU, RAM
from .tools.config_file_parse import config_file_parse
from .tools.config_file_parse import ss_file_parse
from .tools.config_file_parse import use_colors

from .tools.seq import SEQ_TYPE_AA

from .tools.bioio import read_fasta
from .tools.bioio import seq_records_to_dict
from .tools.transl_tables import TranslationTable

from .flow.a_prepare import prepare_output_directories

from .flow.b_process_reads import dnld_sra_fastq_files
from .flow.b_process_reads import dnld_sra_info
from .flow.b_process_reads import file_name_patterns
from .flow.b_process_reads import filtered_fq_to_fa
from .flow.b_process_reads import makeblastdb_fq
from .flow.b_process_reads import min_accept_read_len
from .flow.b_process_reads import run_bt2_fq
from .flow.b_process_reads import run_kraken2
from .flow.b_process_reads import run_rcorrector
from .flow.b_process_reads import run_trimmomatic
from .flow.b_process_reads import user_fastq_files

from .flow.c_process_queries import combine_aa_fasta
from .flow.c_process_queries import dnld_pfam_uniprot_seqs
from .flow.c_process_queries import dnld_prot_seqs
from .flow.c_process_queries import filter_queries
from .flow.c_process_queries import pfam_uniprot_accessions
from .flow.c_process_queries import user_aa_fasta
from .flow.c_process_queries import user_entrez_search
from .flow.c_process_queries import user_protein_accessions

from .flow.d_search_reads import run_tblastn_on_reads
from .flow.d_search_reads import run_vsearch_on_reads

from .flow.e_process_assmbl import combine_assemblies
from .flow.e_process_assmbl import makeblastdb_assemblies
from .flow.e_process_assmbl import run_spades

from .flow.f_search_assmbl import run_tblastn_on_assemblies
from .flow.g_find_orfs import find_orfs_translate
from .flow.h_prepare_gff import gff_from_json
from .flow.i_inter_pro_scan import run_inter_pro_scan
from .flow.j_dnld_aa_query_cds import dnld_cds_for_ncbi_prot_acc

# Silence urllib3 overly verbose logging.
logging.getLogger('urllib3').setLevel(logging.WARNING)

# Command line arguments -----------------------------------------------------
USAGE = CONBLUE + '{} --cfg project_configuration_file ' \
    '--ss search_strategies_file'.format(KKP) + CONSDFL

PARSER = argparse.ArgumentParser(
    prog=KKP,
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
    help='Path to a {} project configuration file.'.format(KKP))

PARSER.add_argument(
    '--ss',
    type=str,
    required=False,
    metavar='path',
    dest='SS_FILE_PATH',
    help='Path to a {} search strategies file.'.format(KKP))

PARSER.add_argument(
    '--ncpu',
    type=int,
    metavar='count',
    required=False,
    dest='NCPU_ARG',
    help='Number of CPUs to use.')

PARSER.add_argument(
    '--stop-after-filter',
    action='store_true',
    required=False,
    dest='STOP_AFTER_FILTER',
    help='Stop {} after Kraken2/Bowtie2 filtering '
         'step.'.format(KKP))

PARSER.add_argument(
    '--force-deps',
    action='store_true',
    default=False,
    required=False,
    dest='FORCE_DEPS',
    help='Force the use of {}-installed dependencies,\neven if they are '
         'already available on the system.'.format(KKP))

PARSER.add_argument(
    '--install-deps',
    action='store_true',
    default=False,
    required=False,
    dest='INSTALL_DEPS',
    help='Install {} dependencies and quit.'.format(KKP))

PARSER.add_argument(
    '--dnld-kraken-dbs',
    action='store_true',
    default=False,
    required=False,
    dest='DNLD_KRAKEN_DBS',
    help='Download Kraken2 databases and quit.')

PARSER.add_argument(
    '--clean-data-dir',
    action='store_true',
    required=False,
    dest='CLEAN_DATA_DIR',
    help='Remove cached NCBI taxonomy data and all software\ndependencies '
         'downloaded by {}.'.format(KKP))

PARSER.add_argument(
    '-v', '--version',
    action='store_true',
    required=False,
    dest='PRINT_VERSION',
    help='Print {} version.'.format(KKP))

PARSER.add_argument(
    '-h', '--help',
    action='store_true',
    required=False,
    dest='PRINT_HELP',
    help='Print {} help information.'.format(KKP))

ARGS = PARSER.parse_args()

CLEAN_DATA_DIR = ARGS.CLEAN_DATA_DIR
CONFIG_FILE_PATH = ARGS.CONFIG_FILE_PATH
SS_FILE_PATH = ARGS.SS_FILE_PATH
NCPU_ARG = ARGS.NCPU_ARG
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
    print(KKP + ' v' + KKPV)
    exit(0)

if CLEAN_DATA_DIR is True and ope(DIR_DAT):
    print(CONSRED + 'Removing ' + KKP + ' data directory: '
          + CONSDFL + DIR_DAT)
    try:
        rmtree(DIR_DAT, ignore_errors=True)
    except Exception:
        try:
            rmtree(DIR_DAT, ignore_errors=True)
        except Exception:
            pass
    exit(0)
elif CLEAN_DATA_DIR is True:
    print(CONSRED + 'The data directory does not exist. Nothing to do.'
          + CONSDFL)
    exit(0)

if CLEAN_DATA_DIR is False and CONFIG_FILE_PATH is not None:
    if not ope(CONFIG_FILE_PATH):
        print(CONSRED + 'Configuration file ' + CONFIG_FILE_PATH
              + ' does not exist.' + CONSDFL)
        exit(0)
elif INSTALL_DEPS is True or DNLD_KRAKEN_DBS is True:
    pass
else:
    print(SCRIPT_INFO)
    print(CONSRED + 'Configuration file was not provided. Nothing to do.'
          + CONSDFL)
    print()
    PARSER.print_help()
    exit(0)

if CLEAN_DATA_DIR is False and SS_FILE_PATH is not None:
    if not ope(SS_FILE_PATH):
        print(CONSRED + 'Search strategies file ' + SS_FILE_PATH
              + ' does not exist.' + CONSDFL)
        exit(0)
elif INSTALL_DEPS is True or DNLD_KRAKEN_DBS is True or SS_FILE_PATH is None:
    pass
else:
    pass

print(SCRIPT_INFO)

NCPU = NCPU - 1
if NCPU_ARG is not None:
    if NCPU_ARG > NCPU:
        print(CONSRED + 'The number of CPU cores requested is larger '
              'than is available on the system (minus one). Will use '
              + str(NCPU) + ' instead.' + CONSDFL)
        NCPU_ARG = NCPU
    else:
        NCPU = NCPU_ARG

COLORS = False
if CONFIG_FILE_PATH is not None:
    COLORS = use_colors(CONFIG_FILE_PATH)


def main():
    """Run the script."""
    # Prepare initial logger (before we know the log file path) --------------
    prj_log_file_suffix = time_stamp()
    log_stream = StringIO()

    Log.set_colors(COLORS)
    Log.set_file(log_stream)
    Log.set_write(True)

    # Prepare kakapo data directory ----------------------------------------
    if ope(DIR_DAT):
        Log.log(inf='Found kakapo data directory:', s=DIR_DAT)
    else:
        Log.wrn('Creating kakapo data directory:', DIR_DAT)
        make_dirs(DIR_DAT)

    print()

    # Check for dependencies -------------------------------------------------
    Log.inf('Checking for dependencies.')
    make_dirs(DIR_DEP)
    make_dirs(DIR_KRK)

    gzip = deps.dep_check_gzip()
    pigz = deps.dep_check_pigz()
    if (gzip is None) and (pigz is None):
        Log.err('Could not find either "gzip" or "pigz". Cannot continue.', '')
        exit(0)
    if pigz is None:
        Log.wrn('Will be using slower "gzip" program for compression. '
                'For much faster compression, please install "pigz".', '')

    java = deps.dep_check_java()
    if (java is None):
        Log.err('Could not find "Java". Cannot continue.', '')
        exit(0)

    seqtk = deps.dep_check_seqtk(DIR_DEP, FORCE_DEPS)
    trimmomatic, adapters = deps.dep_check_trimmomatic(DIR_DEP)
    fasterq_dump = deps.dep_check_sra_toolkit(DIR_DEP, OS_ID, DIST_ID,
                                              DEBIAN_DISTS, REDHAT_DISTS,
                                              FORCE_DEPS)
    makeblastdb, _, tblastn = deps.dep_check_blast(DIR_DEP, OS_ID, DIST_ID,
                                                   DEBIAN_DISTS, REDHAT_DISTS,
                                                   FORCE_DEPS)
    vsearch = deps.dep_check_vsearch(DIR_DEP, OS_ID, DIST_ID, MACHINE_TYPE, DEBIAN_DISTS,
                                     REDHAT_DISTS, FORCE_DEPS)
    spades = deps.dep_check_spades(DIR_DEP, OS_ID, FORCE_DEPS)
    bowtie2, bowtie2_build = deps.dep_check_bowtie2(DIR_DEP, OS_ID, FORCE_DEPS)
    rcorrector = deps.dep_check_rcorrector(DIR_DEP, FORCE_DEPS)
    kraken2, kraken2_build = deps.dep_check_kraken2(DIR_DEP, OS_ID,
                                                    RELEASE_NAME, MACHINE_TYPE,
                                                    FORCE_DEPS)

    kakapolib = deps.dep_check_kakapolib(FORCE_DEPS)
    if kakapolib is None:
        Log.err('Could not compile "kakapolib". Cannot continue.', '')
        exit(0)

    print()

    kraken2_dbs = deps.dnld_kraken2_dbs(DIR_KRK)

    if INSTALL_DEPS is True or DNLD_KRAKEN_DBS is True:
        exit(0)

    print()

    # Initialize NCBI taxonomy database --------------------------------------
    # tax = TaxonomyRAM(data_dir=DIR_TAX, logger=Log, check_for_updates=True)
    tax = TaxonomySQL(backend='SQLite', data_dir=DIR_TAX, logger=Log, check_for_updates=True)

    # Parse configuration file -----------------------------------------------
    Log.log(inf='Reading project configuration file:', s=CONFIG_FILE_PATH)
    _ = config_file_parse(CONFIG_FILE_PATH, tax)

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
    requery_after = _['requery_after']
    fq_pe = _['fq_pe']
    fq_se = _['fq_se']
    should_run_rcorrector = _['should_run_rcorrector']
    rcorrector_before_trimmomatic = _['rcorrector_before_trimmomatic']
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

    os.environ['ENTREZ_KEY'] = _['entrez_api_key']

    taxid_for_queries = tax_group

    print()

    # Parse search strategies file -------------------------------------------
    if SS_FILE_PATH is not None:
        Log.log(inf='Reading search strategies file:', s=SS_FILE_PATH)
        sss = ss_file_parse(SS_FILE_PATH)
    else:
        Log.wrn('Search strategies file was not provided.\n'
                'Will process reads, assemblies and then stop.', '')
        sss = dict()

    print()

    # Create output directory ------------------------------------------------
    if dir_out is not None:
        if ope(dir_out):
            Log.log(inf='Found output directory:', s=dir_out)
        else:
            Log.wrn('Creating output directory:', dir_out)
            make_dirs(dir_out)

    print()

    # Write Kakapo version information to the output directory ---------------
    version_file = opj(dir_out, 'kakapo_version.txt')
    if ope(version_file):
        with open(version_file, 'r') as f:
            version_prev = f.read().strip()
            if KKPV != version_prev:
                Log.wrn('The output directory contains data produced by a '
                        + 'different version of Kakapo: ' + version_prev
                        + '.\nThe currently running version is: ' + KKPV
                        + '.\n'
                        + 'Delete "kakapo_version.txt" file located in the '
                        + 'output directory if you would like to continue.',
                        '')
                exit(0)

    with open(version_file, 'w') as f:
        f.write(KKPV)

    # Create subdirectories in the output directory --------------------------
    _ = prepare_output_directories(dir_out, prj_name,
                                   rcorrector_before_trimmomatic)

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

    globals()['KKP_DIR_TMP'] = dir_temp

    # Archive the configuration and search strategies files used -------------
    run_cfg_file = opj(dir_prj_logs, prj_name + '_' + prj_log_file_suffix
                       + '.cfg.ini')
    copyfile(CONFIG_FILE_PATH, run_cfg_file)

    if SS_FILE_PATH is not None:
        run_ss_file = opj(dir_prj_logs, prj_name + '_' + prj_log_file_suffix
                          + '.ss.ini')
        copyfile(SS_FILE_PATH, run_ss_file)

    # Prepare logger ---------------------------------------------------------
    prj_log_file = opj(dir_prj_logs, prj_name + '_' + prj_log_file_suffix
                       + '.log')
    with open(prj_log_file, 'w') as f:
        f.write(SCRIPT_INFO.strip() + '\n\n' + log_stream.getvalue())

    Log.set_colors(COLORS)
    Log.set_file(prj_log_file)
    Log.set_write(True)

    log_stream.close()

    # Download SRA run metadata if needed ------------------------------------
    sra_runs_info, sras_acceptable = dnld_sra_info(sras, dir_cache_prj)

    # ------------------------------------------------------------------------
    tax_ids_user |= set([int(sra_runs_info[x]['TaxID']) for x in sras_acceptable])
    _ = tax.common_taxid(tax_ids_user)
    _ = tax.higher_rank_for_taxid(_, 'order')
    if _ != '':
        taxid_for_queries = tax.taxid_for_name_and_group_taxid(_, tax_group)

    # Download SRA run FASTQ files if needed ---------------------------------
    se_fastq_files_sra, pe_fastq_files_sra = dnld_sra_fastq_files(
        sras_acceptable, sra_runs_info, dir_fq_data, fasterq_dump, NCPU,
        dir_temp)

    # User provided FASTQ files ----------------------------------------------
    se_fastq_files_usr, pe_fastq_files_usr = user_fastq_files(fq_se, fq_pe)

    # Collate FASTQ file info ------------------------------------------------
    se_fastq_files = se_fastq_files_sra.copy()
    se_fastq_files.update(se_fastq_files_usr)
    pe_fastq_files = pe_fastq_files_sra.copy()
    pe_fastq_files.update(pe_fastq_files_usr)

    def gc_tt(k, d_local, tax_local: Taxonomy):
        taxid = d_local[k]['tax_id']

        gc = tax_local.genetic_code_for_taxid(taxid)

        d_local[k]['gc_id'] = gc
        d_local[k]['gc_tt'] = TranslationTable(gc)

        gc_mito = None
        tt_mito = None

        gc_plastid = None
        tt_plastid = None

        if tax_local.is_eukaryote(taxid) is True:
            gc_mito = tax_local.mito_genetic_code_for_taxid(taxid)
            if gc_mito != 0:
                tt_mito = TranslationTable(gc_mito)

            if tax_local.contains_plastid(taxid) is True:
                gc_plastid = tax_local.plastid_genetic_code_for_taxid(taxid)
                if gc_plastid != 0:
                    tt_plastid = TranslationTable(gc_plastid)

        d_local[k]['gc_id_mito'] = gc_mito
        d_local[k]['gc_tt_mito'] = tt_mito

        d_local[k]['gc_id_plastid'] = gc_plastid
        d_local[k]['gc_tt_plastid'] = tt_plastid

    for se in se_fastq_files:
        gc_tt(se, se_fastq_files, tax)

    for pe in pe_fastq_files:
        gc_tt(pe, pe_fastq_files, tax)

    # Minimum acceptable read length -----------------------------------------
    min_accept_read_len(se_fastq_files, pe_fastq_files, dir_temp,
                        dir_cache_fq_minlen)

    # File name patterns -----------------------------------------------------
    a, b, c, d, e = file_name_patterns()

    pe_trim_fq_file_patterns = a
    pe_trim_fa_file_patterns = b
    pe_blast_db_file_patterns = c
    pe_blast_results_file_patterns = d
    pe_vsearch_results_file_patterns = e

    if rcorrector_before_trimmomatic is True:
        # Run Rcorrector -----------------------------------------------------
        run_rcorrector(se_fastq_files, pe_fastq_files, dir_fq_cor_data,
                       rcorrector, NCPU, dir_temp, pe_trim_fq_file_patterns,
                       should_run_rcorrector, rcorrector_before_trimmomatic)

        # Run Trimmomatic ----------------------------------------------------
        run_trimmomatic(se_fastq_files, pe_fastq_files, dir_fq_trim_data,
                        trimmomatic, adapters, pe_trim_fq_file_patterns,
                        NCPU, rcorrector_before_trimmomatic)
    else:
        # Run Trimmomatic ----------------------------------------------------
        run_trimmomatic(se_fastq_files, pe_fastq_files, dir_fq_trim_data,
                        trimmomatic, adapters, pe_trim_fq_file_patterns,
                        NCPU, rcorrector_before_trimmomatic)

        # Run Rcorrector -----------------------------------------------------
        run_rcorrector(se_fastq_files, pe_fastq_files, dir_fq_cor_data,
                       rcorrector, NCPU, dir_temp, pe_trim_fq_file_patterns,
                       should_run_rcorrector, rcorrector_before_trimmomatic)

    # Run Bowtie 2 -----------------------------------------------------------
    run_bt2_fq(se_fastq_files, pe_fastq_files, dir_fq_filter_bt2_data,
               bowtie2, bowtie2_build, NCPU, dir_temp, bt2_order,
               pe_trim_fq_file_patterns, tax, dir_cache_refseqs,
               rcorrector_before_trimmomatic)

    # Run Kraken 2 -----------------------------------------------------------
    run_kraken2(krkn_order, kraken2_dbs, se_fastq_files, pe_fastq_files,
                dir_fq_filter_krkn2_data, kraken_confidence, kraken2, NCPU,
                dir_temp, pe_trim_fq_file_patterns, gzip, pigz)

    se_fastq_files = OrderedDict(se_fastq_files)
    pe_fastq_files = OrderedDict(pe_fastq_files)

    se_fastq_files = OrderedDict(sorted(se_fastq_files.items(),
                                        key=lambda x: x[1]['filter_path_fq']))

    pe_fastq_files = OrderedDict(sorted(pe_fastq_files.items(),
                                        key=lambda x: x[1]['filter_path_fq']))

    # Stop after filter ------------------------------------------------------
    if STOP_AFTER_FILTER is True:
        print()
        Log.wrn('Stopping after Kraken 2 / Bowtie 2 filtering step as requested.', '')
        exit(0)

    # Convert filtered FASTQ files to FASTA ----------------------------------
    filtered_fq_to_fa(se_fastq_files, pe_fastq_files, dir_fa_trim_data, seqtk,
                      pe_trim_fa_file_patterns, NCPU, gzip, pigz)

    # Run makeblastdb on reads -----------------------------------------------
    if len(se_fastq_files) > 0 or len(pe_fastq_files) > 0:
        print()
        Log.inf('Building BLAST databases for reads.')
        if makeblastdb is None:
            Log.err('makeblastdb is not available. Cannot continue. Exiting.', '')
            exit(0)

    parallel_run_count = NCPU

    def makeblastdb_fq_se_parallel(se_fastq_files_local):
        makeblastdb_fq(se_fastq_files_local, dict(), dir_blast_fa_trim,
                       makeblastdb, pe_blast_db_file_patterns)

        se_fastq_files.update(se_fastq_files_local)

    if parallel_run_count > 0:
        _ = se_fastq_files.keys()
        Parallel(n_jobs=parallel_run_count, verbose=0, require='sharedmem')(
            delayed(makeblastdb_fq_se_parallel)({k: se_fastq_files[k]}) for k in _)

    def makeblastdb_fq_pe_parallel(pe_fastq_files_local):
        makeblastdb_fq(dict(), pe_fastq_files_local, dir_blast_fa_trim,
                       makeblastdb, pe_blast_db_file_patterns)

        pe_fastq_files.update(pe_fastq_files_local)

    if parallel_run_count > 0:
        _ = pe_fastq_files.keys()
        Parallel(n_jobs=parallel_run_count, verbose=0, require='sharedmem')(
            delayed(makeblastdb_fq_pe_parallel)({k: pe_fastq_files[k]}) for k in _)

    if len(sss) > 0 and (len(se_fastq_files) > 0 or len(pe_fastq_files) > 0):
        print()

    # Pfam uniprot accessions ------------------------------------------------
    pfam_uniprot_acc = OrderedDict()
    for ss in sss:
        pfam_acc = sss[ss]['pfam_families']
        pfam_uniprot_acc[ss] = pfam_uniprot_accessions(ss, pfam_acc,
                                                       [taxid_for_queries],
                                                       dir_cache_pfam_acc)

    # Download Pfam uniprot sequences if needed ------------------------------
    aa_uniprot_files = OrderedDict()
    for ss in sss:
        aa_uniprot_files[ss] = opj(dir_prj_queries, 'aa_uniprot__' + ss
                                   + '.fasta')
        dnld_pfam_uniprot_seqs(ss, pfam_uniprot_acc[ss], aa_uniprot_files[ss],
                               dir_cache_prj)

    # User provided entrez query ---------------------------------------------
    prot_acc_user_from_query = OrderedDict()
    for ss in sss:
        entrez_queries = sss[ss]['entrez_search_queries']
        prot_acc_user_from_query[ss] = user_entrez_search(ss, entrez_queries,
                                                          dir_cache_prj,
                                                          requery_after)

    print()

    # User provided protein accessions ---------------------------------------
    prot_acc_user = OrderedDict()
    for ss in sss:
        prot_acc_all = sorted(set(sss[ss]['ncbi_accessions_aa']
                                  + prot_acc_user_from_query[ss]))
        prot_acc_user[ss] = user_protein_accessions(ss, prot_acc_all,
                                                    dir_cache_prj, tax)

    # Download from NCBI if needed -------------------------------------------
    aa_prot_ncbi_files = OrderedDict()
    for ss in sss:
        aa_prot_ncbi_files[ss] = opj(dir_prj_queries, 'aa_prot_ncbi__' + ss
                                     + '.fasta')
        prot_acc_user[ss] = dnld_prot_seqs(ss, prot_acc_user[ss],
                                           aa_prot_ncbi_files[ss])

    print()

    # User provided protein sequences ----------------------------------------
    aa_prot_user_files = OrderedDict()
    for ss in sss:
        user_queries = sss[ss]['fasta_files_aa']
        aa_prot_user_files[ss] = opj(dir_prj_queries, 'aa_prot_user__' + ss
                                     + '.fasta')
        user_aa_fasta(ss, user_queries, aa_prot_user_files[ss])
        print()

    # Combine all AA queries -------------------------------------------------
    aa_queries_files = OrderedDict()
    for ss in sss:
        aa_queries_files[ss] = opj(dir_prj_queries, 'aa_all__' + ss + '.fasta')
        combine_aa_fasta(ss, [aa_uniprot_files[ss], aa_prot_ncbi_files[ss],
                              aa_prot_user_files[ss]], aa_queries_files[ss])

    # Filter AA queries ------------------------------------------------------
    prot_acc_user_filtered = OrderedDict()
    for ss in sss:
        min_query_length = sss[ss]['min_query_length']
        max_query_length = sss[ss]['max_query_length']
        max_query_identity = sss[ss]['max_query_identity']

        # Dereplicate all queries
        filter_queries(ss, aa_queries_files[ss], min_query_length,
                       max_query_length, max_query_identity,
                       vsearch, prot_acc_user[ss], overwrite=True)

        # Dereplicate only NCBI queries. CDS for these will be downloaded
        # later for reference.
        if ope(aa_prot_ncbi_files[ss]):
            prot_acc_user_filtered[ss] = filter_queries(
                ss, aa_prot_ncbi_files[ss], min_query_length, max_query_length,
                max_query_identity, vsearch, prot_acc_user[ss],
                overwrite=False, logging=False)

    # Check if there are any query sequences ---------------------------------
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

        ss_organelle = sss[ss]['organelle']
        changed_blast_1 = run_tblastn_on_reads(
            se_fastq_files, pe_fastq_files, aa_queries_files[ss], tblastn,
            blast_1_evalue, blast_1_max_hsps, blast_1_qcov_hsp_perc,
            blast_1_best_hit_overhang, blast_1_best_hit_score_edge,
            blast_1_max_target_seqs, dir_prj_blast_results_fa_trim,
            pe_blast_results_file_patterns, ss, NCPU, seqtk, vsearch,
            dir_cache_prj, ss_organelle)

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

    prepare_output_directories(dir_out, prj_name,
                               rcorrector_before_trimmomatic)

    # Run vsearch on reads ---------------------------------------------------
    for ss in sss:
        if stat(aa_queries_files[ss]).st_size == 0:
            continue
        print()
        Log.log(inf='Checking if Vsearch should be run:', s=ss)
        run_vsearch_on_reads(se_fastq_files, pe_fastq_files, vsearch,
                             dir_prj_vsearch_results_fa_trim,
                             pe_vsearch_results_file_patterns, ss, seqtk)

    # Run SPAdes -------------------------------------------------------------
    for ss in sss:
        if stat(aa_queries_files[ss]).st_size == 0:
            for se in se_fastq_files:
                se_fastq_files[se]['spades_assembly' + '__' + ss] = None
            for pe in pe_fastq_files:
                pe_fastq_files[pe]['spades_assembly' + '__' + ss] = None
            continue
        print()
        Log.log(inf='Checking if SPAdes should be run:', s=ss)
        run_spades(se_fastq_files, pe_fastq_files, dir_prj_spades_assemblies,
                   spades, dir_temp, ss, NCPU, RAM)

    # Combine SPAdes and user provided assemblies ----------------------------
    assemblies = combine_assemblies(se_fastq_files, pe_fastq_files,
                                    user_assemblies, tax, sss)

    # Run makeblastdb on assemblies  -----------------------------------------
    makeblastdb_assemblies(assemblies, dir_prj_blast_assmbl, makeblastdb)

    if any_queries is False:
        print()
        Log.wrn('No query sequences were provided.', '')

    # Run tblastn on assemblies ----------------------------------------------
    for ss in sss:

        if stat(aa_queries_files[ss]).st_size == 0:
            continue

        should_run_tblastn = False
        for a in assemblies:
            assmbl_src = a['src']
            assmbl_name = a['name']
            if assmbl_src != 'user_fasta':
                if assmbl_name.endswith('__' + ss):
                    should_run_tblastn = True
                    break
            else:
                should_run_tblastn = True
                break

        if should_run_tblastn is False:
            print()
            Log.log(inf='Will not run BLAST. No transcripts exist:', s=ss)
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
                                  blast_2_max_target_seqs_ss, NCPU,
                                  dir_cache_prj, dir_prj_ips)

    # Prepare BLAST hits for analysis: find ORFs, translate ------------------
    for ss in sss:

        if stat(aa_queries_files[ss]).st_size == 0:
            continue

        min_target_orf_len_ss = sss[ss]['min_target_orf_length']
        max_target_orf_len_ss = sss[ss]['max_target_orf_length']
        ss_organelle = sss[ss]['organelle']

        blast_2_qcov_hsp_perc_ss = sss[ss]['blast_2_qcov_hsp_perc']

        if blast_2_qcov_hsp_perc_ss is None:
            blast_2_qcov_hsp_perc_ss = blast_2_qcov_hsp_perc

        find_orfs_translate(ss, assemblies, dir_prj_transcripts, seqtk,
                            dir_temp, prepend_assmbl, min_target_orf_len_ss,
                            max_target_orf_len_ss, allow_non_aug,
                            allow_no_strt_cod,
                            allow_no_stop_cod, tax, tax_group, tax_ids_user,
                            blast_2_qcov_hsp_perc_ss, ss_organelle)

    # GFF3 files from kakapo results JSON files ------------------------------
    for ss in sss:
        if stat(aa_queries_files[ss]).st_size == 0:
            continue
        gff_from_json(ss, assemblies, dir_prj_ips,
                      dir_prj_transcripts_combined, prj_name)

    # Run InterProScan 5 -----------------------------------------------------
    if should_run_ipr is True:
        ss_names = tuple(sss.keys())
        if len(ss_names) > 0:
            print()

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

        ipr_list = list()
        for asmbl in assemblies:
            for ss in ss_names:
                if stat(aa_queries_files[ss]).st_size == 0:
                    continue
                ipr_list.append([ss, asmbl])

        parallel_run_count = min(5, len(ipr_list))

        def run_inter_pro_scan_parallel(ipr_comb):

            run_inter_pro_scan(ipr_comb[0], [ipr_comb[1], ], email, dir_prj_ips,
                               dir_cache_prj, parallel_run_count,
                               max_title_a_len, max_run_id_len)

            # GFF3 files from kakapo and InterProScan 5 results JSON files
            # gff_from_json(ipr_comb[0], [ipr_comb[1], ], dir_prj_ips,
            #               dir_prj_transcripts_combined, prj_name)

        if parallel_run_count > 0:
            Parallel(n_jobs=parallel_run_count, verbose=0, require='sharedmem')(
                delayed(run_inter_pro_scan_parallel)(ipr_comb) for ipr_comb in ipr_list)

        # --------------------------------------------------------------------

        parallel_run_count = NCPU

        def produce_final_gff(ss_local):
            if stat(aa_queries_files[ss_local]).st_size == 0:
                return

            # GFF3 files from kakapo and InterProScan 5 results JSON files
            gff_from_json(ss_local, assemblies, dir_prj_ips,
                          dir_prj_transcripts_combined, prj_name)

        if parallel_run_count > 0:
            Parallel(n_jobs=parallel_run_count, verbose=0, require='sharedmem')(
                delayed(produce_final_gff)(ss) for ss in ss_names)

    # Download CDS for NCBI protein queries ----------------------------------
    prot_cds_ncbi_files = OrderedDict()

    def dnld_cds_for_ncbi_prot_acc_parallel(ss_local):
        if stat(aa_queries_files[ss_local]).st_size == 0:
            return

        if ss_local not in prot_acc_user_filtered:
            return

        prot_cds_ncbi_files[ss_local] = opj(
            dir_prj_transcripts_combined, prj_name + '_ncbi_query_cds__'
            + ss_local + '.fasta')

        if len(prot_acc_user_filtered[ss_local]) > 0:
            dnld_cds_for_ncbi_prot_acc(ss_local, prot_acc_user_filtered[ss_local],
                                       prot_cds_ncbi_files[ss_local], tax,
                                       dir_cache_prj)

    ss_names = tuple(sss.keys())
    if len(ss_names) > 0:
        print()
    Parallel(n_jobs=2, verbose=0, require='sharedmem')(
        delayed(dnld_cds_for_ncbi_prot_acc_parallel)(ss) for ss in ss_names)

    # ------------------------------------------------------------------------
    rmtree(dir_temp)
    # ------------------------------------------------------------------------
    return True
    # ------------------------------------------------------------------------


def run_kakapo():
    try:
        while True:
            stop = main()
            if stop is True:
                break

    except KeyboardInterrupt:
        print('\r', end='')
        print('  ')
        Log.log(wrn=f'{KKP} was interrupted.', s='Good bye.', timestamp=True)

        # ToDo: cleanup. -----------------------------------------------------
        if 'KKP_DIR_TMP' in globals():
            dir_tmp = globals()['KKP_DIR_TMP']
            rmtree(dir_tmp, ignore_errors=True)
        # --------------------------------------------------------------------

        try:
            exit(130)
        except SystemExit:
            os._exit(130)


if __name__ == '__main__':
    run_kakapo()
