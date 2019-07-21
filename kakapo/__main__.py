#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""kakapo"""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import with_statement

from os.path import exists as ope
from os.path import join as opj
from shutil import rmtree
from sys import exit

from ncbi_taxonomy_local import taxonomy

from kakapo import dependencies as deps
from kakapo.config import DIR_CFG, DIR_DEP, DIR_TAX
from kakapo.config import OS_STR, PY_V_STR
from kakapo.config import THREADS
from kakapo.config_file_parse import config_file_parse
from kakapo.helpers import make_dir
from kakapo.workflow import combine_aa_fasta
from kakapo.workflow import descending_tax_ids
from kakapo.workflow import dnld_pfam_uniprot_seqs
from kakapo.workflow import dnld_prot_seqs
from kakapo.workflow import dnld_sra_fastq_files
from kakapo.workflow import dnld_sra_info
from kakapo.workflow import filter_queries
from kakapo.workflow import makeblastdb_fq
from kakapo.workflow import min_accept_read_len
from kakapo.workflow import pfam_uniprot_accessions
from kakapo.workflow import prepare_output_directories
from kakapo.workflow import run_tblastn_on_reads
from kakapo.workflow import run_trimmomatic
from kakapo.workflow import run_vsearch_on_reads
from kakapo.workflow import trimmed_fq_to_fa
from kakapo.workflow import user_aa_fasta
from kakapo.workflow import user_fastq_files
from kakapo.workflow import user_protein_accessions


# Command line arguments -----------------------------------------------------
CLEAN_CONFIG_DIR = False
CONFIG_FILE_PATH = 'tests/data/kakapo.ini'


def main():
    """Run the script."""
    print('\n\nOperating system: {os}'.format(os=OS_STR))
    print('Python version: {pv}\n'.format(pv=PY_V_STR))

    # Remove configuration directory, if requested ---------------------------
    if CLEAN_CONFIG_DIR and ope(DIR_CFG):
        print('Removing configuration directory:\n\t\t' + DIR_CFG + '\n')
        rmtree(DIR_CFG)
        exit(0)

    # Create config directory with all the subdirectories --------------------
    if ope(DIR_CFG):
        print('Found configuration directory:\n\t\t' + DIR_CFG + '\n')
    else:
        print('Creating configuration directory:\n\t\t' + DIR_CFG + '\n')
        make_dir(DIR_CFG)

    # Check for dependencies -------------------------------------------------
    print('Checking for dependencies:\n')
    make_dir(DIR_DEP)
    seqtk = deps.dep_check_seqtk()  # noqa
    trimmomatic, adapters = deps.dep_check_trimmomatic()  # noqa
    fasterq_dump = deps.dep_check_sra_toolkit()  # noqa
    makeblastdb, blastn, tblastn = deps.dep_check_blast()  # noqa
    vsearch = deps.dep_check_vsearch()  # noqa
    spades = deps.dep_check_spades()  # noqa

    # Initialize NCBI taxonomy database --------------------------------------
    tax = taxonomy(DIR_TAX)  # noqa

    # Parse configuration file -----------------------------------------------
    __ = config_file_parse(CONFIG_FILE_PATH, tax)

    prj_name = __['project_name']
    email = __['email']  # noqa
    dir_out = __['output_directory']
    sras = __['sras']
    fq_pe = __['fq_pe']
    fq_se = __['fq_se']
    assmbl = __['assmbl']  # noqa
    min_query_length = __['min_query_length']
    max_query_length = __['max_query_length']
    user_queries = __['user_queries']
    blast_1_culling_limit = __['blast_1_culling_limit']
    blast_1_evalue = __['blast_1_evalue']
    blast_1_max_target_seqs = __['blast_1_max_target_seqs']
    blast_1_qcov_hsp_perc = __['blast_1_qcov_hsp_perc']
    blast_2_culling_limit = __['blast_2_culling_limit']  # noqa
    blast_2_evalue = __['blast_2_evalue']  # noqa
    blast_2_max_target_seqs = __['blast_2_max_target_seqs']  # noqa
    blast_2_qcov_hsp_perc = __['blast_2_qcov_hsp_perc']  # noqa
    tax_group = __['tax_group']
    tax_group_name = __['tax_group_name']
    tax_ids_user = __['tax_ids']
    pfam_acc = __['pfam_acc']
    prot_acc_user = __['prot_acc']

    # Genetic code information and translation tables ------------------------

    print('Loading genetic code information and translation tables for ' +
          tax_group_name + '\n')

    gc = tax.genetic_code_for_taxid(tax_group)
    # gc_mito = tax.mito_genetic_code_for_taxid(tax_group)
    # gc_plastid = tax.plastid_genetic_code()
    # gc_tt = tax.trans_table_for_genetic_code_id(gc)
    # gc_mito_tt = tax.trans_table_for_genetic_code_id(gc_mito)
    # gc_plastid_tt = tax.trans_table_for_genetic_code_id(gc_plastid)

    # Create output directory with all the subdirectories --------------------
    if dir_out is not None:
        if ope(dir_out):
            print('Found output directory:\n\t\t' + dir_out + '\n')
        else:
            print('Creating output directory:\n\t\t' + dir_out + '\n')
            make_dir(dir_out)

    __ = prepare_output_directories(dir_out, prj_name)

    dir_temp = __['dir_temp']
    # dir_cache = __['dir_cache']
    dir_cache_pfam_acc = __['dir_cache_pfam_acc']
    dir_cache_fq_minlen = __['dir_cache_fq_minlen']
    dir_cache_prj = __['dir_cache_prj']
    # dir_prj = __['dir_prj']
    dir_prj_queries = __['dir_prj_queries']
    dir_fq_data = __['dir_fq_data']
    dir_fq_trim_data = __['dir_fq_trim_data']
    dir_fa_trim_data = __['dir_fa_trim_data']
    dir_blast_fa_trim = __['dir_blast_fa_trim']
    dir_blast_results_fa_trim = __['dir_blast_results_fa_trim']
    dir_vsearch_results_fa_trim = __['dir_vsearch_results_fa_trim']

    # Housekeeping done. Start the analyses. ---------------------------------

    # Resolve descending taxonomy nodes --------------------------------------
    tax_ids = descending_tax_ids(tax_ids_user, tax)

    # Pfam uniprot accessions ------------------------------------------------
    pfam_uniprot_acc = pfam_uniprot_accessions(pfam_acc, tax_ids,
                                               dir_cache_pfam_acc)

    # Download Pfam uniprot sequences if needed ------------------------------
    aa_uniprot_file = opj(dir_prj_queries, 'aa_uniprot.fasta')
    dnld_pfam_uniprot_seqs(pfam_uniprot_acc, aa_uniprot_file, dir_cache_prj)

    # User provided protein accessions ---------------------------------------
    prot_acc_user = user_protein_accessions(prot_acc_user)

    # Download from NCBI if needed -------------------------------------------
    aa_prot_ncbi_file = opj(dir_prj_queries, 'aa_prot_ncbi.fasta')
    dnld_prot_seqs(prot_acc_user, aa_prot_ncbi_file, dir_cache_prj)

    # User provided protein sequences ----------------------------------------
    aa_prot_user_file = opj(dir_prj_queries, 'aa_prot_user.fasta')
    user_aa_fasta(user_queries, aa_prot_user_file)

    # Combine all AA queries -------------------------------------------------
    aa_queries_file = opj(dir_prj_queries, 'aa_all.fasta')
    combine_aa_fasta([aa_uniprot_file,
                      aa_prot_ncbi_file,
                      aa_prot_user_file], aa_queries_file)

    # Filter AA queries ------------------------------------------------------
    filter_queries(aa_queries_file, min_query_length, max_query_length)

    # Download SRA run metadata if needed ------------------------------------
    sra_runs_info = dnld_sra_info(sras, dir_cache_prj)

    # Download SRA run FASTQ files if needed ---------------------------------
    se_fastq_files_sra, pe_fastq_files_sra = dnld_sra_fastq_files(
        sras, sra_runs_info, dir_fq_data, fasterq_dump, THREADS, dir_temp)

    # User provided FASTQ files ----------------------------------------------
    se_fastq_files_usr, pe_fastq_files_usr = user_fastq_files(fq_se, fq_pe)

    # Collate FASTQ file info ------------------------------------------------
    se_fastq_files = se_fastq_files_sra.copy()
    se_fastq_files.update(se_fastq_files_usr)
    pe_fastq_files = pe_fastq_files_sra.copy()
    pe_fastq_files.update(pe_fastq_files_usr)

    # Minimum acceptable read length -----------------------------------------
    min_accept_read_len(se_fastq_files, pe_fastq_files, dir_temp,
                        dir_cache_fq_minlen, vsearch)

    # File name patterns -----------------------------------------------------

    pe_trim_pair_1_sfx = '_paired_1'
    pe_trim_pair_2_sfx = '_paired_2'
    pe_trim_unpr_1_sfx = '_unpaired_1'
    pe_trim_unpr_2_sfx = '_unpaired_2'

    pe_trim_suffixes = [pe_trim_pair_1_sfx, pe_trim_pair_2_sfx,
                        pe_trim_unpr_1_sfx, pe_trim_unpr_2_sfx]

    pe_file_pattern = opj('@D@', '@N@')

    pe_trim_fq_file_patterns = list(
        zip([pe_file_pattern] * 4, pe_trim_suffixes, ['.fastq'] * 4))
    pe_trim_fq_file_patterns = [''.join(x) for x in pe_trim_fq_file_patterns]

    pe_trim_fa_file_patterns = [x.replace('.fastq', '.fasta') for x in
                                pe_trim_fq_file_patterns]

    pe_blast_db_file_patterns = list(
        zip([pe_file_pattern] * 4, pe_trim_suffixes))
    pe_blast_db_file_patterns = [''.join(x) for x in pe_blast_db_file_patterns]

    pe_blast_results_file_patterns = [x.replace('.fastq', '.txt') for x in
                                      pe_trim_fq_file_patterns]

    pe_vsearch_results_file_patterns = pe_blast_results_file_patterns

    # Run Trimmomatic --------------------------------------------------------
    run_trimmomatic(se_fastq_files, pe_fastq_files, dir_fq_trim_data,
                    trimmomatic, adapters, pe_trim_fq_file_patterns, THREADS)

    # Convert trimmed FASTQ files to FASTA -----------------------------------
    trimmed_fq_to_fa(se_fastq_files, pe_fastq_files, dir_fa_trim_data, seqtk,
                     pe_trim_fa_file_patterns)

    # Run makeblastdb --------------------------------------------------------
    makeblastdb_fq(se_fastq_files, pe_fastq_files, dir_blast_fa_trim,
                   makeblastdb, pe_blast_db_file_patterns)

    # Run tblastn ------------------------------------------------------------
    run_tblastn_on_reads(se_fastq_files, pe_fastq_files, aa_queries_file,
                         tblastn, blast_1_evalue, blast_1_max_target_seqs,
                         blast_1_culling_limit, blast_1_qcov_hsp_perc,
                         dir_blast_results_fa_trim,
                         pe_blast_results_file_patterns, THREADS, gc, seqtk,
                         vsearch)

    # Run vsearch ------------------------------------------------------------
    run_vsearch_on_reads(se_fastq_files, pe_fastq_files, vsearch,
                         dir_vsearch_results_fa_trim,
                         pe_vsearch_results_file_patterns, seqtk)

    # print('SE:')
    # for k in se_fastq_files:
    #     print(se_fastq_files[k])

    # print('PE:')
    # for k in pe_fastq_files:
    #     print(pe_fastq_files[k])


##############################################################################


if __name__ == '__main__':
    main()
