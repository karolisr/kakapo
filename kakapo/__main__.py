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
from kakapo.helpers import make_dir
from kakapo.config_file_parse import config_file_parse
from kakapo.workflow import descending_tax_ids
from kakapo.workflow import prepare_output_directories
from kakapo.workflow import pfam_uniprot_accessions
from kakapo.workflow import user_protein_accessions
from kakapo.workflow import dnld_pfam_uniprot_seqs
from kakapo.workflow import dnld_prot_seqs
from kakapo.workflow import user_aa_fasta
from kakapo.workflow import combine_aa_fasta
from kakapo.workflow import filter_queries
from kakapo.workflow import dnld_sra_info

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

    prj_name = __['project_name']  # noqa
    email = __['email']  # noqa
    dir_out = __['output_directory']  # noqa
    sras = __['sras']  # noqa
    fq_pe = __['fq_pe']  # noqa
    fq_se = __['fq_se']  # noqa
    assmbl = __['assmbl']  # noqa
    min_query_length = __['min_query_length']  # noqa
    max_query_length = __['max_query_length']  # noqa
    user_queries = __['user_queries']  # noqa
    blast_1_culling_limit = __['blast_1_culling_limit']  # noqa
    blast_1_evalue = __['blast_1_evalue']  # noqa
    blast_1_max_target_seqs = __['blast_1_max_target_seqs']  # noqa
    blast_1_qcov_hsp_perc = __['blast_1_qcov_hsp_perc']  # noqa
    blast_2_culling_limit = __['blast_2_culling_limit']  # noqa
    blast_2_evalue = __['blast_2_evalue']  # noqa
    blast_2_max_target_seqs = __['blast_2_max_target_seqs']  # noqa
    blast_2_qcov_hsp_perc = __['blast_2_qcov_hsp_perc']  # noqa
    tax_group = __['tax_group']  # noqa
    tax_group_name = __['tax_group_name']  # noqa
    tax_ids_user = __['tax_ids']  # noqa
    pfam_acc = __['pfam_acc']  # noqa
    prot_acc_user = __['prot_acc']  # noqa

    # Create output directory with all the subdirectories --------------------
    if dir_out is not None:
        if ope(dir_out):
            print('Found output directory:\n\t\t' + dir_out + '\n')
        else:
            print('Creating output directory:\n\t\t' + dir_out + '\n')
            make_dir(dir_out)

    __ = prepare_output_directories(dir_out, prj_name)

    dir_cache = __['dir_cache']  # noqa
    dir_cache_pfam_acc = __['dir_cache_pfam_acc']  # noqa
    dir_cache_prj = __['dir_cache_prj']  # noqa
    dir_prj = __['dir_prj']  # noqa
    dir_prj_queries = __['dir_prj_queries']  # noqa

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

    print(sra_runs_info)


##############################################################################


if __name__ == '__main__':
    main()
