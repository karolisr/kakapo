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
from shutil import rmtree
from sys import exit

from ncbi_taxonomy_local import taxonomy

from kakapo import dependencies as deps
from kakapo.config import DIR_CFG, DIR_DEP, DIR_TAX
from kakapo.config import OS_STR, PY_V_STR
from kakapo.helpers import make_dir
from kakapo.config_file_parse import config_file_parse

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

    project_name = __['project_name']  # noqa
    email = __['email']  # noqa
    output_directory = __['output_directory']  # noqa
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
    tax_ids = __['tax_ids']  # noqa
    pfam_acc = __['pfam_acc']  # noqa
    prot_acc = __['prot_acc']  # noqa

    # Housekeeping done. Start the analyses. ---------------------------------


##############################################################################


if __name__ == '__main__':
    main()
