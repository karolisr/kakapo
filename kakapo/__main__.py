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


##############################################################################


# Command line arguments -----------------------------------------------------
CLEAN_CONFIG_DIR = False
CONFIG_FILE_PATH = None


##############################################################################


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


##############################################################################


if __name__ == '__main__':
    main()


##############################################################################
