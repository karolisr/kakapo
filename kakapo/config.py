# -*- coding: utf-8 -*-

"""
Provides "global" variables.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import with_statement

import os

from kakapo import __script_name__, __version__
from kakapo.helpers import sys_ram
from kakapo.os_diffs import check_os
from kakapo.py_v_diffs import python_version
from multiprocessing import cpu_count

DEBUG_MODE = True
DEBUG_PROCESSES = False

# System information ---------------------------------------------------------
OS_ID, OS_STR, DIST_ID = check_os()
_, PY_V_STR = python_version()
THREADS = cpu_count()
RAM = sys_ram()

# Other ----------------------------------------------------------------------
PICKLE_PROTOCOL = 2

CONSRED = '\033[0;91m'
CONGREE = '\033[0;92m'
CONYELL = '\033[0;93m'
CONBLUE = '\033[0;94m'
CONSDFL = '\033[0m'

MT_PT_KRKN_DB = 'mitochondrion_and_plastid'

# Filesystem paths for configuration directory -------------------------------
DIR_USR = os.path.expanduser('~')
DIR_CFG = os.path.join(DIR_USR, '.config', __script_name__)
DIR_DEP = os.path.join(DIR_CFG, 'dependencies')
DIR_TAX = os.path.join(DIR_CFG, 'ncbi-taxonomy')
DIR_KRK = os.path.join(DIR_CFG, 'kraken2_dbs')

# Script Info ----------------------------------------------------------------
SCRIPT_INFO = ('\n' +
               '{s} version: {v}\n'.format(s=__script_name__.title(),
                                           v=__version__) +
               'Python version: {pv}\n'.format(pv=PY_V_STR) +
               'Operating system: {os}\n'.format(os=OS_STR) +
               'System info: {cpus} CPUs, {ram} GB RAM\n'.format(
                   cpus=THREADS, ram='{0:.2f}'.format(RAM)))
