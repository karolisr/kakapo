"""Provides global variables."""

import os

from os import cpu_count

from kakapo import __script_name__, __version__
from kakapo.utils.misc import sys_ram, python_version
from kakapo.utils.os_diffs import check_os

DEBUG_MODE = False
DEBUG_PROCESSES = False

# System information ---------------------------------------------------------
_ = check_os()
OS_ID = _['os_id']
OS_STR = _['os_str']
DIST_ID = _['dist_id']

_, PY_V_STR = python_version()
THREADS = cpu_count()
RAM = sys_ram(OS_ID)

# Other ----------------------------------------------------------------------
PICKLE_PROTOCOL = 2

CONSRED = '\033[0;91m'
CONGREE = '\033[0;92m'
CONYELL = '\033[0;93m'
CONBLUE = '\033[0;94m'
CONSDFL = '\033[0m'

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
