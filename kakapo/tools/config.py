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
MACHINE_TYPE = _['machine_type']
OS_ID = _['os_id']
OS_STR = _['os_str']
RELEASE_ID = _['release_id']
RELEASE_NAME = _['release_name']
DIST_ID = _['dist_id']
DEBIAN_DISTS = _['debian_dists']
REDHAT_DISTS = _['redhat_dists']
SUPPORTED_DISTS = _['supported_dists']
OS_INFO = OS_STR

if OS_ID == 'mac':
    OS_INFO += ' ' + RELEASE_ID
    if RELEASE_NAME != '':
        OS_INFO += ' (' + RELEASE_NAME + ')'

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

# Filesystem paths for kakapo data directories -------------------------------
DIR_USR = os.path.expanduser('~')
DIR_DAT = os.path.join(DIR_USR, '.local', 'share', __script_name__)  # XDG_DATA_HOME
DIR_DEP = os.path.join(DIR_DAT, 'dependencies')
DIR_TAX = os.path.join(DIR_DAT, 'ncbi-taxonomy')
DIR_KRK = os.path.join(DIR_DAT, 'kraken2_dbs')

# Script Info ----------------------------------------------------------------
SCRIPT_INFO = ('\n' +
               '{s} version: {v}\n'.format(s=__script_name__.title(),
                                           v=__version__) +
               'Python version: {pv}\n'.format(pv=PY_V_STR) +
               'Operating system: {os}\n'.format(os=OS_INFO) +
               'System info: {cpus} CPUs, {ram} GB RAM ({mt})\n'.format(
                   cpus=THREADS, ram='{0:.2f}'.format(RAM), mt=MACHINE_TYPE))
