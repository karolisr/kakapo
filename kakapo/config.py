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

from multiprocessing import cpu_count

from kakapo.os_diffs import check_os
from kakapo.py_v_diffs import py_v_str

DEBUG_MODE = True
DEBUG_PROCESSES = False

SCRIPT_NAME = 'kakapo'
OS_ID, OS_STR, DIST_ID = check_os()
PY_V_STR = py_v_str

# System information ---------------------------------------------------------
THREADS = cpu_count()

# Filesystem paths for configuration directory -------------------------------
DIR_USR = os.path.expanduser('~')
DIR_CFG = os.path.join(DIR_USR, '.config', SCRIPT_NAME)
DIR_DEP = os.path.join(DIR_CFG, 'dependencies')
