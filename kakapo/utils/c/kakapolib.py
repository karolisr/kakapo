"""Kakapo C libraries."""

import ctypes

from inspect import currentframe
from inspect import getfile
from os import environ
from os.path import abspath
from os.path import dirname
from os.path import exists as ope
from os.path import join as opj

from kakapo.utils.subp import run
from kakapo.utils.logging import Log

_ = currentframe()
assert _ is not None
FILE_SCRIPT = getfile(_)
DIR_SCRIPT = dirname(abspath(FILE_SCRIPT))

DIR_C_SRC = opj(DIR_SCRIPT, 'src')
DIR_C_LIB = opj(DIR_SCRIPT, 'lib')


def dep_check_kakapolib(force=False, quiet=False):
    kkpl_suffix = environ['KKP_OS_ID'] + '_' + environ['KKP_MACHINE_TYPE']
    so_name = 'kakapolib_' + kkpl_suffix + '.so'
    kkpl = opj(DIR_C_LIB, so_name)
    if not ope(kkpl):
        if quiet is False:
            Log.wrn('Compiling kakapolib.', '')
        run(['make', 'install'], cwd=DIR_C_SRC)
    if ope(kkpl):
        if quiet is False:
            Log.msg('kakapolib is available:', kkpl)
    else:
        Log.err('Compilation of kakapolib failed.', '')
        return None
    return ctypes.CDLL(kkpl)


def _c_str(py_str):
    return ctypes.create_string_buffer(str.encode(py_str))


def fq_avg_read_len(fp):
    kkpl = dep_check_kakapolib(quiet=True)
    if kkpl is None:
        return None
    fp = _c_str(fp)
    return kkpl.fq_avg_read_len(fp)


def fq_avg_read_len_gz(fp):
    kkpl = dep_check_kakapolib(quiet=True)
    if kkpl is None:
        return None
    fp = _c_str(fp)
    return int(kkpl.fq_avg_read_len_gz(fp))
