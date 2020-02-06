#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Generate start codon context Python files"""

import autopep8
import inspect
import os
import sys

from csv import reader
from itertools import starmap
from operator import truediv
from pprint import pprint

##############################################################################
SCRIPT_FILE_PATH = inspect.getfile(inspect.currentframe())
SCRIPT_DIR_PATH = os.path.dirname(os.path.abspath(SCRIPT_FILE_PATH))
KAKAPO_DIR_PATH = os.path.sep.join(SCRIPT_DIR_PATH.split(os.path.sep)[0:-2])
sys.path.insert(0, KAKAPO_DIR_PATH)
SAVED_WD = os.getcwd()
##############################################################################

def _read_start_codon_context_csv(file_path):  # noqa
    with open(file_path, 'r') as f:
        rdr = reader(f)
        cntx_prcnt = tuple(zip(*map(lambda x: tuple(map(int, x)), rdr)))
        cntx_with_sums = tuple(zip(map(sum, cntx_prcnt), cntx_prcnt))
        paired = map(lambda x: tuple(zip(x[1], [x[0]] * len(x[1]))),
                     cntx_with_sums)
        cntx_frac = tuple(map(lambda x: tuple(starmap(truediv, x)),
                              paired))
        return cntx_frac


def read_l(file_path):  # noqa
    return tuple(reversed(_read_start_codon_context_csv(file_path)))


def read_r(file_path):  # noqa
    return _read_start_codon_context_csv(file_path)


def generate(dir_path, f_out):  # noqa

    file_names = tuple(sorted(filter(lambda x: x.endswith('.csv'),
                                     os.listdir(dir_path))))

    file_paths = tuple(map(lambda x: os.path.join(dir_path, x), file_names))
    var_names = tuple(map(lambda x: x.split('.csv')[0], file_names))
    round2 = lambda x: tuple(map(lambda y: round(y, 2), x))  # noqa

    context_dict = dict()

    for i in range(len(var_names)):
        file_path = file_paths[i]
        var_name = var_names[i]
        var_name_split = var_name.split('_')
        taxid = var_name_split[0].split('taxid')[1]
        side = var_name_split[1]

        context = None
        if side == 'L':
            context = read_l(file_path)
        elif side == 'R':
            context = read_r(file_path)
        else:
            continue

        context = tuple(map(round2, context))
        context_dict[taxid + '_' + side] = context

    f = open(f_out, 'w')
    print('"""Start Codon Context"""\n\n', end='', file=f)
    print('contexts' + ' = ', end='', file=f)
    pprint(context_dict, f)
    f.close()

    txt = autopep8.fix_file(f_out)
    f = open(f_out, 'w')
    f.write(txt)
    f.close()


if __name__ == '__main__':

    os.chdir(SCRIPT_DIR_PATH)

    d_data = 'start_codon_context_csv'
    f_out = 'start_codon_context.py'

    generate(d_data, f_out)

    os.chdir(SAVED_WD)
