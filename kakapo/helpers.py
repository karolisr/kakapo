# -*- coding: utf-8 -*-

"""
Basic input/output operations.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import with_statement

import fileinput
import hashlib
import os

from datetime import datetime


def debug_print(msg=''): # noqa
    from kakapo.config import DEBUG_MODE
    if DEBUG_MODE:
        import pprint
        PP = pprint.PrettyPrinter(indent=1, width=110, compact=True)
        PP.pprint(msg)
    else:
        pass


def replace_line_in_file(file_path, line_str, replace_str):  # noqa
    for line in fileinput.input(file_path, inplace=1):
        line_strip = line.strip('\n')
        if line_strip == line_str:
            print(replace_str)
        else:
            print(line_strip)


def make_dir(path):  # noqa
    if not os.path.exists(path):
        os.makedirs(path)
    return path


def list_of_dirs(path):  # noqa
    ld = [x for x in os.listdir(path) if os.path.isdir(os.path.join(path, x))]
    return ld


def list_of_files(path):  # noqa
    lf = [x for x in os.listdir(path) if os.path.isfile(os.path.join(path, x))]
    return lf


def generate_md5_hash_for_file(file_path):  # noqa
    return_value = None
    with open(file_path, 'rb') as f:
        return_value = hashlib.md5(f.read()).hexdigest()
    return return_value


def extract_md5_hash(file_path):  # noqa
    md5_reported = None
    with open(file_path, 'r') as f_md5:
        line = f_md5.readline()
        md5_reported = line.split(' ')[0]
    return md5_reported


def unique_lines_in_file(path):  # noqa
    with open(path, 'r') as f:
        x = f.readlines()
    x = list(set(x))
    x.sort()
    return x


def keep_unique_lines_in_file(path):  # noqa
    x = unique_lines_in_file(path)
    with open(path, 'w') as f:
        f.write(''.join(x))


def combine_text_files(paths, out_path):  # noqa
    ret = ''
    for p in paths:
        with open(p, 'r') as f:
            x = f.read()
            ret = ret + x

    with open(out_path, 'w') as f:
        f.write(ret)


def sys_ram():  # noqa
    ram_b = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')
    ram_g = ram_b / (1024**3)
    return ram_g


def time_stamp():  # noqa
    return datetime.now().strftime(format='%Y%m%d%H%M%S')
