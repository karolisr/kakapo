# -*- coding: utf-8 -*-

"""
Basic input/output operations.
"""

import fileinput
import os
import sys
import hashlib

from kakapo.debug import debug_print
from kakapo.py_v_diffs import urlretrieve
from kakapo.shell import call


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


def download_file(url, local_path):  # noqa

    try:
        debug_print('Trying urlretrieve')
        urlretrieve(url, local_path)
    except Exception:
        try:
            debug_print('Trying curl')
            call(['/usr/bin/curl', '-L', '-o', local_path, url])
        except Exception:
            try:
                debug_print('Trying wget')
                call(['wget', '-O', local_path, url])
            except Exception:
                print("\nDownload operation failed.")
                sys.exit(1)


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
