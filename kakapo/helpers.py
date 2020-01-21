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

from datetime import datetime
from ftplib import FTP
from operator import itemgetter
from os.path import splitext
from urllib.parse import urlparse
import fileinput
import gzip
import hashlib
import os

from kakapo.http_k import download_file as download_file_http
from kakapo.py_v_diffs import zip_longest

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


def list_of_dirs(path, linfo=print):  # noqa
    try:
        ld = [x for x in os.listdir(path) if os.path.isdir(
            os.path.join(path, x))]
    except FileNotFoundError:
        linfo('Directory "{}" does not exist.'.format(path))
    return ld


def list_of_files(path, linfo=print):  # noqa
    try:
        lf = [x for x in os.listdir(path) if os.path.isfile(
            os.path.join(path, x))]
    except FileNotFoundError:
        linfo('Directory "{}" does not exist.'.format(path))
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


def sys_ram(os_id):  # noqa
    ram_b = 0
    try:
        page_size = os.sysconf('SC_PAGE_SIZE')
        page_count = os.sysconf('SC_PHYS_PAGES')
        if page_size < 0 or page_count < 0:
            raise SystemError
        ram_b = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')

    except ValueError:
        if os_id == 'mac':
            ram_b = int(float(os.popen("sysctl hw.memsize").readlines()[0].split()[1]))
        elif os_id == 'linux':
            ram_b = int(float(os.popen("free").readlines()[1].split()[1]) * 1024)
        else:
            print('The OS:' + os_id + ' was not recognized.')
            raise NotImplementedError

    ram_g = ram_b / (1024 ** 3)
    return ram_g


def time_stamp():  # noqa
    return datetime.now().strftime(format='%Y%m%d%H%M%S')


def overlap(a, b):  # noqa
    a = sorted(a)
    b = sorted(b)
    ab = tuple(sorted((a, b), key=itemgetter(0)))
    if ab[0][1] <= ab[1][0]:
        overlap = 0
    else:
        overlap = ab[0][1] - ab[1][0]
        if ab[1][1] < ab[0][1]:
            overlap -= ab[0][1] - ab[1][1]
    return overlap


def download_file(url, local_path, protocol='http'):  # noqa
    assert protocol in ('http', 'ftp')

    if protocol == 'http':
        download_file_http(url, local_path)
    elif protocol == 'ftp':
        url_parsed = urlparse(url)
        netloc = url_parsed.netloc
        path = url_parsed.path

        with FTP(netloc) as ftp:
            ftp.login()
            with open(local_path, 'wb') as f:
                ftp.retrbinary(cmd='RETR ' + path,
                               callback=lambda x: f.write(x),
                               blocksize=8192 * 32)


def splitext_gz(path):
    """
    Split extension for files that may have a double extension: x.y.(gz|gzip)

    Return ('x', '.y', '.(gz|gzip)')
    """
    ext = splitext(path)
    ext_gz = None
    ext_fq = None

    if ext[1] in ('.gz', '.gzip'):
        ext_gz = ext[1]
        ext = splitext(ext[0])

    ext_fq = ext[1]
    base = ext[0]

    return base, ext_fq, ext_gz


def _gzip_open(filename, mode='rt', compresslevel=5, encoding=None,
               errors=None, newline=None):
    return gzip.open(filename, mode, compresslevel, encoding, errors, newline)


def plain_or_gzip(in_file):  # noqa
    read_mode = 'r'
    write_mode = 'w'
    append_mode = 'a'
    fqopen = open
    ext = ''

    ext_info = splitext_gz(in_file)
    if ext_info[2] is not None:
        read_mode = 'rt'
        write_mode = 'wt'
        append_mode = 'at'
        fqopen = _gzip_open
        ext = ext_info[2]

    return read_mode, write_mode, append_mode, fqopen, ext


def grouper(iterable, n, fillvalue=None):
    """
    Break the data into fixed-length chunks or blocks.

    Allows parsing of a (FASTQ) file n (4) lines at a time.

    Edited from "Transcriptome workshop at Botany 2018"
    Ya Yang and Stephen A. Smith, which was modified from the code in:
    FilterUncorrectabledPEfastq.py file by Adam H. Freedman:
    https://github.com/harvardinformatics/TranscriptomeAssemblyTools
    """
    args = [iter(iterable)] * n
    return zip_longest(fillvalue=fillvalue, *args)


def split_mixed_fq(in_file, out_file_1, out_file_2):  # noqa

    r_mode, w_mode, a_mode, fqopen, ext = plain_or_gzip(in_file)

    with fqopen(in_file, r_mode) as in_f, \
            fqopen(out_file_1, w_mode) as out_f_1, \
            fqopen(out_file_2, w_mode) as out_f_2:
        entries = grouper(in_f, 4)
        for entry in entries:
            head, seq, plhld, qual = [i.strip() for i in entry]
            entry_str = '\n'.join([head, seq, plhld, qual])
            if ' 1:N:' in head:
                out_f_1.write('{}\n'.format(entry_str))
            elif ' 2:N:' in head:
                out_f_2.write('{}\n'.format(entry_str))


def split_seq_defn_for_printing(defn):  # noqa
    defn_split = defn.split(' ')
    defn_a = defn_split[0]
    defn_b = ' '.join(defn_split[1:]).split('|')[0]
    defn_b = defn_b.replace('_', ' ').replace('-', ' ')
    defn_b = defn_b[0].upper() + defn_b[1:]
    defn_b = defn_b.split('; RevComp')[0]
    return defn_a, defn_b
