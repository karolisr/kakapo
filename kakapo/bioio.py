# -*- coding: utf-8 -*-
"""Read and write biological sequence and alignment files."""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import with_statement

import re
from functools import reduce
from operator import add
from collections import OrderedDict

from kakapo.py_v_diffs import HANDLE_TYPES, StringIO


def read_fasta(f, return_type='dict', upper=True):
    """Read FASTA file."""
    assert return_type in ('dict', 'list', 'tuple')

    if upper is False:
        def process_seq(seq):
            return seq
    else:
        def process_seq(seq):
            return seq.upper()

    handle = False
    if isinstance(f, HANDLE_TYPES):
        handle = True
    if handle is False:
        f = open(f, 'r')

    lines = f.read().splitlines()

    if handle is False:
        f.close()

    records = list()
    for line in lines:
        if line.startswith('>'):
            seq_lines = list()
            records.append([line[1:], seq_lines])
        else:
            seq_lines.append(line)

    for rec in records:
        rec[1] = process_seq(''.join(rec[1]))

    if return_type == 'dict':
        return_object = dict()
        for rec in records:
            return_object[rec[0]] = rec[1]
    elif return_type == 'tuple':
        return_object = tuple(records)
    else:
        return_object = records

    return return_object


def dict_to_fasta(d):  # noqa
    if len(d) == 0:
        return ''
    return reduce(add, ['>' + k + '\n' + v + '\n' for (k, v) in d.items()])


def _no_spaces(name, sep='_'):
    return sep.join(re.findall(r'([^\s]+)', name))


def write_fasta(data, f, handle_types=HANDLE_TYPES):
    """Write FASTA"""
    def rec_name(r):
        org = r['organism']
        dfn = r['definition'].replace(org, '').strip(' []')
        x = r['accession'] + '.' + r['version'] + '|' + \
            dfn.title() + '|' + org
        return x

    def rec_seq(r):
        return r['seq'].upper()

    handle = False

    if isinstance(f, handle_types):
        handle = True

    if handle is False:
        f = open(f, 'w')

    if type(data) in (dict, OrderedDict):
        text = dict_to_fasta(data)

    elif type(data) in (list, tuple):
        names = map(_no_spaces, map(rec_name, data))
        seqs = map(rec_seq, data)
        fasta_dict = {k: v for (k, v) in zip(names, seqs)}
        text = dict_to_fasta(fasta_dict)

    f.write(text)

    if handle is False:
        f.close()


def standardize_fasta_text(text):  # noqa
    parsed_fasta = read_fasta(StringIO(text))
    names = map(_no_spaces, parsed_fasta)
    seqs = parsed_fasta.values()
    fasta_dict = {k: v for (k, v) in zip(names, seqs)}
    return dict_to_fasta(fasta_dict)


def trim_desc_to_first_space_in_fasta_text(text):  # noqa
    parsed_fasta = read_fasta(StringIO(text))
    names = map(lambda x: x.split(' ')[0], parsed_fasta)
    seqs = parsed_fasta.values()
    fasta_dict = {k: v for (k, v) in zip(names, seqs)}
    return dict_to_fasta(fasta_dict)


def filter_fasta_text_by_length(text, min_len, max_len):  # noqa
    x = tuple(read_fasta(StringIO(text)).items())
    seqs = filter(lambda y: min_len <= len(y[1]) <= max_len, x)
    fasta_dict = {k: v for (k, v) in seqs}
    return dict_to_fasta(fasta_dict)
