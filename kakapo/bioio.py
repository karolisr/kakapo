# -*- coding: utf-8 -*-
"""Read and write biological sequence and alignment files."""


from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import with_statement

import re
from collections import OrderedDict

from kakapo.py_v_diffs import handle_types


def write_fasta_file(records, file_path_or_handle):  # noqa
    handle = False

    for h in handle_types:
        if issubclass(type(file_path_or_handle), h):
            handle = True
            break

    if not handle:
        file_path_or_handle = open(file_path_or_handle, 'w')

    for rec in records:

        description = None
        seq = None

        if hasattr(rec, 'description'):
            description = rec.description
        elif hasattr(rec, 'name'):
            description = rec.name
        elif isinstance(rec, dict):
            if 'description' in rec:
                description = rec['description']
            elif 'definition' in rec:
                description = rec['definition']
            elif 'name' in rec:
                description = rec['name']
        else:
            raise Exception('No name, definition, or description attributes'
                            'or keys in record.')

        if hasattr(rec, 'accession'):
            description = rec.accession + '|' + description
        elif isinstance(rec, dict):
            if 'accession' in rec:
                description = rec['accession'] + '|' + description

        if hasattr(rec, 'seq'):
            seq = rec.seq.seq
        elif isinstance(rec, dict):
            if 'seq' in rec:
                seq = rec['seq']
                if hasattr(seq, 'seq'):
                    seq = seq.seq
        else:
            raise Exception('No "seq" attribute or key in record.')

        fasta_entry = '>' + description + '\n' + seq
        file_path_or_handle.write(fasta_entry + '\n')

    if not handle:
        file_path_or_handle.close()


def read_fasta_file(file_path_or_handle):  # noqa
    handle = False

    for h in handle_types:
        if issubclass(type(file_path_or_handle), h):
            handle = True
            break

    if not handle:
        file_path_or_handle = open(file_path_or_handle, 'r')

    seq_names = list()
    seqs = list()
    seq_name = None
    seq = ''
    for ln in file_path_or_handle:
        ln = ln.strip('\n')
        if ln.startswith('>'):
            if seq_name is not None:
                seqs.append(seq)
                seq = ''
            seq_name = ln.strip('>')
            seq_names.append(seq_name)
        else:
            seq = seq + ln

    seqs.append(seq)

    if not handle:
        file_path_or_handle.close()

    return_list = list()
    seq_list = list(zip(seq_names, seqs))
    for s in seq_list:
        rec = {'description': s[0], 'seq': s[1].upper()}
        return_list.append(rec)

    return return_list


def parse_fasta_text(text):  # noqa
    desc_lines = re.findall('\>.*', text)
    lines = text.split('\n')
    data = OrderedDict()
    for l in lines:
        if l in desc_lines:
            desc = l.strip('>')
            data[desc] = ''
        else:
            data[desc] = data[desc] + l.upper()

    return data


def read_fasta_file_dict(in_file):  # noqa
    with open(in_file, 'r') as f:
        fasta_text = f.read()
    parsed_fasta = parse_fasta_text(fasta_text)
    return parsed_fasta


def standardize_fasta_text(text):  # noqa
    parsed_fasta = parse_fasta_text(text)
    t = ''
    for k in parsed_fasta:
        desc = k.replace(' ', '_')
        desc = '>' + desc
        seq = parsed_fasta[k]
        t = t + desc + '\n' + seq + '\n'
    return t


def trim_desc_to_first_space_in_fasta_text(text):  # noqa
    parsed_fasta = parse_fasta_text(text)
    t = ''
    for k in parsed_fasta:
        desc = k.split(' ')[0]
        desc = '>' + desc
        seq = parsed_fasta[k]
        t = t + desc + '\n' + seq + '\n'
    return t


def filter_fasta_text_by_length(fasta_text, min_len, max_len):  # noqa
    parsed_fasta = parse_fasta_text(fasta_text)

    filtered_text = ''

    for k in parsed_fasta:
        desc = '>' + k
        seq = parsed_fasta[k]
        if len(seq) < min_len:
            continue
        if len(seq) > max_len:
            continue

        filtered_text = filtered_text + desc + '\n' + seq + '\n'

    return filtered_text
