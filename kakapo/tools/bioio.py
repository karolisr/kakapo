"""Read and write biological sequence and alignment files."""

import io
import re
from collections import OrderedDict
from functools import reduce
from io import StringIO
from math import ceil
from operator import add
from kakapo.tools.seq import Seq, SeqRecord
from kakapo.tools.seq import SEQ_TYPES

HANDLE_TYPES = (io.IOBase, StringIO)


def read_fasta(f, seq_type, upper=True, def_to_first_space=False) -> list:
    """Read FASTA file."""

    assert seq_type.upper() in SEQ_TYPES

    if upper is False:
        def process_seq(seq):
            return seq
    else:
        def process_seq(seq):
            return seq.upper()

    if def_to_first_space is False:
        def process_defn(defn):
            return defn[1:]
    else:
        def process_defn(defn):
            return defn[1:].split(sep=None, maxsplit=1)[0]

    handle = False
    if isinstance(f, HANDLE_TYPES):
        handle = True
    if handle is False:
        f = open(f)

    lines = f.read().splitlines()

    if handle is False:
        f.close()

    records = list()
    seq_lines = None
    for line in lines:
        if line.startswith('>'):
            seq_lines = list()
            records.append([process_defn(line), seq_lines])
        else:
            seq_lines.append(line)

    for rec in records:
        rec[1] = process_seq(''.join(rec[1]))

    return_object = list()
    for rec in records:
        seq_record = SeqRecord(rec[0], Seq(rec[1], seq_type))
        return_object.append(seq_record)

    return return_object


def dict_to_fasta(d, max_line_len=None):
    if len(d) == 0:
        return ''

    if max_line_len is not None:
        def split_seq(seq):
            seq = str(seq)
            l = max_line_len
            s = ''
            for i in range(0, ceil(len(seq) / l)):
                s += seq[i * l:(i + 1) * l] + '\n'
            return s.strip()
    else:
        def split_seq(seq):
            seq = str(seq)
            return seq

    fasta = reduce(
        add,
        ['>' + k + '\n' + split_seq(v) + '\n' for (k, v) in d.items()])
    return fasta


def write_fasta(data, f, max_line_len=None):
    """Write FASTA file."""

    handle = False
    if isinstance(f, HANDLE_TYPES):
        handle = True
    if handle is False:
        f = open(f, 'w')

    if type(data) in (dict, OrderedDict):
        text = dict_to_fasta(data, max_line_len)
        f.write(text)

    if handle is False:
        f.close()


def _no_spaces(name, sep='_'):
    return sep.join(re.findall(r'([^\s]+)', name))


def standardize_fasta_text(text):
    parsed_fasta = read_fasta(StringIO(text))
    names = map(_no_spaces, parsed_fasta)
    seqs = parsed_fasta.values()
    fasta_dict = {k: v for (k, v) in zip(names, seqs)}
    return dict_to_fasta(fasta_dict)


def trim_desc_to_first_space_in_fasta_text(text):
    parsed_fasta = read_fasta(StringIO(text))
    names = map(lambda x: x.split(' ')[0], parsed_fasta)
    seqs = parsed_fasta.values()
    fasta_dict = {k: v for (k, v) in zip(names, seqs)}
    return dict_to_fasta(fasta_dict)


def filter_fasta_text_by_length(text, min_len, max_len):
    x = tuple(read_fasta(StringIO(text)).items())
    seqs = filter(lambda y: min_len <= len(y[1]) <= max_len, x)
    fasta_dict = {k: v for (k, v) in seqs}
    return dict_to_fasta(fasta_dict)
