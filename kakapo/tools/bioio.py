"""Read and write biological sequence and alignment files."""

import io
import re
from collections import OrderedDict
from functools import reduce
from io import StringIO
from math import ceil
from operator import add
from multipledispatch import dispatch
from kakapo.tools.seq import Seq, SeqRecord, SeqRecordCDS
from kakapo.tools.seq import SEQ_TYPES

from typing import List

BIOIO_NS = dict()
HANDLE_TYPES = (io.IOBase, StringIO)


def read_fasta(f, seq_type, upper=True, def_to_first_space=False,
               parse_def=False, simple_return=False) -> List[SeqRecord] | List[str]:
    """Read a FASTA file."""
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
            if seq_lines is not None:
                seq_lines.append(line)

    for rec in records:
        rec[1] = process_seq(''.join(rec[1]))

    return_object = list()
    for rec in records:
        if simple_return is True:
            seq_record = (rec[0], rec[1])
        else:
            seq_record = SeqRecord(rec[0], Seq(rec[1], seq_type))
            if parse_def is True:
                defn = seq_record.definition
                re_loc_pattern = '^(.*)(\\:)(c*\\d+\\-\\d+,*\\S*)(\\s)(.*$)'
                defn_re_loc = re.findall(re_loc_pattern, defn)
                if len(defn_re_loc) > 0 and len(defn_re_loc[0]) == 5:
                    defn = defn_re_loc[0][0] + ' ' + defn_re_loc[0][-1]
                defn_split = defn.split(sep='|', maxsplit=1)
                if len(defn_split) != 2:
                    defn_split = defn.split(sep=' ', maxsplit=1)
                if len(defn_split) == 2:
                    acc_ver = defn_split[0].strip()
                    acc_ver_split = acc_ver.split('.', maxsplit=1)

                    if len(acc_ver_split) == 1:
                        seq_record.definition = defn_split[1]
                        seq_record.accession = acc_ver_split[0]

                    if len(acc_ver_split) == 2:
                        seq_record.definition = defn_split[1]
                        seq_record.accession = acc_ver_split[0]
                        seq_record.version = acc_ver_split[1]

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


def seq_records_to_dict(records, prepend_acc=False, prepend_org=False):
    d = dict()
    for rec in records:
        if type(rec) in (SeqRecord, SeqRecordCDS):
            dfn = rec.definition
            if prepend_org is True:
                org = rec.organism
                if org is not None:
                    dfn = org.replace(' ', '_') + ' ' + dfn
            if prepend_acc is True:
                acc_ver = rec.accession_version
                if acc_ver is not None:
                    dfn = acc_ver + ' ' + dfn
            d[dfn] = str(rec)
    return d


def seq_records_to_fasta(records, max_line_len=None, prepend_acc=False, prepend_org=False):
    d = seq_records_to_dict(records, prepend_acc, prepend_org)
    return dict_to_fasta(d, max_line_len)


def write_fasta(data, f, max_line_len=None, prepend_acc=False, prepend_org=False):
    """Write a FASTA file."""
    handle = False
    if isinstance(f, HANDLE_TYPES):
        handle = True
    if handle is False:
        f = open(f, 'w')

    if type(data) in (dict, OrderedDict):
        text = dict_to_fasta(data, max_line_len)
        f.write(text)

    if type(data) in (list, tuple):
        text = seq_records_to_fasta(data, max_line_len, prepend_acc, prepend_org)
        f.write(text)

    if handle is False:
        f.close()


@dispatch(str, str, namespace=BIOIO_NS)
def _no_spaces(name, sep='_'):
    return sep.join(re.findall(r'([^\s]+)', name))


@dispatch(SeqRecord, namespace=BIOIO_NS)
def _no_spaces(seq_record, sep='_'):
    return _no_spaces(seq_record.definition, sep)


def standardize_fasta(f, seq_type, pfam=False):
    parsed_fasta = read_fasta(f, seq_type)
    names = tuple(map(_no_spaces, parsed_fasta))
    if pfam is True:
        names = tuple(['|'.join(x.split('|')[1:]) for x in names])
    for i, rec in enumerate(parsed_fasta):
        rec.definition = names[i]
    return parsed_fasta


def standardize_fasta_text(text, seq_type, pfam=False):
    return standardize_fasta(StringIO(text), seq_type, pfam)


def trim_desc_to_first_space_in_fasta_text(text, seq_type):
    parsed_fasta = read_fasta(StringIO(text), seq_type, def_to_first_space=True)
    return parsed_fasta


def filter_fasta_by_length(f, seq_type, min_len, max_len) -> List[SeqRecord]:
    parsed_fasta = read_fasta(f, seq_type)
    recs = filter(lambda y: min_len <= y.length <= max_len, parsed_fasta)
    return list(recs)


def filter_fasta_text_by_length(text, seq_type, min_len, max_len):
    return filter_fasta_by_length(StringIO(text), seq_type, min_len, max_len)
