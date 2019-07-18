# -*- coding: utf-8 -*-

"""
Defines various character sets used in biological sequence representations.

Reference: http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import with_statement

from kakapo.py_v_diffs import maketrans


NT_SHARED_CHARS = set('ACG')
NT_AMBIGUOUS_CHARS = set('BDHKMNRSVWY')
NT_GAPS_STRING = '-.'
NT_GAPS_CHARS = set(NT_GAPS_STRING)
DNA_ONLY_CHARS = set('T')
RNA_ONLY_CHARS = set('U')

AA_CHARS = set('GAVLIPFYCMHKRWSTDENQ.')
AA_AMBIGUOUS_CHARS = set('BZX')
AA_GAPS_CHARS = set('-')

DNA_UNAMBIGUOUS = NT_SHARED_CHARS | DNA_ONLY_CHARS
DNA_UNAMBIGUOUS_GAPS = DNA_UNAMBIGUOUS | NT_GAPS_CHARS
DNA_AMBIGUOUS = DNA_UNAMBIGUOUS | NT_AMBIGUOUS_CHARS
DNA_AMBIGUOUS_GAPS = DNA_AMBIGUOUS | NT_GAPS_CHARS

RNA_UNAMBIGUOUS = NT_SHARED_CHARS | RNA_ONLY_CHARS
RNA_UNAMBIGUOUS_GAPS = RNA_UNAMBIGUOUS | NT_GAPS_CHARS
RNA_AMBIGUOUS = RNA_UNAMBIGUOUS | NT_AMBIGUOUS_CHARS
RNA_AMBIGUOUS_GAPS = RNA_AMBIGUOUS | NT_GAPS_CHARS

NT_UNAMBIGUOUS = DNA_UNAMBIGUOUS | RNA_UNAMBIGUOUS
NT_UNAMBIGUOUS_GAPS = DNA_UNAMBIGUOUS_GAPS | RNA_UNAMBIGUOUS_GAPS
NT_AMBIGUOUS = DNA_AMBIGUOUS | RNA_AMBIGUOUS
NT_AMBIGUOUS_GAPS = DNA_AMBIGUOUS_GAPS | RNA_AMBIGUOUS_GAPS

AA_UNAMBIGUOUS = AA_CHARS
AA_UNAMBIGUOUS_GAPS = AA_UNAMBIGUOUS | AA_GAPS_CHARS
AA_AMBIGUOUS = AA_UNAMBIGUOUS | AA_AMBIGUOUS_CHARS
AA_AMBIGUOUS_GAPS = AA_AMBIGUOUS | AA_GAPS_CHARS

DNA_COMPLEMENT_CHARS_FROM = 'ACGTRYMKWSBDHV'
DNA_COMPLEMENT_CHARS_TO = 'TGCAYRKMWSVHDB'
DNA_COMPLEMENT_TABLE = maketrans(DNA_COMPLEMENT_CHARS_FROM,
                                 DNA_COMPLEMENT_CHARS_TO)

IUPAC_DNA_DICT = {
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'T',
    'AG': 'R',
    'CT': 'Y',
    'AC': 'M',
    'GT': 'K',
    'AT': 'W',
    'CG': 'S',
    'CGT': 'B',
    'AGT': 'D',
    'ACT': 'H',
    'ACG': 'V'
}

IUPAC_DNA_DICT_REVERSE = {
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'T',
    'R': 'AG',
    'Y': 'CT',
    'M': 'AC',
    'K': 'GT',
    'W': 'AT',
    'S': 'CG',
    'B': 'CGT',
    'D': 'AGT',
    'H': 'ACT',
    'V': 'ACG'
}

IUPAC_AMBIGUOUS_DNA_DICT = {
    'AG': 'R',
    'CT': 'Y',
    'AC': 'M',
    'GT': 'K',
    'AT': 'W',
    'CG': 'S',
    'CGT': 'B',
    'AGT': 'D',
    'ACT': 'H',
    'ACG': 'V'
}
