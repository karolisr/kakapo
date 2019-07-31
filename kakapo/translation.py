# -*- coding: utf-8 -*-

"""
Translation table manipulation.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import with_statement


from collections import defaultdict

from kakapo.iupac import IUPAC_DNA_DICT


def invert_tt(tt):  # noqa

    if 'trans_table' not in tt.keys():
        msg = 'No trans_table key.'
        raise Exception(msg)

    trans_table = tt['trans_table']
    start_codons = tt['start_codons']
    stop_codons = tt['stop_codons']

    trans_table_inv = defaultdict(list)

    {trans_table_inv[v].append(k) for k, v in trans_table.items()}
    trans_table_inv = dict(trans_table_inv)

    tt_inv = {'trans_table_inv': trans_table_inv,
              'start_codons': start_codons,
              'stop_codons': stop_codons}

    return tt_inv


def ambiguous_codons(codons):  # noqa

    def _ac(codons):

        base_1 = list()
        base_2 = list()
        base_3 = list()

        for c in codons:
            base_1.append(c[0])
            base_2.append(c[1])
            base_3.append(c[2])

        bases_123 = dict()

        for b1 in base_1:
            if b1 not in bases_123:
                bases_123[b1] = dict()
            for c in codons:
                if b1 == c[0]:
                    if c[1] not in bases_123[b1]:
                        bases_123[b1][c[1]] = set()
                    for b2 in base_2:
                        if b2 == c[1]:
                            bases_123[b1][b2].add(c[2])

        amb_codons = []

        for b1 in bases_123:
            for b2 in bases_123[b1]:
                b3 = bases_123[b1][b2]
                b3 = ''.join(sorted(list(b3)))
                b3 = IUPAC_DNA_DICT[b3]

                amb_codons.append(b1 + b2 + b3)

        return amb_codons

    amb_codons = _ac(codons)

    # Repeat with all codons reversed
    codons_rev = list()
    for c in codons:
        c_rev = c[::-1]
        codons_rev.append(c_rev)

    _amb_codons_rev = _ac(codons_rev)

    # Reverse the codons back
    amb_codons_rev = list()
    for c in _amb_codons_rev:
        c_rev = c[::-1]
        amb_codons_rev.append(c_rev)

    amb_codons = sorted(list(set(amb_codons + amb_codons_rev + codons)))

    return amb_codons


def ambiguous_tt(tt):  # noqa
    tt_inv = invert_tt(tt)

    trans_table_inv = tt_inv['trans_table_inv']
    start_codons = tt['start_codons']
    stop_codons = tt['stop_codons']

    trans_table_amb = dict()
    for aa in trans_table_inv:
        codons = trans_table_inv[aa]
        codons_amb = ambiguous_codons(codons)
        for ca in codons_amb:
            trans_table_amb[ca] = aa

    start_codons = list(set(ambiguous_codons(start_codons) + start_codons))
    stop_codons = list(set(ambiguous_codons(stop_codons) + stop_codons))

    tt_amb = {'trans_table': trans_table_amb,
              'start_codons': start_codons,
              'stop_codons': stop_codons}

    return tt_amb
