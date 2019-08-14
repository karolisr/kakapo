# -*- coding: utf-8 -*-
"""Translation Tables"""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import with_statement

from collections import OrderedDict
from collections import defaultdict
from difflib import SequenceMatcher

from kakapo.iupac import IUPAC_DNA_DICT
from kakapo.iupac import IUPAC_AMBIGUOUS_SECOND_ORDER_DNA_DICT_REVERSE


GC_ID_NAME_MAP = {

    '01': 'Standard',
    '02': 'Vertebrate Mitochondrial',
    '03': 'Yeast Mitochondrial',
    '04': 'Mold Mitochondrial; Protozoan Mitochondrial; '
          'Coelenterate Mitochondrial; Mycoplasma; Spiroplasma',
    '05': 'Invertebrate Mitochondrial',
    '06': 'Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear',
    '09': 'Echinoderm Mitochondrial; Flatworm Mitochondrial',
    '10': 'Euplotid Nuclear',
    '11': 'Bacterial, Archaeal and Plant Plastid',
    '12': 'Alternative Yeast Nuclear',
    '13': 'Ascidian Mitochondrial',
    '14': 'Alternative Flatworm Mitochondrial',
    '15': 'Blepharisma Macronuclear',
    '16': 'Chlorophycean Mitochondrial',
    '21': 'Trematode Mitochondrial',
    '22': 'Scenedesmus obliquus mitochondrial',
    '23': 'Thraustochytrium mitochondrial code',
    '24': 'Pterobranchia Mitochondrial',
    '25': 'Candidate Division SR1 and Gracilibacteria',
    '26': 'Pachysolen tannophilus Nuclear',
    '27': 'Karyorelict Nuclear',
    '28': 'Condylostoma Nuclear',
    '29': 'Mesodinium Nuclear',
    '30': 'Peritrich Nuclear',
    '31': 'Blastocrithidia Nuclear',
    '32': 'Balanophoraceae Plastid',
    '33': 'Cephalodiscidae Mitochondrial'
}

GC_ID_AA_MAP = {

    '01': 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    '02': 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG',
    '03': 'FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    '04': 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    '05': 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG',
    '06': 'FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    '09': 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
    '10': 'FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    '11': 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    '12': 'FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    '13': 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG',
    '14': 'FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
    '15': 'FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    '16': 'FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    '21': 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
    '22': 'FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    '23': 'FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    '24': 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG',
    '25': 'FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    '26': 'FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    '27': 'FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    '28': 'FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    '29': 'FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    '30': 'FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    '31': 'FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    '32': 'FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    '33': 'FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG'
}

GC_ID_START_STOP_MAP = {

    '01': '---M------**--*----M---------------M----------------------------',
    '02': '----------**--------------------MMMM----------**---M------------',
    '03': '----------**----------------------MM---------------M------------',
    '04': '--MM------**-------M------------MMMM---------------M------------',
    '05': '---M------**--------------------MMMM---------------M------------',
    '06': '--------------*--------------------M----------------------------',
    '09': '----------**-----------------------M---------------M------------',
    '10': '----------**-----------------------M----------------------------',
    '11': '---M------**--*----M------------MMMM---------------M------------',
    '12': '----------**--*----M---------------M----------------------------',
    '13': '---M------**----------------------MM---------------M------------',
    '14': '-----------*-----------------------M----------------------------',
    '15': '----------*---*--------------------M----------------------------',
    '16': '----------*---*--------------------M----------------------------',
    '21': '----------**-----------------------M---------------M------------',
    '22': '------*---*---*--------------------M----------------------------',
    '23': '--*-------**--*-----------------M--M---------------M------------',
    '24': '---M------**-------M---------------M---------------M------------',
    '25': '---M------**-----------------------M---------------M------------',
    '26': '----------**--*----M---------------M----------------------------',
    '27': '--------------*--------------------M----------------------------',
    '28': '----------**--*--------------------M----------------------------',
    '29': '--------------*--------------------M----------------------------',
    '30': '--------------*--------------------M----------------------------',
    '31': '----------**-----------------------M----------------------------',
    '32': '---M------*---*----M------------MMMM---------------M------------',
    '33': '---M-------*-------M---------------M---------------M------------'
}

CODONS = [

    'TTT', 'TTC', 'TTA', 'TTG', 'TCT', 'TCC', 'TCA', 'TCG',
    'TAT', 'TAC', 'TAA', 'TAG', 'TGT', 'TGC', 'TGA', 'TGG',
    'CTT', 'CTC', 'CTA', 'CTG', 'CCT', 'CCC', 'CCA', 'CCG',
    'CAT', 'CAC', 'CAA', 'CAG', 'CGT', 'CGC', 'CGA', 'CGG',
    'ATT', 'ATC', 'ATA', 'ATG', 'ACT', 'ACC', 'ACA', 'ACG',
    'AAT', 'AAC', 'AAA', 'AAG', 'AGT', 'AGC', 'AGA', 'AGG',
    'GTT', 'GTC', 'GTA', 'GTG', 'GCT', 'GCC', 'GCA', 'GCG',
    'GAT', 'GAC', 'GAA', 'GAG', 'GGT', 'GGC', 'GGA', 'GGG'
]


def validate_gc_id(gc_id):  # noqa
    gc_id = int(gc_id)
    gc_id = '{:02d}'.format(gc_id)

    if gc_id in GC_ID_NAME_MAP:
        return gc_id
    else:
        return None


def get_gen_code_name(gc_id):  # noqa
    gc_id = validate_gc_id(gc_id)
    if gc_id is None:
        return None
    gc_name = GC_ID_NAME_MAP[gc_id]
    return gc_name


def get_trans_table(gc_id):  # noqa
    gc_id = validate_gc_id(gc_id)
    if gc_id is None:
        return None

    amino_acids = tuple(GC_ID_AA_MAP[gc_id])
    trans_table = OrderedDict(zip(CODONS, amino_acids))
    trans_table = OrderedDict(sorted(trans_table.items(), key=lambda x: x[0]))

    return trans_table


def get_start_codons(gc_id):  # noqa
    gc_id = validate_gc_id(gc_id)
    if gc_id is None:
        return None

    start_stop = tuple(GC_ID_START_STOP_MAP[gc_id])
    start_idx = [i for i, x in enumerate(start_stop) if x == 'M']
    start_codons = tuple(sorted([CODONS[i] for i in start_idx]))

    return start_codons


def get_stop_codons(gc_id):  # noqa
    gc_id = validate_gc_id(gc_id)
    if gc_id is None:
        return None

    start_stop = tuple(GC_ID_START_STOP_MAP[gc_id])
    stop_idx = [i for i, x in enumerate(start_stop) if x == '*']
    stop_codons = tuple(sorted([CODONS[i] for i in stop_idx]))

    return stop_codons


def ambiguous_codons(codons):  # noqa
    IASODDR = IUPAC_AMBIGUOUS_SECOND_ORDER_DNA_DICT_REVERSE

    def _ac(codons):
        match_dict = dict()
        for i in range(1, len(codons)):
            a = codons[i - 1]
            b = codons[i]
            sm = SequenceMatcher(a=a, b=b)
            match = sm.find_longest_match(0, 2, 0, 2)
            if match.a == 0 and match.b == 0 and match.size == 2:
                matching_prefix = a[0:2]
                if matching_prefix not in match_dict:
                    match_dict[matching_prefix] = set()
                match_dict[matching_prefix].add(a[2])
                match_dict[matching_prefix].add(b[2])

        codons_amb = list()
        for prefix in match_dict:
            suffixes = ''.join(sorted(match_dict[prefix]))
            suffix_amb = IUPAC_DNA_DICT[suffixes]
            codons_amb.append(prefix + suffix_amb)
        return codons_amb

    c1 = _ac(codons)
    # repeat with each codon reversed
    codons_reversed = tuple(map(lambda x: x[::-1], codons))
    c2_rev = _ac(codons_reversed)
    # reverse each codon back
    c2 = list(map(lambda x: x[::-1], c2_rev))
    c_grouped = set(c1 + c2 + list(codons))

    c_grouped_temp = set()
    for c in c_grouped:
        b2 = c[1]
        if b2 in IASODDR:
            for b2_alt in IASODDR[b2]:
                c_grouped_temp.add(c[0] + b2_alt + c[2])
    c_grouped = c_grouped | c_grouped_temp

    c_grouped_temp = set()
    for c in c_grouped:
        b1 = c[0]
        if b1 in IASODDR:
            for b1_alt in IASODDR[b1]:
                c_grouped_temp.add(b1_alt + c[1:3])
    c_grouped = c_grouped | c_grouped_temp

    c_grouped_temp = set()
    for c in c_grouped:
        b3 = c[2]
        if b3 in IASODDR:
            for b3_alt in IASODDR[b3]:
                c_grouped_temp.add(c[0:2] + b3_alt)
    c_grouped = c_grouped | c_grouped_temp

    codons_ambiguous = tuple(sorted(c_grouped))
    return codons_ambiguous


def invert_dict(d):  # noqa
    """Will only work if the values in the provided dict are hashable."""
    d_type = type(d)
    d_inv = defaultdict(list)
    {d_inv[v].append(k) for k, v in d.items()}
    for k in d_inv:
        d_inv[k] = tuple(d_inv[k])
    d_inv = d_type(d_inv)
    if d_type is OrderedDict:
        d_inv = OrderedDict(sorted(d_inv.items(), key=lambda x: x[0]))
    return d_inv


def ambiguous_table(trans_table):  # noqa
    d_type = type(trans_table)
    tt_inv = invert_dict(trans_table)
    tt_amb = d_type()
    for k in tt_inv:
        codons = tt_inv[k]
        codons_amb = ambiguous_codons(codons)
        for ca in codons_amb:
            tt_amb[ca] = k
    if d_type is OrderedDict:
        tt_amb = OrderedDict(sorted(tt_amb.items(), key=lambda x: x[0]))
    return tt_amb


class TranslationTable(object):
    """Translation Table"""

    def __init__(self, gc_id):  # noqa

        gc_id = validate_gc_id(gc_id)

        if gc_id is None:
            message = 'Invalid genetic code ID.'
            raise Exception(message)

        self._gc_id = int(gc_id)
        self._gc_name = get_gen_code_name(self._gc_id)
        self._table = get_trans_table(self._gc_id)
        self._table_inv = invert_dict(self._table)
        self._table_ambiguous = ambiguous_table(self._table)
        self._table_ambiguous_inv = invert_dict(self._table_ambiguous)
        self._start_codons = get_start_codons(self._gc_id)
        self._stop_codons = get_stop_codons(self._gc_id)
        self._start_codons_ambiguous = ambiguous_codons(self._start_codons)
        self._stop_codons_ambiguous = ambiguous_codons(self._stop_codons)

    @property
    def gc_id(self):  # noqa
        return self._gc_id

    @property
    def gc_name(self):  # noqa
        return self._gc_name

    @property
    def table(self):  # noqa
        return dict(self._table)

    @property
    def table_inv(self):  # noqa
        return dict(self._table_inv)

    @property
    def table_ambiguous(self):  # noqa
        return dict(self._table_ambiguous)

    @property
    def table_ambiguous_inv(self):  # noqa
        return dict(self._table_ambiguous_inv)

    @property
    def start_codons(self):  # noqa
        return self._start_codons

    @property
    def stop_codons(self):  # noqa
        return self._stop_codons

    @property
    def start_codons_ambiguous(self):  # noqa
        return self._start_codons_ambiguous

    @property
    def stop_codons_ambiguous(self):  # noqa
        return self._stop_codons_ambiguous
