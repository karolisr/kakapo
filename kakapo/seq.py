# -*- coding: utf-8 -*-

"""

Classes that deal with amino acid and nucleotide sequences.

"""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import with_statement

from datetime import date
from random import randint
from re import match

from kakapo.iupac import *
from kakapo.py_v_diffs import basestring
from kakapo.translation import ambiguous_tt

SEQ_TYPE_NT = 'NT'
SEQ_TYPE_DNA = 'DNA'
SEQ_TYPE_RNA = 'RNA'
SEQ_TYPE_AA = 'AA'

SEQ_TYPES = [SEQ_TYPE_NT, SEQ_TYPE_DNA, SEQ_TYPE_RNA, SEQ_TYPE_AA]

MOL_TO_SEQ_TYPE_MAP = {
    'NT': SEQ_TYPE_NT,
    'DNA': SEQ_TYPE_DNA,
    'RNA': SEQ_TYPE_RNA,
    'mRNA': SEQ_TYPE_RNA,
    'rRNA': SEQ_TYPE_RNA,
    'tRNA': SEQ_TYPE_RNA,
    'AA': SEQ_TYPE_AA}


def reverse_complement(seq):  # noqa
    seq = seq.upper()
    comp = seq.translate(DNA_COMPLEMENT_TABLE)
    comp_rev = comp[::-1]
    return(comp_rev)


def resolve_ambiguities(seq):  # noqa
    seq = seq.upper()
    for k in IUPAC_AMBIGUOUS_DNA_DICT:
        for i in range(0, seq.count(IUPAC_AMBIGUOUS_DNA_DICT[k])):
            rand = randint(0, len(k) - 1)
            seq = seq.replace(IUPAC_AMBIGUOUS_DNA_DICT[k], k[rand], 1)
    return(seq)


def translate(seq, trans_table):  # noqa
    trans_table = ambiguous_tt(trans_table)
    seq = seq.upper()
    strtc = trans_table['start_codons']
    tbl = trans_table['trans_table']
    seq_codons = [(seq[i:i + 3]) for i in range(0, len(seq), 3)]

    idx = 0
    len_s = len(seq_codons)
    trans = ''

    if seq_codons[0] in strtc:
        trans = 'M'
        idx = 1

    trans = ''

    for t in seq_codons[idx:len_s]:
        aa = 'X'
        if t in tbl:
            aa = tbl[t]

        trans = trans + aa

    return trans


class Seq(object):
    """
    Seq

    This class can be used to create new amino acid and nucleotide
    sequence objects of types: :class:`NTSeq`,
    :class:`DNASeq`, :class:`RNASeq`, :class:`AASeq`.

    :param seq: Nucleotide or amino acid sequence.
    :type seq: str

    :param seq_type: The type of the sequence. One of the constants
        defined in module :mod:`krpy.seq` should be used:
        ``SEQ_TYPE_NT``, ``SEQ_TYPE_DNA``, ``SEQ_TYPE_RNA``,
        ``SEQ_TYPE_AA``. These constants have string values ``NT``,
        ``DNA``, ``RNA``, ``AA``, respectively.
    :type seq_type: str

    """

    def __new__(self, seq, seq_type):  # noqa

        if seq_type in SEQ_TYPES:

            seq = seq.upper()

            if seq_type is SEQ_TYPE_NT:
                return NTSeq(seq=seq)
            elif seq_type is SEQ_TYPE_DNA:
                return DNASeq(seq=seq)
            elif seq_type is SEQ_TYPE_RNA:
                return RNASeq(seq=seq)
            elif seq_type is SEQ_TYPE_AA:
                return AASeq(seq=seq)

        else:
            # ToDo: Report Exception
            pass

    # __init__ declaration is exactly the same as __new__ so Sphinx
    # docstring parser picks it up.
    def __init__(self, seq, seq_type):  # noqa
        pass


class _Seq(object):
    """

    An abstract class.

    """

    def __init__(self, seq):
        self._seq = seq

    @property
    def seq(self):
        return self._seq

    @property
    def length(self):
        return len(self.seq)


class NTSeq(_Seq):
    """

    NTSeq

    This class represents any nucleotide sequence and does not
    distinguish between RNA and DNA. Both thymine and uracil can occur
    in the instances of this class. :class:`DNASeq` and :class:`RNASeq`
    classes inherit from this class.

    :param seq: Nucleotide sequence.
    :type seq: str

    """

    def __init__(self, seq):  # noqa
        if set(seq) <= NT_AMBIGUOUS:
            super(NTSeq, self).__init__(seq)
        else:
            message = ('NT sequence should not contain these characters: {s}.')
            message = message.format(s=', '.join(
                str(s) for s in set(seq) - NT_AMBIGUOUS))
            raise Exception(message)


class DNASeq(NTSeq):
    """

    Represents DNA sequence.

    :param seq: DNA sequence.
    :type seq: str

    """

    def __init__(self, seq):  # noqa
        if set(seq) <= DNA_AMBIGUOUS:
            super(DNASeq, self).__init__(seq)
        elif set(seq) - DNA_AMBIGUOUS == RNA_ONLY_CHARS:
            seq = seq.replace('U', 'T')
            super(DNASeq, self).__init__(seq)
        else:
            message = (
                'DNA sequence should not contain these characters: {s}.')
            message = message.format(s=', '.join(
                str(s) for s in set(seq) - DNA_AMBIGUOUS))
            raise Exception(message)


class RNASeq(NTSeq):
    """

    Represents RNA sequence.

    :param seq: RNA sequence.
    :type seq: str

    """

    def __init__(self, seq):  # noqa
        if set(seq) <= RNA_AMBIGUOUS:
            super(RNASeq, self).__init__(seq)
        elif set(seq) - RNA_AMBIGUOUS == DNA_ONLY_CHARS:
            seq = seq.replace('T', 'U')
            super(RNASeq, self).__init__(seq)
        else:
            message = (
                'RNA sequence should not contain these characters: {s}.')
            message = message.format(s=', '.join(
                str(s) for s in set(seq) - RNA_AMBIGUOUS))
            raise Exception(message)


class AASeq(_Seq):
    """

    Represents amino acid sequence.

    :param seq: Amino acid sequence.
    :type seq: str

    """

    def __init__(self, seq):  # noqa
        if set(seq) <= AA_AMBIGUOUS:
            super(AASeq, self).__init__(seq)
        else:
            message = ('AA sequence should not contain these characters: {s}.')
            message = message.format(s=', '.join(
                str(s) for s in set(seq) - AA_AMBIGUOUS))
            raise Exception(message)


class SeqRecord(object):
    """

    SeqRecord

    This class stores a sequence and its associated meta-information. It
    is designed to accomodate GenBank records. For reference see:
    http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html

    """

    def __init__(self,  # noqa
                 seq,
                 mol_type,
                 accession=None,
                 version=None,
                 # gi=None,
                 description=None,
                 strandedness=None,
                 topology=None,
                 division=None,
                 date_create=None,
                 date_update=None,
                 taxid=None,
                 organism=None,
                 # taxonomy=None,
                 features=None):

        # immutable setters
        self._set_mol_type(mol_type)
        self._set_seq_type(self.mol_type)
        self._set_seq(seq)

        # mutable setters
        self.accession = accession
        self.version = version
        self.description = description

        self.strandedness = strandedness
        self.topology = topology
        self.division = division

        self.date_create = date_create
        self.date_update = date_update

        self.taxid = taxid
        self.organism = organism
        # self.taxonomy = taxonomy

        self.features = features

    # seq_type and mol_type are not the same. mol_type is more specific.
    # For example: mol_type mRNA, seq_type RNA

    # mol_type
    @property
    def mol_type(self):  # noqa
        return self._mol_type

    def _set_mol_type(self, value):
        self._mol_type = value

    # seq_type
    def _set_seq_type(self, value):

        if value in MOL_TO_SEQ_TYPE_MAP:
            self._seq_type = MOL_TO_SEQ_TYPE_MAP[value]
        elif value in SEQ_TYPES:
            self._seq_type = value
        else:
            message = ('Molecule type not supported: {s}.')
            message = message.format(s=value)
            raise Exception(message)

    # seq
    @property
    def seq(self):  # noqa
        return self._seq

    def _set_seq(self, value):
        if issubclass(type(value), _Seq):
            self._seq = value
        elif issubclass(type(value), basestring):
            value = value.upper()
            self._seq = Seq(seq=value, seq_type=self._seq_type)

    # accession
    @property
    def accession(self):  # noqa
        return self._accession

    @accession.setter
    def accession(self, value):
        self._accession = value

    # version
    @property
    def version(self):  # noqa
        return self._accession + '.' + str(self._version)

    @version.setter
    def version(self, value):
        if value is not None:
            try:
                self._version = int(value)
            except ValueException:
                print('Version should be an integer.')

    # description (definition)
    @property
    def description(self):  # noqa
        return self._description

    @description.setter
    def description(self, value):
        self._description = value

    # strandedness
    @property
    def strandedness(self):  # noqa
        return self._strandedness

    @strandedness.setter
    def strandedness(self, value):
        self._strandedness = value

    # topology
    @property
    def topology(self):  # noqa
        return self._topology

    @topology.setter
    def topology(self, value):
        self._topology = value

    # division
    @property
    def division(self):  # noqa
        return self._division

    @division.setter
    def division(self, value):
        self._division = value

    # Deal with dates
    def _process_date_string(self, date_str):
        return_value = None
        split_date_str = date_str.split('-')
        if len(split_date_str) != 3:
            raise Exception('Date should be a string of format: YYYY-MM-DD')
        else:
            try:
                return_value = date(
                    int(split_date_str[0]),
                    int(split_date_str[1]),
                    int(split_date_str[2]))
            except Exception as e:
                raise e

        return return_value

    # date_create
    @property
    def date_create(self):  # noqa
        return self._date_create

    @date_create.setter
    def date_create(self, value):
        if value is not None:
            self._date_create = self._process_date_string(date_str=value)

    # date_update
    @property
    def date_update(self):  # noqa
        return self._date_update

    @date_update.setter
    def date_update(self, value):
        if value is not None:
            self._date_update = self._process_date_string(date_str=value)

    # taxid
    @property
    def taxid(self):  # noqa
        return self._taxid

    @taxid.setter
    def taxid(self, value):
        if value is not None:
            try:
                self._taxid = int(value)
            except ValueException:
                print('taxid should be an integer.')

    # organism
    @property
    def organism(self):  # noqa
        return self._organism

    @organism.setter
    def organism(self, value):
        self._organism = value

    # features
    @property
    def features(self):  # noqa
        return self._features

    @features.setter
    def features(self, value):
        self._features = value

    def intervals_for_feature(self, feature=None, qualifier_label=None,  # noqa
                              qualifier_value=None, strict=True, regex=False):

        # Prepare an empty list to hold the returned intervals.
        ints = {'intervals': list(),
                'interval_directions': list()}
        # Iterate over all features.
        for f in self._features:
            if feature is None:
                pass
            # If the current feature is not what we are looking for, continue.
            elif f['key'] != feature:
                continue
            # Otherwise, iterate over all qualifiers.
            for qi in f['qualifiers']:
                ql = None
                # If we do not care for the specific qualifier, consider all
                # qualifiers.
                if qualifier_label is None:
                    ql = list(qi.keys())[0]
                    if not issubclass(type(qi[ql]), basestring):
                        continue
                # If we care about a specific qualifier, see if this is the one
                # we are looking for.
                elif qualifier_label not in qi.keys():
                    continue
                else:
                    ql = qualifier_label
                # If we do not care for the specific qualifier value, add the
                # interval to the list to be returned
                if qualifier_value is None:
                    ints['intervals'].append(f['intervals'])
                    ints['interval_directions'].append(
                        f['interval_directions'])
                else:
                    q = qi[ql]
                    if regex:
                        re_match = match(qualifier_value, q)
                        if re_match is not None:
                            ints['intervals'].append(f['intervals'])
                            ints['interval_directions'].append(
                                f['interval_directions'])
                            break
                    elif strict:
                        if qualifier_value == q:
                            ints['intervals'].append(f['intervals'])
                            ints['interval_directions'].append(
                                f['interval_directions'])
                            break
                    else:
                        if qualifier_value in q:
                            ints['intervals'].append(f['intervals'])
                            ints['interval_directions'].append(
                                f['interval_directions'])
                            break

        return ints
