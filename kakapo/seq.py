# -*- coding: utf-8 -*-
"""seq"""

from kakapo.iupac import AA_AMBIGUOUS
from kakapo.iupac import DNA_AMBIGUOUS
from kakapo.iupac import DNA_COMPLEMENT_TABLE
from kakapo.iupac import DNA_ONLY_CHARS
from kakapo.iupac import NT_AMBIGUOUS
from kakapo.iupac import RNA_AMBIGUOUS
from kakapo.iupac import RNA_ONLY_CHARS

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


def reverse(seq):  # noqa
    seq = seq.upper()
    return seq[::-1]


def complement(seq):  # noqa
    seq = seq.upper()

    seq_contains_uracil = False
    if 'U' in seq:
        seq_contains_uracil = True
        seq = seq.replace('U', 'T')

    seq_complemented = seq.translate(DNA_COMPLEMENT_TABLE)

    if seq_contains_uracil is True:
        seq_complemented = seq_complemented.replace('T', 'U')

    return seq_complemented


def reverse_complement(seq):  # noqa
    return reverse(complement(seq))


def translate(seq, trans_table, start_codons):  # noqa
    seq = seq.upper()

    if 'U' in seq:
        seq = seq.replace('U', 'T')

    codons = [(seq[i:i + 3]) for i in range(0, len(seq), 3)]

    # Clip the last item if it consists of less than three nucleotides
    codons = [cod for cod in codons if len(cod) == 3]

    idx = 0
    codon_count = len(codons)
    translated = ''

    if codons[0] in start_codons:
        translated = 'M'
        idx = 1

    for t in codons[idx:codon_count]:
        residue = 'X'
        if t in trans_table:
            residue = trans_table[t]

        translated = translated + residue

    return translated


def untranslate(seq, trans_table_inv):  # noqa
    seq = seq.upper()
    tt = trans_table_inv
    tt['X'] = ('NNN',)
    unt = ''.join(tuple(map(lambda x: tt[x][0], seq)))
    return unt


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

        seq_type = seq_type.upper()

        if seq_type in SEQ_TYPES:

            if seq_type is SEQ_TYPE_NT:
                return NTSeq(seq=seq)
            elif seq_type is SEQ_TYPE_DNA:
                return DNASeq(seq=seq)
            elif seq_type is SEQ_TYPE_RNA:
                return RNASeq(seq=seq)
            elif seq_type is SEQ_TYPE_AA:
                return AASeq(seq=seq)

        else:
            message = 'Invalid sequence type: {}'.format(seq_type)
            raise Exception(message)

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
        seq = seq.upper()
        if set(seq) <= NT_AMBIGUOUS:
            super(NTSeq, self).__init__(seq)
        else:
            message = ('NT sequence should not contain these characters: {s}.')
            message = message.format(s=', '.join(
                str(s) for s in set(seq) - NT_AMBIGUOUS))
            raise Exception(message)

    def translate(self, trans_table, start_codons):  # noqa
        raw = translate(self._seq, trans_table, start_codons)
        return AASeq(raw)

    def reverse(self):  # noqa
        return type(self)(reverse(self._seq))

    def complement(self):  # noqa
        return type(self)(complement(self._seq))

    def reverse_complement(self):  # noqa
        return type(self)(reverse_complement(self._seq))


class DNASeq(NTSeq):
    """

    Represents DNA sequence.

    :param seq: DNA sequence.
    :type seq: str

    """

    def __init__(self, seq):  # noqa
        seq = seq.upper()
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
        seq = seq.upper()
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
        seq = seq.upper()
        if set(seq) <= AA_AMBIGUOUS:
            super(AASeq, self).__init__(seq)
        else:
            message = ('AA sequence should not contain these characters: {s}.')
            message = message.format(s=', '.join(
                str(s) for s in set(seq) - AA_AMBIGUOUS))
            raise Exception(message)
