"""seq."""

from kakapo.tools.transl_tables import TranslationTable
from kakapo.tools.iupac import AA_AMBIGUOUS
from kakapo.tools.iupac import DNA_AMBIGUOUS
from kakapo.tools.iupac import DNA_COMPLEMENT_TABLE
from kakapo.tools.iupac import DNA_ONLY_CHARS
from kakapo.tools.iupac import NT_AMBIGUOUS
from kakapo.tools.iupac import RNA_AMBIGUOUS
from kakapo.tools.iupac import RNA_ONLY_CHARS

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


def reverse(seq):
    seq = str(seq).upper()
    return seq[::-1]


def complement(seq):
    seq = str(seq).upper()

    seq_contains_uracil = False
    if 'U' in seq:
        seq_contains_uracil = True
        seq = seq.replace('U', 'T')

    seq_complemented = seq.translate(DNA_COMPLEMENT_TABLE)

    if seq_contains_uracil is True:
        seq_complemented = seq_complemented.replace('T', 'U')

    return seq_complemented


def reverse_complement(seq):
    return reverse(complement(seq))


def translate(seq, trans_table, start_codons):
    seq = str(seq).upper()

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


def untranslate(seq, trans_table_inv):
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

    def __new__(cls, seq, seq_type):

        seq_type = seq_type.upper()

        if seq_type in MOL_TO_SEQ_TYPE_MAP:
            seq_type = MOL_TO_SEQ_TYPE_MAP[seq_type]

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
    def __init__(self, seq, seq_type):
        pass


class _Seq(object):
    """
    An abstract class.
    """

    def __init__(self, seq):
        self._seq = seq
        self._gc_code = None
        self._transl_table = None

    def __repr__(self):
        return '_Seq(\'' + self._seq + '\')'

    def __str__(self):
        return self._seq

    @property
    def gc_code(self):
        return self._gc_code

    @gc_code.setter
    def gc_code(self, value):
        if self._gc_code is None:
            self._gc_code = int(value)
            self._transl_table = TranslationTable(self._gc_code)
        else:
            print('Cannot change the value of "gc_code".')

    @property
    def transl_table(self):
        return self._transl_table

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

    def __init__(self, seq):

        self._gc_code = None
        self._transl_table = None

        seq = str(seq).upper()
        if set(seq) <= NT_AMBIGUOUS:
            super(NTSeq, self).__init__(seq)
        else:
            message = 'NT sequence should not contain these characters: {s}.'
            message = message.format(s=', '.join(
                str(s) for s in set(seq) - NT_AMBIGUOUS))
            raise Exception(message)

    def __repr__(self):
        return 'NTSeq(\'' + self._seq + '\')'

    def translate(self, start_codons=tuple('ATG')):
        tt = self._transl_table
        if start_codons is None:
            start_codons = tt.start_codons
        raw = translate(self._seq, tt.table, start_codons)
        return AASeq(raw)

    @property
    def reversed(self):
        return type(self)(reverse(self._seq))

    @property
    def complemented(self):
        return type(self)(complement(self._seq))

    @property
    def reversed_complemented(self):
        return type(self)(reverse_complement(self._seq))


class DNASeq(NTSeq):
    """
    Represents a DNA sequence.

    :param seq: DNA sequence.
    :type seq: str
    """

    def __init__(self, seq):
        seq = str(seq).upper()
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

    def __repr__(self):
        return 'DNASeq(\'' + self._seq + '\')'


class RNASeq(NTSeq):
    """
    Represents an RNA sequence.

    :param seq: RNA sequence.
    :type seq: str
    """

    def __init__(self, seq):
        seq = str(seq).upper()
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

    def __repr__(self):
        return 'RNASeq(\'' + self._seq + '\')'


class AASeq(_Seq):
    """
    Represents an amino acid sequence.

    :param seq: Amino acid sequence.
    :type seq: str
    """

    def __init__(self, seq):
        seq = str(seq).upper()
        if set(seq) <= AA_AMBIGUOUS:
            super(AASeq, self).__init__(seq)
        else:
            message = 'AA sequence should not contain these characters: {s}.'
            message = message.format(s=', '.join(
                str(s) for s in set(seq) - AA_AMBIGUOUS))
            raise Exception(message)

    def __repr__(self):
        return 'AASeq(\'' + self._seq + '\')'


from typing import Union


class SeqRecord(object):
    """
    Represents a sequence record.

    :param seq: Nucleotide or amino acid sequence::class:`NTSeq`,
        :class:`DNASeq`, :class:`RNASeq`, :class:`AASeq`.

    :param definition: Sequence definition
    """

    def __init__(self, definition: str,
                 seq: Union[Seq, NTSeq, DNASeq, RNASeq, AASeq]):
        self._definition = definition
        self._seq = seq
        self._accession = None
        self._version = None
        self._taxid = None
        self._organism = None
        self._features = None
        self._organelle = None
        self._coded_by = None
        self._gc_code = None
        self._transl_table = None

    def __repr__(self):
        name = self.accession_version
        if name is None:
            name = self._definition
        return 'SeqRecord(\'' + name + '\', ' + self._seq.__repr__() + ')'

    def __str__(self):
        return self._seq.__str__()

    @property
    def seq(self):
        return self._seq

    @property
    def length(self):
        return self._seq.length

    @property
    def definition(self):
        return self._definition

    @definition.setter
    def definition(self, value):
        self._definition = value

    @property
    def accession(self):
        return self._accession

    @accession.setter
    def accession(self, value):
        self._accession = value

    @property
    def version(self):
        return int(self._version)

    @version.setter
    def version(self, value):
        self._version = value

    @property
    def accession_version(self):
        version = ''
        if self._accession is None:
            return None
        if self._version is not None:
            version = '.' + self._version
        return self._accession + version

    @property
    def taxid(self):
        return self._taxid

    @taxid.setter
    def taxid(self, value):
        self._taxid = int(value)

    @property
    def organism(self):
        return self._organism

    @organism.setter
    def organism(self, value):
        self._organism = value

    @property
    def features(self):
        return self._features

    @features.setter
    def features(self, value):
        self._features = value

    @property
    def organelle(self):
        return self._organelle

    @organelle.setter
    def organelle(self, value):
        if self._organelle is None:
            self._organelle = value
        else:
            print('Cannot change the value of "organelle".')

    @property
    def coded_by(self):
        return self._coded_by

    @coded_by.setter
    def coded_by(self, value):
        if self._coded_by is None:
            self._coded_by = value
        else:
            print('Cannot change the value of "coded_by".')

    @property
    def gc_code(self):
        return self._seq.gc_code

    @gc_code.setter
    def gc_code(self, value):
        self._seq.gc_code = value

    @property
    def transl_table(self):
        return self._seq.transl_table
