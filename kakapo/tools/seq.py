"""seq."""

from copy import deepcopy
from typing import Union, List, Iterable, Dict

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


def translate(seq, trans_table, start_codons, strip_stop_codon=False):
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

    if strip_stop_codon is True:
        translated = translated.strip('*')

    return translated


def untranslate(seq, trans_table_inv):
    seq = str(seq).upper()
    tt = trans_table_inv
    tt['X'] = ('NNN',)
    unt = ''.join(tuple(map(lambda x: tt[x][0], seq)))
    return unt


def parse_qualifiers(qualifiers: Iterable[Dict]) -> Dict:
    return_value = dict()
    for d in qualifiers:
        ks = tuple(d.keys())
        vs = tuple(d.values())
        assert len(ks) == 1
        assert len(vs) == 1
        k = ks[0]
        v = vs[0]
        return_value.update({k: v})
    return return_value


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

    def __new__(cls, seq, seq_type, gc_id=None):

        seq_type = seq_type.upper()

        if seq_type in MOL_TO_SEQ_TYPE_MAP:
            seq_type = MOL_TO_SEQ_TYPE_MAP[seq_type]

        if seq_type in SEQ_TYPES:

            if seq_type is SEQ_TYPE_NT:
                return NTSeq(seq=seq, gc_id=gc_id)
            elif seq_type is SEQ_TYPE_DNA:
                return DNASeq(seq=seq, gc_id=gc_id)
            elif seq_type is SEQ_TYPE_RNA:
                return RNASeq(seq=seq, gc_id=gc_id)
            elif seq_type is SEQ_TYPE_AA:
                return AASeq(seq=seq, gc_id=gc_id)

        else:
            message = 'Invalid sequence type: {}'.format(seq_type)
            raise Exception(message)

    # __init__ declaration is exactly the same as __new__ so Sphinx
    # docstring parser picks it up.
    def __init__(self, seq, seq_type, gc_id):
        pass


class _Seq(object):
    """An abstract sequence class."""

    def __init__(self, seq, gc_id=None):
        self._seq = seq
        self._transl_table = None
        self.gc_id = gc_id

    def __repr__(self):
        return '_Seq(\'' + self._seq + '\')'

    def __str__(self):
        return self._seq

    def __len__(self):
        return self.length

    def __add__(self, other):
        return type(self)(str(self) + str(other))

    def __lt__(self, other):
        return self.length < other.length

    def __getitem__(self, key):
        s = str(self._seq)
        if isinstance(key, int) and key >= 0:
            return type(self)(s[key:key + 1])
        elif isinstance(key, slice):
            return type(self)(s[key.start:key.stop:key.step])
        else:
            raise KeyError('Keys must be integers or slices, '
                           'not {}'.format(type(key)))

    @property
    def gc_id(self):
        if self.transl_table is None:
            print('No genetic code ID (gc_id) was set for this '
                  'sequence object.')
            return None
        return self.transl_table.gc_id

    @gc_id.setter
    def gc_id(self, value):
        if self._transl_table is None:
            if value is None:
                self._transl_table = None
            else:
                self._transl_table = TranslationTable(int(value))
        else:
            print('Cannot change the value of "gc_id".')

    @property
    def transl_table(self):
        if self._transl_table is None:
            print('No genetic code ID (gc_id) was set for this sequence object'
                  '. Therefore, translation table is not present.')
        return self._transl_table

    @property
    def length(self):
        return len(self._seq)


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

    def __init__(self, seq, gc_id=None):
        seq = str(seq).upper()
        if set(seq) <= NT_AMBIGUOUS:
            super(NTSeq, self).__init__(seq=seq, gc_id=gc_id)
        else:
            message = 'NT sequence should not contain these characters: {s}.'
            message = message.format(s=', '.join(
                str(s) for s in set(seq) - NT_AMBIGUOUS))
            raise Exception(message)

    def __repr__(self):
        return 'NTSeq(\'' + self._seq + '\')'

    def translate(self, start_codons=None, strip_stop_codon=False):
        if self.transl_table is None:
            translated = None
        else:
            tt = self.transl_table
            if start_codons is None:
                start_codons = tt.start_codons
            translated = AASeq(
                translate(self._seq, tt.table, start_codons, strip_stop_codon), tt.gc_id)
        return translated

    def untranslate(self):
        return self

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

    def __init__(self, seq, gc_id=None):
        seq = str(seq).upper()
        if set(seq) <= DNA_AMBIGUOUS:
            super(DNASeq, self).__init__(seq=seq, gc_id=gc_id)
        elif set(seq) - DNA_AMBIGUOUS == RNA_ONLY_CHARS:
            seq = seq.replace('U', 'T')
            super(DNASeq, self).__init__(seq=seq, gc_id=gc_id)
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

    def __init__(self, seq, gc_id=None):
        seq = str(seq).upper()
        if set(seq) <= RNA_AMBIGUOUS:
            super(RNASeq, self).__init__(seq=seq, gc_id=gc_id)
        elif set(seq) - RNA_AMBIGUOUS == DNA_ONLY_CHARS:
            seq = seq.replace('T', 'U')
            super(RNASeq, self).__init__(seq=seq, gc_id=gc_id)
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

    def __init__(self, seq, gc_id=None):
        seq = str(seq).upper()
        if set(seq) <= AA_AMBIGUOUS:
            super(AASeq, self).__init__(seq=seq, gc_id=gc_id)
        else:
            message = 'AA sequence should not contain these characters: {s}.'
            message = message.format(s=', '.join(
                str(s) for s in set(seq) - AA_AMBIGUOUS))
            raise Exception(message)

    def __repr__(self):
        return 'AASeq(\'' + self._seq + '\')'

    def translate(self):
        return self

    def untranslate(self):
        if self.transl_table is None:
            untranslated = None
        else:
            tt = self.transl_table
            untranslated = NTSeq(untranslate(self._seq, tt.table_inv),
                                 tt.gc_id)
        return untranslated


class SeqRecord(object):
    """
    Represents a sequence record.

    :param seq: Nucleotide or amino acid sequence::class:`NTSeq`,
        :class:`DNASeq`, :class:`RNASeq`, :class:`AASeq`.

    :param definition: Sequence definition
    """

    def __init__(self, definition: str,
                 seq: Union[NTSeq, DNASeq, RNASeq, AASeq]):
        self._definition = definition
        self._seq = seq
        self._accession = None
        self._version = None
        self._taxid = None
        self._organism = None
        self._features = None
        self._organelle = None
        self._coded_by = None
        self._parent = None

    def __repr__(self):
        name = self.accession_version
        if name is None:
            name = self._definition
        return 'SeqRecord(\'' + name + '\', ' + self._seq.__repr__() + ')'

    def __str__(self):
        return self._seq.__str__()

    def __len__(self):
        return self.length

    def __lt__(self, other):
        return self.length < other.length

    def __getitem__(self, key):
        s = str(self._seq)
        if isinstance(key, int) and key >= 0:
            return type(self._seq)(s[key:key + 1])
        elif isinstance(key, slice):
            return type(self._seq)(s[key.start:key.stop:key.step])
        else:
            raise KeyError('Keys must be integers or slices, '
                           'not {}'.format(type(key)))

    def translate(self, start_codons=None, strip_stop_codon=False):
        return self.seq.translate(start_codons=start_codons, strip_stop_codon=strip_stop_codon)

    def untranslate(self):
        return self.seq.untranslate()

    @property
    def seq(self):
        return self._seq

    @seq.setter
    def seq(self, value):
        self._seq = value

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
        return self._version

    @version.setter
    def version(self, value):
        if value is not None:
            self._version = int(value)

    @property
    def accession_version(self):
        version = ''
        if self._accession is None:
            return None
        if self._version is not None:
            version = '.' + str(self._version)
        return self._accession + version

    @property
    def taxid(self):
        return self._taxid

    @taxid.setter
    def taxid(self, value):
        if self._taxid is None:
            self._taxid = int(value)
        else:
            print('Cannot change the value of "taxid".')

    @property
    def organism(self):
        return self._organism

    @organism.setter
    def organism(self, value):
        if self._organism is None:
            self._organism = value
        else:
            print('Cannot change the value of "organism".')

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
    def parent(self):
        return self._parent

    @parent.setter
    def parent(self, value):
        if self._parent is None:
            self._parent = value
        else:
            print('Cannot change the value of "parent".')

    @property
    def cds(self):
        return cds(self)


class SeqRecordCDS(SeqRecord):
    """Represents a CDS sequence record."""

    def __init__(self, seq: Union[NTSeq, DNASeq, RNASeq]):
        self._gene = None
        self._product = None
        self._parent_intervals = None
        self._parent_interval_directions = None
        self._parent_codon_start = None

        super(SeqRecordCDS, self).__init__(definition='DEFINITION', seq=seq)

    def __repr__(self):
        name = self.accession_version
        if name is None:
            name = self.definition
        return 'SeqRecordCDS(\'' + name + '\', ' + self._seq.__repr__() + ')'

    @property
    def gene(self):
        return self._gene

    @gene.setter
    def gene(self, value):
        self._gene = value

    @property
    def product(self):
        return self._product

    @product.setter
    def product(self, value):
        self._product = value

    @property
    def parent_intervals(self):
        return self._parent_intervals

    @parent_intervals.setter
    def parent_intervals(self, value):
        self._parent_intervals = value

    @property
    def parent_interval_directions(self):
        return self._parent_interval_directions

    @parent_interval_directions.setter
    def parent_interval_directions(self, value):
        self._parent_interval_directions = value

    @property
    def parent_codon_start(self):
        return self._parent_codon_start

    @parent_codon_start.setter
    def parent_codon_start(self, value):
        self._parent_codon_start = value

    @property
    def definition(self):
        name = self.gene
        if name is None:
            name = self.product
        # definition = f'{self.organism.replace(" ", "_")} {self.organism.replace(" ", "_")} {name} located on {self.parent.accession_version} at {str(self.parent_intervals)} {str(self.parent_interval_directions)} {self.organelle}'
        definition = f'{name}||{self.organism.replace(" ", "_")}||{self.parent.accession_version}||{str(self.parent_intervals).replace(", ", "_")}||{str(self.parent_interval_directions).replace(", ", "_")}||{self.organelle}'
        return definition


def cds(sr: SeqRecord) -> List[SeqRecordCDS]:
    # K:
    # print(sr.organism, sr.accession_version)
    return_value: List[SeqRecordCDS] = list()
    fts = deepcopy(sr.features)
    fts_cds = [x for x in fts if x['key'] == 'CDS']
    for ft_cds in fts_cds:

        # K:
        # print(ft_cds)

        qualifiers = parse_qualifiers(ft_cds['qualifiers'])
        codon_start = int(qualifiers['codon_start']) - 1
        intervals = tuple(
            zip(ft_cds['intervals'], ft_cds['interval_directions']))

        cds_seq = DNASeq('', sr.seq.gc_id)

        for interval in intervals:
            # K:
            # print(interval)
            b = interval[0][0]
            e = interval[0][1]
            direction = interval[1]
            if direction == -1:
                b = interval[0][1]
                e = interval[0][0]
            cds_seq_temp = sr.seq[b:e]
            if direction == -1:
                cds_seq_temp = cds_seq_temp.reversed_complemented
            cds_seq += cds_seq_temp

        cds_seq = cds_seq[codon_start:]
        cds_seq.gc_id = sr.seq.gc_id

        # translation = cds_seq.translate(start_codons=None)
        # translation_annotated = qualifiers['translation']
        # mRNA editing will fuck up the test below:
        # assert str(translation[:-1]) == translation_annotated

        product = ''
        if 'product' in qualifiers:
            product = qualifiers['product']

        gene = ''
        if 'gene' in qualifiers:
            gene = qualifiers['gene']

        # note = ''
        # if 'note' in qualifiers:
        #     note = qualifiers['note']

        # trans_splicing = ''
        # if 'trans_splicing' in qualifiers:
        #     trans_splicing = qualifiers['trans_splicing']

        # locus_tag = qualifiers['locus_tag']
        # protein_id = qualifiers['protein_id']
        # db_xref = qualifiers['db_xref']

        ft_intervals = ft_cds['intervals']
        ft_interval_directions = ft_cds['interval_directions']

        ft_intervals_new = list()
        fti_start = 0
        for i, fti in enumerate(ft_intervals):
            fti_len = abs(fti[0] - fti[1])
            if i == 0:
                fti_len -= codon_start
            fti_new = [fti_start, fti_len + fti_start]
            ft_intervals_new.append(fti_new)
            fti_start += fti_len

        ft_interval_directions_new = [abs(x) for x in ft_interval_directions]

        ft_cds['intervals'] = ft_intervals_new
        ft_cds['interval_directions'] = ft_interval_directions_new

        # definition = f'{sr.accession_version}||{str(ft_intervals)}||{str(ft_interval_directions)}||{locus_tag}||{protein_id}||{db_xref}||{product}||{gene}||{note}||{trans_splicing}'

        # name = gene
        # if name is None:
        #     name = product
        # definition = f'{name} located on {sr.accession_version} at {str(ft_intervals)} {str(ft_interval_directions)}'

        # definition = ''
        # if sr.accession_version is None:
        #     if sr.parent is not None:
        #         definition = sr.definition

        cds_sr = SeqRecordCDS(seq=cds_seq)

        cds_sr.taxid = sr.taxid
        cds_sr.organism = sr.organism
        cds_sr.organelle = sr.organelle
        cds_sr.gene = gene
        cds_sr.product = product
        cds_sr.parent = sr
        cds_sr.parent_intervals = ft_intervals
        cds_sr.parent_interval_directions = ft_interval_directions
        cds_sr.parent_codon_start = codon_start

        cds_sr.features = [ft_cds]
        return_value.append(cds_sr)
    return return_value
