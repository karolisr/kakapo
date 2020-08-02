from kakapo.tools import seq

from kakapo.tools.transl_tables import TranslationTable

TT = TranslationTable(1)
TRANS_TABLE = TT.table
START_CODONS = TT.start_codons


def test_reverse():
    assert seq.reverse('ACGTAAA') == 'AAATGCA'


def test_complement():
    assert seq.complement('ACGTRYMKWSBDHV') == 'TGCAYRKMWSVHDB'


def test_reverse_complement():
    assert seq.reverse_complement('ACGTRYMKWSBDHV') == 'BDHVSWMKRYACGT'


def test_translate():
    assert seq.translate('ATGGTTAAACCACAACTC',
                         TRANS_TABLE,
                         START_CODONS) == 'MVKPQL'

# ToDo: Test below

# def test_untranslate():
#     assert True
#
#
# def test_seq():
#     assert True
#
#
# def test_ntseq():
#     assert True
#
#
# def test_dnaseq():
#     assert True
#
#
# def test_rnaseq():
#     assert True
#
#
# def test_aaseq():
#     assert True
