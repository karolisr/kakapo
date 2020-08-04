from pprint import pprint

from kakapo.tools import bioio

F1 = 'kakapo_test_data/kakapo_input/F-boxes.fasta'
F2 = 'kakapo_test_data/kakapo_input/T2-RNases.fasta'

T1 = '>NAME1a name1b\nFTIHGLWPDN\n>NAME2a name2b\nFTIHGLWPDNFTIHGLWPDN\n>NAME3a name3b\nWPDN\n'


def test_read_fasta():
    parsed = bioio.read_fasta(F1, seq_type='AA', upper=True,
                              def_to_first_space=False)
    print()
    pprint(parsed)


# def test_dict_to_fasta():
#     parsed = bioio.read_fasta(F1, seq_type='AA', upper=True,
#                               def_to_first_space=False)
#     fasta = bioio.dict_to_fasta(parsed, max_line_len=None)
#     # print()
#     # pprint(fasta)
#
#     fasta = bioio.dict_to_fasta(parsed, max_line_len=60)
#     # print()
#     # pprint(fasta)


# def test_write_fasta():
#     parsed = bioio.read_fasta(F1, seq_type='AA', upper=True,
#                               def_to_first_space=False)
#     bioio.write_fasta(parsed, 'x.fasta', max_line_len=20)


def test_standardize_fasta_text():
    stndr = bioio.standardize_fasta_text(T1, 'AA')

    print()
    pprint(stndr)


def test_trim_desc_to_first_space_in_fasta_text():
    trimmed = bioio.trim_desc_to_first_space_in_fasta_text(T1, 'AA')

    print()
    pprint(trimmed)


def test_filter_fasta_text_by_length():
    filt1 = bioio.filter_fasta_text_by_length(T1, 'AA', 1, 20)
    filt2 = bioio.filter_fasta_text_by_length(T1, 'AA', 9, 12)
    filt3 = bioio.filter_fasta_text_by_length(T1, 'AA', 5, 20)
    print()
    pprint(filt1)
    pprint(filt2)
    pprint(filt3)
