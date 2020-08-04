# from pprint import pprint
#
# from click import open_file
#
# from kakapo.flow.parse_ss import parse_ss_file
#
#
# def test_parse_ss_file():
#     ss_fp = 'kakapo_test_data/kakapo_search_strategies.yaml'
#     ss_file = open_file(ss_fp)
#     print()
#     ss = parse_ss_file(ss_file)
#     pprint(ss)
#     assert 'atpB' in ss
#     assert ss['atpB']['evalue'] == 1e-20
#     # ToDo: more tests
