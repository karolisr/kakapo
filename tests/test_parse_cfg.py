# from pprint import pprint
#
# from click import open_file
#
# from kakapo.flow.parse_cfg import parse_config_file
#
#
# class Taxonomy:
#     """A dummy class for testing only."""
#
#     @classmethod
#     def taxid_valid(cls, taxid):
#         return True
#
#     @classmethod
#     def tax_id_for_name_and_group_tax_id(cls, name, group_tax_id):
#         return 0
#
#
# def test_parse_config_file():
#     cfg_fp = 'kakapo_test_data/kakapo_config.yaml'
#     config_file = open_file(cfg_fp)
#     print()
#     config = parse_config_file(config_file, Taxonomy)
#     pprint(config)
#     assert config['tax_ids'] == (0,)
#     assert config['group_tax_id'] == 33090
#     # ToDo: more tests
