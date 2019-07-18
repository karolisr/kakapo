# -*- coding: utf-8 -*-

from __future__ import print_function

from os.path import join as opj

import unittest

from kakapo.bioio import _parse_gbseq_xml_handle

from tests import test_data_dir_path


def setUpModule():
    print('\nsetUpModule kakapoBioioTests')


def tearDownModule():
    print('\n\ntearDownModule kakapoBioioTests')


class kakapoBioioTests(unittest.TestCase):

    print('kakapoBioioTests')

    def test_parse_gbseq_xml_handle(self):

        print('\ntest_parse_gbseq_xml_handle')

        with open(opj(test_data_dir_path, 'gbseq_xml_sample.xml'), 'r') as f:
            results = _parse_gbseq_xml_handle(f)

        self.assertEqual(len(results), 6)

        for r in results:

            self.assertEqual(len(r.keys()), 14)

            self.assertTrue('accession' in r)
            self.assertTrue('date_create' in r)
            self.assertTrue('date_update' in r)
            self.assertTrue('definition' in r)
            self.assertTrue('division' in r)
            self.assertTrue('features' in r)
            self.assertTrue('length' in r)
            self.assertTrue('mol_type' in r)
            self.assertTrue('organism' in r)
            self.assertTrue('seq' in r)
            self.assertTrue('strandedness' in r)
            self.assertTrue('taxid' in r)
            self.assertTrue('topology' in r)
            self.assertTrue('version' in r)
