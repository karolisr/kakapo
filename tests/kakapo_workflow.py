# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import with_statement

import unittest

from tempfile import TemporaryDirectory

# from tests import test_data_dir_path

from kakapo.workflow import dnld_sra_info
from kakapo.workflow import pfam_uniprot_accessions


def setUpModule():
    print('\nsetUpModule kakapoWorkflowTests')


def tearDownModule():
    print('\n\ntearDownModule kakapoWorkflowTests')


class kakapoWorkflowTests(unittest.TestCase):
    print('kakapoWorkflowTests')

    def test_dnld_sra_info(self):

        print('\ntest_dnld_sra_info')

        sras = ['SRR4099901', 'SRR1985008', 'SRR3087077', 'SRR4434436']
        with TemporaryDirectory() as temp_dir:
            sra_runs_info = dnld_sra_info(sras=sras, dir_cache_prj=temp_dir)
        self.assertEqual(len(sra_runs_info), 4)

    def test_pfam_uniprot_accessions(self):

        pfam_acc = ['PF00445', 'PF00992']
        tax_ids = [3702, 3701, 60711]
        expected_results = ['P42813', 'P42814', 'P42815', 'Q9XI64', 'F4HUG9',
                            'F4IW05', 'A0A0D9RI88', 'A0A0D9RBV6', 'A0A0D9RBX1',
                            'A0A0D9RN56', 'A0A0D9RNF3', 'A0A0D9S0D5',
                            'A0A0D9S5X7', 'A0A0D9S5X8', 'A0A0D9S5X9']

        with TemporaryDirectory() as temp_dir:
            pfam_uniprot_acc = pfam_uniprot_accessions(
                pfam_acc=pfam_acc,
                tax_ids=tax_ids,
                dir_cache_pfam_acc=temp_dir)

            self.assertEqual(pfam_uniprot_acc, expected_results)
