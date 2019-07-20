# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import with_statement

import unittest

# from tests import test_data_dir_path

from kakapo.os_diffs import check_os


TEMP_DIR = 'temp'

def setUpModule():
    print('\nsetUpModule kakapoOSDiffsTests')


def tearDownModule():
    print('\n\ntearDownModule kakapoOSDiffsTests')

class kakapoOSDiffsTests(unittest.TestCase):
    print('kakapoOSDiffsTests')

    def test_check_os(self):

        print('\ntest_check_os')

        os_id, os_str, dist_id = check_os()

        print(os_id)
        print(os_str)
        print(dist_id)
