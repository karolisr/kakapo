#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import unittest
import datetime
from tests import *

from kakapo.py_v_diffs import py_v_str
py_ver_msg = '\nPython version: {pv}\n'.format(pv=py_v_str)
print(py_ver_msg)


def main():

    unittest.main()

    start_time = datetime.datetime.now()
    end_time = datetime.datetime.now()
    time_taken = end_time - start_time
    print('Time taken by tests:', time_taken)


if __name__ == '__main__':
    main()
