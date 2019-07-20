# -*- coding: utf-8 -*-
"""Configuration for unit tests for kakapo."""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import with_statement

import sys

from os.path import join
from os.path import dirname
from os.path import abspath

test_data_dir_path = join(dirname(abspath(sys.argv[0])), 'tests', 'data')
