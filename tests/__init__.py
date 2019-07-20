# -*- coding: utf-8 -*-
"""Unit tests for kakapo."""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import with_statement

from tests.config import test_data_dir_path

from tests.kakapo_seq import kakapoSeqTests
from tests.kakapo_bioio import kakapoBioioTests
from tests.kakapo_entrez import kakapoEntrezTests
from tests.kakapo_workflow import kakapoWorkflowTests

__all__ = ['test_data_dir_path', 'kakapoSeqTests', 'kakapoBioioTests',
           'kakapoEntrezTests', 'kakapoWorkflowTests']
