# -*- coding: utf-8 -*-

"""
Calls to other programs.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import with_statement

from subprocess import PIPE
from subprocess import Popen


def call(cmd, stdout=PIPE, stderr=PIPE, cwd=None):  # noqa
    from kakapo.config import DEBUG_MODE, DEBUG_PROCESSES
    if DEBUG_MODE and DEBUG_PROCESSES:
        stdout = None
        stderr = None

    p = Popen(cmd, stdout=stdout, stderr=stderr, cwd=cwd)
    out, err = p.communicate()

    return out, err
