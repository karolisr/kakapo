# -*- coding: utf-8 -*-

"""Manages Python 2 and Python 3 differences."""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import with_statement

import io
import re
import sys


def python_version():
    """
    Determine the Python version.
    """
    py_v_hex = sys.hexversion

    py_v_1 = sys.version_info[0]
    py_v_2 = sys.version_info[1]
    py_v_3 = sys.version_info[2]

    py_v_str = '{v0}.{v1}.{v2}'.format(v0=py_v_1, v1=py_v_2, v2=py_v_3)

    return py_v_hex, py_v_str


_py_v_hex, _ = python_version()

if _py_v_hex >= 0x03000000:
    from configparser import ConfigParser as _CP # noqa
    from io import StringIO # noqa
    from functools import singledispatch # noqa

    maketrans = str.maketrans

    unicode = str
    str = str
    bytes = bytes
    basestring = (str, bytes)

    file = io.IOBase
    HANDLE_TYPES = (file, io.IOBase, StringIO)

elif _py_v_hex < 0x03000000:
    from string import maketrans # noqa
    from ConfigParser import SafeConfigParser as _CP # noqa
    from StringIO import StringIO # noqa
    from singledispatch import singledispatch # noqa

    unicode = unicode
    bytes = str
    str = unicode
    basestring = basestring

    HANDLE_TYPES = (file, io.IOBase, StringIO) # noqa


class ConfigParser(_CP):
    """
    ConfigParser that does not allow ':'. as a separator.

    This is not necessary in Python 3 because it is possible to set the
    separators as optional arguments. In Python 2, this is impossible.
    """

    # OPTCRE = re.compile(
    #     r'(?P<option>[^:=\s][^:=]*)'          # very permissive!
    #     r'\s*(?P<vi>[:=])\s*'                 # any number of space/tab,
    #                                           # followed by separator
    #                                           # (either : or =), followed
    #                                           # by any # space/tab
    #     r'(?P<value>.*)$'                     # everything up to eol
    #     )

    # OPTCRE_NV = re.compile(
    #     r'(?P<option>[^:=\s][^:=]*)'          # very permissive!
    #     r'\s*(?:'                             # any number of space/tab,
    #     r'(?P<vi>[:=])\s*'                    # optionally followed by
    #                                           # separator (either : or
    #                                           # =), followed by any #
    #                                           # space/tab
    #     r'(?P<value>.*))?$'                   # everything up to eol
    #     )

    OPTCRE_NV = re.compile(
        r'(?P<option>[^=\s][^=]*)'              # very permissive!
        r'\s*(?:'                               # any number of space/tab,
        r'(?P<vi>[=])\s*'                       # optionally followed by
                                                # separator (=),
                                                # followed by any #
                                                # space/tab
        r'(?P<value>.*))?$')                    # everything up to eol

    OPTCRE = OPTCRE_NV
