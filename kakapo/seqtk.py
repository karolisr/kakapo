# -*- coding: utf-8 -*-

"""seqtk"""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import with_statement

from kakapo.shell import call


def seqtk_fq_to_fa(seqtk, in_file, out_file):  # noqa
    cmd = [seqtk, 'seq', '-A', in_file]
    # out is stored in memory, could use a lot of RAM
    out, err = call(cmd)
    with open(out_file, mode='wb') as f:
        f.write(out)
    call(cmd)


def seqtk_extract_reads(seqtk, in_file, out_file, ids_file):  # noqa
    cmd = [seqtk, 'subseq', in_file, ids_file]
    # out is stored in memory, could use a lot of RAM
    out, err = call(cmd)
    with open(out_file, mode='wb') as f:
        f.write(out)
    call(cmd)
