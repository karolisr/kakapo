# -*- coding: utf-8 -*-

"""seqtk"""

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
