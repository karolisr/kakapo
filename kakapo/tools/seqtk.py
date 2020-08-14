"""seqtk."""

from kakapo.utils.subp import run


def seqtk_fq_to_fa(seqtk, in_file, out_file):
    cmd = [seqtk, 'seq', '-A', in_file]
    # out is stored in memory, could use a lot of RAM
    out = run(cmd, do_not_raise=True)
    with open(out_file, mode='w') as f:
        f.write(out.stdout)


def seqtk_extract_reads(seqtk, in_file, out_file, ids_file):
    cmd = [seqtk, 'subseq', in_file, ids_file]
    # out is stored in memory, could use a lot of RAM
    out = run(cmd, do_not_raise=True)
    with open(out_file, mode='w') as f:
        f.write(out.stdout)


def seqtk_sample_reads(seqtk, in_file, out_file, n, seed=11):
    # n can be a fraction or a number of sequences to sample
    cmd = [seqtk, 'sample', '-2', '-s', str(seed), in_file, n]
    # out is stored in memory, could use a lot of RAM
    out = run(cmd, do_not_raise=True)
    with open(out_file, mode='w') as f:
        f.write(out.stdout)
