"""SPAdes."""

from os import stat

from kakapo.utils.subp import run
from kakapo.utils.subp import which


PY3 = which('python3')


def run_spades_se(spades, out_dir, input_file, threads, memory, rna):

    memory = str(memory).split('.')[0]

    cmd = [spades,
           '-o', out_dir,
           '-s', input_file,
           '--only-assembler',
           '--threads', str(threads),
           '--memory', memory,
           '--phred-offset', '33']

    if rna:
        cmd.append('--rna')

    cmd = [PY3] + cmd

    run(cmd, do_not_raise=True)


def run_spades_pe(spades, out_dir, input_files, threads, memory, rna):

    memory = str(memory).split('.')[0]

    cmd = [spades,
           '-o', out_dir,
           '--only-assembler',
           '--threads', str(threads),
           '--memory', memory,
           '--phred-offset', '33']

    if stat(input_files[0]).st_size > 512 and stat(input_files[1]).st_size > 512:
        cmd.append('--pe1-1')      # paired_1.fastq
        cmd.append(input_files[0])
        cmd.append('--pe1-2')      # paired_2.fastq
        cmd.append(input_files[1])

    if stat(input_files[2]).st_size > 512:  # unpaired_1.fastq
        cmd.append('--s1')
        cmd.append(input_files[2])

    if stat(input_files[3]).st_size > 512:  # unpaired_2.fastq
        cmd.append('--s2')
        cmd.append(input_files[3])

    if rna:
        cmd.append('--rna')

    cmd = [PY3] + cmd

    run(cmd, do_not_raise=True)
