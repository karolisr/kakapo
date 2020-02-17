# -*- coding: utf-8 -*-

"""SPAdes."""

from kakapo.shell import call


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

    call(cmd)


def run_spades_pe(spades, out_dir, input_files, threads, memory, rna):

    memory = str(memory).split('.')[0]

    cmd = [spades,
           '-o', out_dir,
           '--pe1-1', input_files[0],  # paired_1.fastq
           '--pe1-2', input_files[1],  # paired_2.fastq
           '--s1', input_files[2],     # unpaired_1.fastq
           '--s2', input_files[3],     # unpaired_2.fastq
           '--only-assembler',
           '--threads', str(threads),
           '--memory', memory,
           '--phred-offset', '33']

    if rna:
        cmd.append('--rna')

    call(cmd)
