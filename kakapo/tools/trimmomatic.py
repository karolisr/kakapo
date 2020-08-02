"""Trimmomatic."""

from kakapo.utils.subp import run


def trimmomatic_se(trimmomatic, adapters, in_file, out_file, stats_file,
                   threads, minlen):

    cmd = ['java', '-jar', trimmomatic, 'SE', '-threads', str(threads),
           '-phred33',
           '-summary', stats_file, in_file, out_file,
           'ILLUMINACLIP:' + adapters + ':2:30:10',
           'SLIDINGWINDOW:4:20',
           'LEADING:20',
           'TRAILING:20',
           'MINLEN:' + str(minlen)]

    run(cmd)


def trimmomatic_pe(trimmomatic, adapters, in_file_1, in_file_2,
                   out_file_paired_1, out_file_unpaired_1,
                   out_file_paired_2, out_file_unpaired_2,
                   stats_file, threads, minlen):

    cmd = ['java', '-jar', trimmomatic, 'PE', '-threads', str(threads),
           '-phred33',
           '-summary', stats_file,
           in_file_1, in_file_2,
           out_file_paired_1, out_file_unpaired_1,
           out_file_paired_2, out_file_unpaired_2,
           'ILLUMINACLIP:' + adapters + ':2:30:10',
           'SLIDINGWINDOW:4:20',
           'LEADING:20',
           'TRAILING:20',
           'MINLEN:' + str(minlen)]

    run(cmd)
