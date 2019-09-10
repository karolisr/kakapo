# -*- coding: utf-8 -*-

"""bowtie2"""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import with_statement

from os import remove
from os.path import join as opj

from kakapo.helpers import split_mixed_fq
from kakapo.shell import call


def build_bt2_index(bowtie2_build, input_files, output_path, threads):  # noqa
    cmd = [bowtie2_build, '--threads', str(threads),
           ','.join(input_files), output_path]

    call(cmd)


def run_bowtie2_se(bowtie2, input_file,
                   output_file, sam_output_file,
                   index, threads, dir_temp):  # noqa

    cmd = [bowtie2, '--threads', str(threads), '--very-sensitive-local',
           '--phred33', '--no-unal', '--no-mixed', '--no-discordant',
           '-x', index,
           '-U', input_file,
           '--al', output_file,
           '-S', sam_output_file]

    call(cmd, cwd=dir_temp)


def run_bowtie2_pe(bowtie2, input_files, paired_out_pattern,
                   unpaired_out_1, unpaired_out_2, sam_output_file,
                   index, threads, dir_temp):  # noqa

    temp_unpaired_file = opj(dir_temp, 'temp_unpaired.fastq')

    cmd = [bowtie2, '--threads', str(threads), '--very-sensitive-local',
           '--phred33', '--no-unal', '--no-mixed', '--no-discordant',
           '-x', index,
           '-1', input_files[0], '-2', input_files[1],
           '-U', input_files[2] + ',' + input_files[3],
           '--al-conc', paired_out_pattern, '--al', temp_unpaired_file,
           '-S', sam_output_file]

    call(cmd, cwd=dir_temp)
    split_mixed_fq(temp_unpaired_file, unpaired_out_1, unpaired_out_2)
    remove(temp_unpaired_file)


# --threads 4
# --minins 0
# --maxins 700
# --very-sensitive-local
# --phred33
# --ma 3
# --mp 7
# --np 2
# --rdg 6,4
# --rfg 6,4
# --score-min G,20,8
# --no-unal
# --no-mixed
# --no-discordant
