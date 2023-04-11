"""bowtie2."""

from os import remove
from os.path import join as opj

from kakapo.utils.misc import split_mixed_fq
from kakapo.utils.subp import run


def build_bt2_index(bowtie2_build, input_files, output_path, threads):
    cmd = [bowtie2_build, '--threads', str(threads),
           ','.join(input_files), output_path]

    run(cmd, do_not_raise=True)


def run_bowtie2_se(bowtie2, input_file,
                   output_file, output_file_un, sam_output_file,
                   index, threads, dir_temp, gz_out=True):

    par_al = '--al'
    par_un = '--un'
    if gz_out is True:
        par_al = '--al-gz'
        par_un = '--un-gz'

    cmd = [bowtie2, '--threads', str(threads), '--very-sensitive',
           '--phred33', '--no-unal', '--no-mixed', '--no-discordant',
           # '-a',
           '--rdg', '1000,1000',
           '--rfg', '1000,1000',
           '-x', index,
           '-U', input_file,
           par_al, output_file,
           par_un, output_file_un,
           '-S', sam_output_file]

    run(cmd, cwd=dir_temp, do_not_raise=True)


def run_bowtie2_pe(bowtie2, input_files, paired_out_pattern,
                   paired_out_pattern_un,
                   unpaired_out_1, unpaired_out_2,
                   unpaired_out_1_un, unpaired_out_2_un, sam_output_file,
                   index, threads, dir_temp, gz_out=True):

    par_al_conc = '--al-conc'
    par_un_conc = '--un-conc'
    par_al = '--al'
    par_un = '--un'
    ext_out = ''
    if gz_out is True:
        par_al_conc = '--al-conc-gz'
        par_un_conc = '--un-conc-gz'
        par_al = '--al-gz'
        par_un = '--un-gz'
        ext_out = '.gz'

    temp_unpaired_file = opj(dir_temp, 'temp_unpaired.fastq' + ext_out)
    temp_unpaired_file_un = opj(dir_temp, 'temp_unpaired_un.fastq' + ext_out)

    cmd = [bowtie2, '--threads', str(threads), '--very-sensitive',
           '--phred33', '--no-unal', '--no-mixed', '--no-discordant',
           # '-a',
           '--rdg', '1000,1000',
           '--rfg', '1000,1000',
           '-x', index,
           '-1', input_files[0], '-2', input_files[1],
           '-U', input_files[2] + ',' + input_files[3],
           par_al_conc, paired_out_pattern, par_al, temp_unpaired_file,
           par_un_conc, paired_out_pattern_un, par_un, temp_unpaired_file_un,
           '-S', sam_output_file]

    run(cmd, cwd=dir_temp, do_not_raise=True)
    split_mixed_fq(temp_unpaired_file, unpaired_out_1, unpaired_out_2)
    remove(temp_unpaired_file)
    split_mixed_fq(temp_unpaired_file_un, unpaired_out_1_un, unpaired_out_2_un)
    remove(temp_unpaired_file_un)


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
