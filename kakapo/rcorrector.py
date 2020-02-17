# -*- coding: utf-8 -*-

"""rcorrector."""

from kakapo.helpers import grouper
from kakapo.helpers import plain_or_gzip
from kakapo.shell import call


def run_rcorrector_se(rcorrector, in_file, out_dir, threads, dir_temp):

    cmd = [rcorrector,
           '-t', str(threads),
           '-s', in_file,
           '-k', '23',
           '-od', out_dir]

    call(cmd, cwd=dir_temp)


def run_rcorrector_pe(rcorrector, in_file_1, in_file_2, out_dir, threads,
                      dir_temp):

    cmd = [rcorrector,
           '-t', str(threads),
           '-1', in_file_1,
           '-2', in_file_2,
           '-k', '23',
           '-od', out_dir]

    call(cmd, cwd=dir_temp)


# The code below is edited from "Transcriptome workshop at Botany 2018"
# Ya Yang and Stephen A. Smith, which was modified from the code in:
# FilterUncorrectabledPEfastq.py file by Adam H. Freedman:
# https://github.com/harvardinformatics/TranscriptomeAssemblyTools


def filter_unc_se(in_file, out_file, log_file=None):

    r_mode, w_mode, a_mode, fqopen, ext = plain_or_gzip(in_file)

    counter = 0
    cor_count = 0
    unc_count = 0

    with fqopen(in_file, r_mode) as in_f, fqopen(out_file, w_mode) as out_f:
        entries = grouper(in_f, 4)
        for entry in entries:
            counter += 1
            # if counter % 100000 == 0:
            #     print('{} reads processed.'.format(counter))
            head, seq, plhld, qual = [i.strip() for i in entry]

            if 'unfixable' in head:
                unc_count += 1

            else:
                if 'cor' in head:
                    cor_count += 1

                # Keep the label information before the Rcorrector flags
                # (low kmer stat, 'cor' and 'unfixable error')
                head = head.split('l:')[0][:-1]
                entry_corrected = '\n'.join([head, seq, plhld, qual])
                out_f.write('{}\n'.format(entry_corrected))

    if log_file is not None:

        log_str = ('Total SE reads: {}\n'
                   'Uncorrectable SE reads removed: {} - {:.2f}%\n'
                   'Corrected SE reads retained: {} - {:.2f}%\n'
                   'Total SE reads retained: {} - {:.2f}%\n')

        with open(log_file, 'w') as f:
            f.write(log_str.format(
                counter,
                unc_count, (unc_count / counter) * 100,
                cor_count, (cor_count / counter) * 100,
                counter - unc_count, ((counter - unc_count) / counter) * 100))


def filter_unc_pe(in_file_1, in_file_2, out_file_1, out_file_2, log_file=None):

    r_mode, w_mode, a_mode, fqopen, ext = plain_or_gzip(in_file_1)

    counter = 0
    cor_1_count = 0
    cor_2_count = 0
    cor_12_count = 0
    unc_count = 0

    with fqopen(in_file_1, r_mode) as in_f_1, \
            fqopen(in_file_2, r_mode) as in_f_2, \
            fqopen(out_file_1, w_mode) as out_f_1, \
            fqopen(out_file_2, w_mode) as out_f_2:

        entries_1 = grouper(in_f_1, 4)
        entries_2 = grouper(in_f_2, 4)

        for entry_1 in entries_1:
            entry_2 = next(entries_2)
            counter += 1
            # if counter % 100000 == 0:
            #     print('{} reads processed.'.format(counter))
            head_1, seq_1, plhld_1, qual_1 = [i.strip() for i in entry_1]
            head_2, seq_2, plhld_2, qual_2 = [j.strip() for j in entry_2]

            if 'unfixable' in head_1 or 'unfixable' in head_2:
                unc_count += 1
            else:
                if 'cor' in head_1 and 'cor' in head_2:
                    cor_12_count += 1
                else:
                    if 'cor' in head_1:
                        cor_1_count += 1
                    elif 'cor' in head_2:
                        cor_2_count += 1

                # Keep the label information before the Rcorrector flags
                # (low kmer stat, 'cor' and 'unfixable error')
                head_1 = head_1.split('l:')[0][:-1]
                entry_1_corrected = '\n'.join([head_1, seq_1, plhld_1, qual_1])
                out_f_1.write('{}\n'.format(entry_1_corrected))

                head_2 = head_2.split('l:')[0][:-1]
                entry_2_corrected = '\n'.join([head_2, seq_2, plhld_2, qual_2])
                out_f_2.write('{}\n'.format(entry_2_corrected))

    if log_file is not None:

        log_str = ('Total PE reads: {}\n'
                   'Uncorrectable PE reads removed: {} - {:.2f}%\n'
                   'Corrected R1 reads retained: {} - {:.2f}%\n'
                   'Corrected R2 reads retained: {} - {:.2f}%\n'
                   'Corrected PE reads retained: {} - {:.2f}%\n'
                   'Total PE reads retained: {} - {:.2f}%\n')

        with open(log_file, 'w') as f:
            f.write(log_str.format(
                counter,
                unc_count, (unc_count / counter) * 100,
                cor_1_count, (cor_1_count / counter) * 100,
                cor_2_count, (cor_2_count / counter) * 100,
                cor_12_count, (cor_12_count / counter) * 100,
                counter - unc_count, ((counter - unc_count) / counter) * 100))
