# -*- coding: utf-8 -*-

"""kraken2"""

from collections import OrderedDict
from os import remove, rename, stat
from os.path import join as opj
from os.path import splitext, basename
from shutil import move

from kakapo.config import RAM
from kakapo.helpers import splitext_gz, plain_or_gzip, make_dir
from kakapo.shell import call


def _use_memory_mapping(db_path, linfo=print):
    db_size = stat(opj(db_path, 'hash.k2d')).st_size / (1024 ** 3)
    mem_max = RAM / 3
    if mem_max < db_size:
        db_name = splitext(basename(db_path))[0]
        linfo('Not enough memory for Kraken2 database {}. '
              'Switching to a slower memory-mapping mode.'.format(db_name))
        return '--memory-mapping'
    else:
        return None


def run_kraken_se(kraken, db, in_file, out_class_file, out_unclass_file,
                  report_file, confidence, threads, dir_temp, linfo=print):  # noqa

    cmd = [kraken, '--db', db, '--threads', str(threads),
           '--confidence', str(confidence), '--output', '-',
           '--report', report_file, '--use-names',
           '--classified-out', out_class_file,
           '--unclassified-out', out_unclass_file, in_file]

    mm = _use_memory_mapping(db, linfo=linfo)
    if mm is not None:
        cmd.insert(1, mm)

    _, _, ext = splitext_gz(in_file)
    if ext is not None:
        cmd.insert(1, '--gzip-compressed')

    call(cmd, cwd=dir_temp)


def run_kraken_pe(kraken, db, in_file_1, in_file_2, out_class_file,
                  out_unclass_file, report_file, confidence, threads,
                  dir_temp, linfo=print):  # noqa

    cmd = [kraken, '--db', db, '--threads', str(threads), '--paired',
           '--confidence', str(confidence), '--output', '-',
           '--report', report_file, '--use-names',
           '--classified-out', out_class_file,
           '--unclassified-out', out_unclass_file, in_file_1, in_file_2]

    mm = _use_memory_mapping(db, linfo=linfo)
    if mm is not None:
        cmd.insert(1, mm)

    _, _, ext = splitext_gz(in_file_1)
    if ext is not None:
        cmd.insert(1, '--gzip-compressed')

    call(cmd, cwd=dir_temp)


def run_kraken_filters(order, dbs, base_name, in_files, dir_out, confidence,
                       kraken2, threads, dir_temp, linfo=print):  # noqa

    dbs_ordered = OrderedDict()
    for dbn in order:
        db_name = dbn[0]
        if db_name in dbs:
            dbs_ordered[db_name] = dbs[db_name]

    # SE
    if isinstance(in_files, (str, bytes)):
        in_file = in_files
        _, in_file_ext, _ = splitext_gz(in_file)
        for i, db in enumerate(dbs_ordered):
            linfo('Filtering SE reads using Kraken2 database: ' + db)
            dir_out_db = opj(dir_out, db)
            make_dir(dir_out_db)
            report_file = opj(dir_out_db, base_name + '.txt')
            out_class_file = opj(dir_out_db, base_name + in_file_ext)
            out_unclass_file = opj(dir_temp, 'zzztemp_' + db + in_file_ext)

            run_kraken_se(
                kraken=kraken2,
                db=dbs_ordered[db],
                in_file=in_file,
                out_class_file=out_class_file,
                out_unclass_file=out_unclass_file,
                report_file=report_file,
                confidence=confidence,
                threads=threads,
                dir_temp=dir_temp,
                linfo=linfo)

            if i > 0:
                remove(in_file)
            in_file = out_unclass_file

        move(in_file, opj(dir_out, base_name + in_file_ext))

    # PE
    elif isinstance(in_files, (list, tuple, set)):

        assert len(in_files) > 1

        _, in_file_ext, _ = splitext_gz(in_files[0])

        in_file_R1 = in_files[0]
        in_file_R2 = in_files[1]

        if len(in_files) > 2:
            in_file = in_files[2]

            for i, db in enumerate(dbs_ordered):
                linfo('Filtering unpaired forward reads using Kraken2 database: ' + db)
                dir_out_db = opj(dir_out, db)
                make_dir(dir_out_db)
                report_file = opj(dir_out_db, base_name + '_unpaired_1.txt')
                out_class_file = opj(dir_out_db, base_name + '_unpaired_1' + in_file_ext)
                out_unclass_file = opj(dir_temp, 'zzztemp_' + db + in_file_ext)

                run_kraken_se(
                    kraken=kraken2,
                    db=dbs_ordered[db],
                    in_file=in_file,
                    out_class_file=out_class_file,
                    out_unclass_file=out_unclass_file,
                    report_file=report_file,
                    confidence=confidence,
                    threads=threads,
                    dir_temp=dir_temp,
                    linfo=linfo)

                if i > 0:
                    remove(in_file)
                in_file = out_unclass_file

            move(in_file, opj(dir_out, base_name + '_unpaired_1' + in_file_ext))

        if len(in_files) == 4:
            in_file = in_files[3]

            for i, db in enumerate(dbs_ordered):
                linfo('Filtering unpaired reverse reads using Kraken2 database: ' + db)
                dir_out_db = opj(dir_out, db)
                make_dir(dir_out_db)
                report_file = opj(dir_out_db, base_name + '_unpaired_2.txt')
                out_class_file = opj(dir_out_db, base_name + '_unpaired_2' + in_file_ext)
                out_unclass_file = opj(dir_temp, 'zzztemp_' + db + in_file_ext)

                run_kraken_se(
                    kraken=kraken2,
                    db=dbs_ordered[db],
                    in_file=in_file,
                    out_class_file=out_class_file,
                    out_unclass_file=out_unclass_file,
                    report_file=report_file,
                    confidence=confidence,
                    threads=threads,
                    dir_temp=dir_temp,
                    linfo=linfo)

                if i > 0:
                    remove(in_file)
                in_file = out_unclass_file

            move(in_file, opj(dir_out, base_name + '_unpaired_2' +
                              in_file_ext))

        for i, db in enumerate(dbs_ordered):
            linfo('Filtering paired reads using Kraken2 database: ' + db)
            dir_out_db = opj(dir_out, db)
            make_dir(dir_out_db)
            report_file = opj(dir_out_db, base_name + '_paired.txt')
            out_class_file = opj(dir_out_db, base_name + '_paired#' + in_file_ext)
            out_unclass_file = opj(dir_temp, 'zzztemp_' + db + '_paired#' + in_file_ext)

            run_kraken_pe(
                kraken=kraken2,
                db=dbs_ordered[db],
                in_file_1=in_file_R1,
                in_file_2=in_file_R2,
                out_class_file=out_class_file,
                out_unclass_file=out_unclass_file,
                report_file=report_file,
                confidence=confidence,
                threads=threads,
                dir_temp=dir_temp,
                linfo=linfo)

            if i > 0:
                remove(in_file_R1)
                remove(in_file_R2)

            in_file_R1 = out_unclass_file.replace('#', '_1')
            in_file_R2 = out_unclass_file.replace('#', '_2')

        move(in_file_R1, opj(dir_out, base_name + '_paired_1' + in_file_ext))
        move(in_file_R2, opj(dir_out, base_name + '_paired_2' + in_file_ext))
