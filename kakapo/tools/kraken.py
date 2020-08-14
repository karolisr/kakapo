"""Kraken2."""

from collections import OrderedDict
from os import remove, stat
from os.path import join as opj
from os.path import splitext, basename
from shutil import copyfile
from shutil import move

from kakapo.tools.config import RAM
from kakapo.utils.logging import Log
from kakapo.utils.misc import splitext_gz, make_dirs
from kakapo.utils.subp import run


def _use_memory_mapping(db_path):
    db_size = stat(opj(db_path, 'hash.k2d')).st_size / (1024 ** 3)
    mem_max = RAM / 3
    if mem_max < db_size:
        db_name = splitext(basename(db_path))[0]
        Log.err('Not enough memory for Kraken2 database {}. '
                'Switching to a slower memory-mapping mode.'.format(db_name))
        return '--memory-mapping'
    else:
        return None


def run_kraken_se(kraken, db, in_file, out_class_file, out_unclass_file,
                  report_file, confidence, threads, dir_temp):

    cmd = [kraken, '--db', db, '--threads', str(threads),
           '--confidence', str(confidence), '--output', '-',
           '--report', report_file, '--use-names',
           '--classified-out', out_class_file,
           '--unclassified-out', out_unclass_file, in_file]

    mm = _use_memory_mapping(db)
    if mm is not None:
        cmd.insert(1, mm)

    _, _, ext = splitext_gz(in_file)
    if ext is not None:
        cmd.insert(1, '--gzip-compressed')

    run(cmd, cwd=dir_temp, do_not_raise=True)


def run_kraken_pe(kraken, db, in_file_1, in_file_2, out_class_file,
                  out_unclass_file, report_file, confidence, threads,
                  dir_temp):

    cmd = [kraken, '--db', db, '--threads', str(threads), '--paired',
           '--confidence', str(confidence), '--output', '-',
           '--report', report_file, '--use-names',
           '--classified-out', out_class_file,
           '--unclassified-out', out_unclass_file, in_file_1, in_file_2]

    mm = _use_memory_mapping(db)
    if mm is not None:
        cmd.insert(1, mm)

    _, _, ext = splitext_gz(in_file_1)
    if ext is not None:
        cmd.insert(1, '--gzip-compressed')

    run(cmd, cwd=dir_temp, do_not_raise=True)


def run_kraken_filters(order, dbs, base_name, in_files, dir_out, confidence,
                       kraken2, threads, dir_temp):

    dbs_ordered = OrderedDict()
    for dbn in order:
        db_name = dbn[0]
        if db_name in dbs:
            dbs_ordered[db_name] = dbs[db_name]
        else:
            Log.wrn('Kraken2 database not found:', db_name)

    # SE
    if isinstance(in_files, (str, bytes)):
        in_file = in_files
        _, in_file_ext, _ = splitext_gz(in_file)
        for i, db in enumerate(dbs_ordered):
            Log.msg('Filtering SE reads using Kraken2 database:', db)
            dir_out_db = opj(dir_out, db)
            make_dirs(dir_out_db)
            report_file = opj(dir_out_db, base_name + '.txt')
            out_class_file = opj(dir_out_db, base_name + in_file_ext)
            out_unclass_file = opj(dir_temp, base_name + '_' + db + '_kraken2_unclassified' + in_file_ext)

            if stat(in_file).st_size > 0:  # Kraken2 freaks out if the file is empty.
                run_kraken_se(
                    kraken=kraken2,
                    db=dbs_ordered[db],
                    in_file=in_file,
                    out_class_file=out_class_file,
                    out_unclass_file=out_unclass_file,
                    report_file=report_file,
                    confidence=confidence,
                    threads=threads,
                    dir_temp=dir_temp)
            else:
                copyfile(in_file, out_class_file)
                copyfile(in_file, out_unclass_file)

            if i > 0:
                remove(in_file)
            in_file = out_unclass_file

        move(in_file, opj(dir_out, base_name + in_file_ext))

    # PE
    elif isinstance(in_files, (list, tuple)):

        assert len(in_files) > 1

        _, in_file_ext, _ = splitext_gz(in_files[0])

        in_file_R1 = in_files[0]
        in_file_R2 = in_files[1]

        if len(in_files) > 2:
            in_file = in_files[2]

            for i, db in enumerate(dbs_ordered):
                Log.msg('Filtering unpaired forward reads using Kraken2 database:', db)
                dir_out_db = opj(dir_out, db)
                make_dirs(dir_out_db)
                report_file = opj(dir_out_db, base_name + '_unpaired_1.txt')
                out_class_file = opj(dir_out_db, base_name + '_unpaired_1' + in_file_ext)
                out_unclass_file = opj(dir_temp, base_name + '_' + db + '_kraken2_unclassified' + in_file_ext)

                if stat(in_file).st_size > 0:
                    run_kraken_se(
                        kraken=kraken2,
                        db=dbs_ordered[db],
                        in_file=in_file,
                        out_class_file=out_class_file,
                        out_unclass_file=out_unclass_file,
                        report_file=report_file,
                        confidence=confidence,
                        threads=threads,
                        dir_temp=dir_temp)
                else:
                    copyfile(in_file, out_class_file)
                    copyfile(in_file, out_unclass_file)

                if i > 0:
                    remove(in_file)
                in_file = out_unclass_file

            move(in_file, opj(dir_out, base_name + '_unpaired_1' + in_file_ext))

        if len(in_files) == 4:
            in_file = in_files[3]

            for i, db in enumerate(dbs_ordered):
                Log.msg('Filtering unpaired reverse reads using Kraken2 database:', db)
                dir_out_db = opj(dir_out, db)
                make_dirs(dir_out_db)
                report_file = opj(dir_out_db, base_name + '_unpaired_2.txt')
                out_class_file = opj(dir_out_db, base_name + '_unpaired_2' + in_file_ext)
                out_unclass_file = opj(dir_temp, base_name + '_' + db + '_kraken2_unclassified' + in_file_ext)

                if stat(in_file).st_size > 0:
                    run_kraken_se(
                        kraken=kraken2,
                        db=dbs_ordered[db],
                        in_file=in_file,
                        out_class_file=out_class_file,
                        out_unclass_file=out_unclass_file,
                        report_file=report_file,
                        confidence=confidence,
                        threads=threads,
                        dir_temp=dir_temp)
                else:
                    copyfile(in_file, out_class_file)
                    copyfile(in_file, out_unclass_file)

                if i > 0:
                    remove(in_file)
                in_file = out_unclass_file

            move(in_file, opj(dir_out, base_name + '_unpaired_2' +
                              in_file_ext))

        for i, db in enumerate(dbs_ordered):
            Log.msg('Filtering paired reads using Kraken2 database:', db)
            dir_out_db = opj(dir_out, db)
            make_dirs(dir_out_db)
            report_file = opj(dir_out_db, base_name + '_paired.txt')
            out_class_file = opj(dir_out_db, base_name + '_paired#' + in_file_ext)
            out_unclass_file = opj(dir_temp, base_name + '_' + db + '_kraken2_unclassified' + '_paired#' + in_file_ext)

            if stat(in_file_R1).st_size > 0 and stat(in_file_R2).st_size > 0:
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
                    dir_temp=dir_temp)
            else:
                copyfile(in_file_R1, copyfile(in_file_R1, out_class_file.replace('#', '_1')))
                copyfile(in_file_R2, copyfile(in_file_R2, out_class_file.replace('#', '_2')))
                copyfile(in_file_R1, copyfile(in_file_R1, out_unclass_file.replace('#', '_1')))
                copyfile(in_file_R2, copyfile(in_file_R2, out_unclass_file.replace('#', '_2')))

            if i > 0:
                remove(in_file_R1)
                remove(in_file_R2)

            in_file_R1 = out_unclass_file.replace('#', '_1')
            in_file_R2 = out_unclass_file.replace('#', '_2')

        move(in_file_R1, opj(dir_out, base_name + '_paired_1' + in_file_ext))
        move(in_file_R2, opj(dir_out, base_name + '_paired_2' + in_file_ext))
