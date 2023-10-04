"""Kakapo workflow: Assemble Reads."""

from operator import itemgetter
from os import remove as osremove
from os import stat as osstat
from os.path import basename
from os.path import exists as ope
from os.path import join as opj
from os.path import splitext
from sys import exit

from kakapo.tools.bioio import read_fasta
from kakapo.tools.blast import make_blast_db
from kakapo.tools.seq import SEQ_TYPE_NT
from kakapo.tools.spades import run_spades_pe, run_spades_se
from kakapo.tools.transl_tables import TranslationTable
from kakapo.utils.logging import Log
from kakapo.utils.misc import combine_text_files, make_dirs


def run_spades(se_fastq_files, pe_fastq_files, dir_spades_assemblies,
               spades, dir_temp, ss, threads, ram):

    if len(se_fastq_files) > 0 or len(pe_fastq_files) > 0:
        if spades is None:
            Log.err('SPAdes is not available. Cannot continue. Exiting.', '')
            exit(0)

    for se in se_fastq_files:
        dir_results = opj(dir_spades_assemblies, se + '__' + ss)
        fq_path = se_fastq_files[se]['vsearch_results_path' + '__' + ss]
        se_fastq_files[se]['spades_assembly' + '__' + ss] = None

        if fq_path is None:
            continue

        if ope(dir_results):
            Log.msg('SPAdes assembly already exists:', se)
        else:
            make_dirs(dir_results)
            Log.msg('Running SPAdes on:', se)
            run_spades_se(spades,
                          out_dir=dir_results,
                          input_file=fq_path,
                          threads=threads,
                          memory=ram,
                          rna=True)

        assmbl_path = opj(dir_results, 'transcripts.fasta')
        if ope(assmbl_path):
            count = len(read_fasta(assmbl_path, SEQ_TYPE_NT))
            tr_str = ' transcripts.'
            if count == 1:
                tr_str = ' transcript.'
            Log.log(msg='SPAdes produced ' + str(count) + tr_str, timestamp=False)
            se_fastq_files[se]['spades_assembly' + '__' + ss] = assmbl_path
        else:
            Log.log(wrn='SPAdes produced no transcripts.', timestamp=False)

    for pe in pe_fastq_files:
        dir_results = opj(dir_spades_assemblies, pe + '__' + ss)
        fq_paths = pe_fastq_files[pe]['vsearch_results_path' + '__' + ss]
        pe_fastq_files[pe]['spades_assembly' + '__' + ss] = None

        if fq_paths is None:
            continue

        if ope(dir_results):
            Log.msg('SPAdes assembly already exists:', pe)
        else:
            make_dirs(dir_results)
            Log.msg('Running SPAdes on: ' + pe, '')

            if osstat(fq_paths[0]).st_size > 0 and \
               osstat(fq_paths[1]).st_size > 0:

                run_spades_pe(spades,
                              out_dir=dir_results,
                              input_files=fq_paths,
                              threads=threads,
                              memory=ram,
                              rna=True)

            else:
                _ = opj(dir_temp, 'temp.fasta')
                combine_text_files(fq_paths, _)
                run_spades_se(spades,
                              out_dir=dir_results,
                              input_file=_,
                              threads=threads,
                              memory=ram,
                              rna=True)
                osremove(_)

        assmbl_path = opj(dir_results, 'transcripts.fasta')
        if ope(assmbl_path):
            count = len(read_fasta(assmbl_path, SEQ_TYPE_NT))
            tr_str = ' transcripts.'
            if count == 1:
                tr_str = ' transcript.'
            Log.log(msg='SPAdes produced ' + str(count) + tr_str, timestamp=False)
            pe_fastq_files[pe]['spades_assembly' + '__' + ss] = assmbl_path
        else:
            Log.log(wrn='SPAdes produced no transcripts.', timestamp=False)


def combine_assemblies(se_fastq_files, pe_fastq_files, user_assemblies, tax,
                       search_strategies):
    assemblies = []

    for ss in search_strategies:
        for se in se_fastq_files:
            a_path = se_fastq_files[se]['spades_assembly' + '__' + ss]
            if a_path is None:
                continue
            a = {}
            a['path'] = a_path
            a['name'] = se + '__' + ss
            a['src'] = se_fastq_files[se]['src']
            a['tax_id'] = se_fastq_files[se]['tax_id']
            a['gc_id'] = se_fastq_files[se]['gc_id']
            a['gc_tt'] = se_fastq_files[se]['gc_tt']
            a['gc_id_mito'] = se_fastq_files[se]['gc_id_mito']
            a['gc_tt_mito'] = se_fastq_files[se]['gc_tt_mito']
            a['gc_id_plastid'] = se_fastq_files[se]['gc_id_plastid']
            a['gc_tt_plastid'] = se_fastq_files[se]['gc_tt_plastid']

            assemblies.append(a)

        for pe in pe_fastq_files:
            a_path = pe_fastq_files[pe]['spades_assembly' + '__' + ss]
            if a_path is None:
                continue
            a = {}
            a['path'] = a_path
            a['name'] = pe + '__' + ss
            a['src'] = pe_fastq_files[pe]['src']
            a['tax_id'] = pe_fastq_files[pe]['tax_id']
            a['gc_id'] = pe_fastq_files[pe]['gc_id']
            a['gc_tt'] = pe_fastq_files[pe]['gc_tt']
            a['gc_id_mito'] = pe_fastq_files[pe]['gc_id_mito']
            a['gc_tt_mito'] = pe_fastq_files[pe]['gc_tt_mito']
            a['gc_id_plastid'] = pe_fastq_files[pe]['gc_id_plastid']
            a['gc_tt_plastid'] = pe_fastq_files[pe]['gc_tt_plastid']

            assemblies.append(a)

    for us in user_assemblies:
        a_path = us[1]
        gc = tax.genetic_code_for_taxid(us[0])
        a = {}
        a['path'] = a_path
        a['name'] = splitext(basename(a_path))[0]
        a['src'] = 'user_fasta'
        a['tax_id'] = us[0]
        a['gc_id'] = gc
        a['gc_tt'] = TranslationTable(gc)

        # ToDo: Refactor this with the code from __main__.py -----------------
        gc_mito = None
        tt_mito = None

        gc_plastid = None
        tt_plastid = None

        if tax.is_eukaryote(a['tax_id']) is True:
            gc_mito = tax.mito_genetic_code_for_taxid(a['tax_id'])
            if gc_mito != '0':
                tt_mito = TranslationTable(gc_mito)

            if tax.contains_plastid(a['tax_id']) is True:
                gc_plastid = tax.plastid_genetic_code_for_taxid(a['tax_id'])
                if gc_plastid != '0':
                    tt_plastid = TranslationTable(gc_plastid)

        a['gc_id_mito'] = gc_mito
        a['gc_tt_mito'] = tt_mito
        a['gc_id_plastid'] = gc_plastid
        a['gc_tt_plastid'] = tt_plastid
        # --------------------------------------------------------------------

        assemblies.append(a)

    assemblies = sorted(assemblies, key=itemgetter('name'), reverse=False)

    return assemblies


def makeblastdb_assemblies(assemblies, dir_prj_blast_assmbl, makeblastdb):
    if len(assemblies) > 0:
        print()
        Log.inf('Building BLAST databases for assemblies.')
        if makeblastdb is None:
            Log.err('makeblastdb is not available. Cannot continue. Exiting.', '')
            exit(0)
    for a in assemblies:
        assmbl_name = a['name']

        assmbl_blast_db_dir = opj(dir_prj_blast_assmbl, assmbl_name)
        assmbl_blast_db_file = opj(assmbl_blast_db_dir, assmbl_name)

        a['blast_db_path'] = assmbl_blast_db_file

        if ope(assmbl_blast_db_dir):
            Log.msg('BLAST database already exists:', assmbl_name)
        else:
            Log.msg(assmbl_name, '')
            make_dirs(assmbl_blast_db_dir)
            make_blast_db(exec_file=makeblastdb,
                          in_file=a['path'],
                          out_file=assmbl_blast_db_file,
                          title=assmbl_name)
