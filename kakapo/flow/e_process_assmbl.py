# -*- coding: utf-8 -*-

"""Kakapo workflow: Assemble Reads."""

from operator import itemgetter
from os import remove as osremove
from os import stat as osstat
from os.path import basename
from os.path import exists as ope
from os.path import join as opj
from os.path import splitext
from sys import exit

from kakapo.bioio import read_fasta
from kakapo.blast import make_blast_db
from kakapo.helpers import combine_text_files
from kakapo.helpers import make_dir
from kakapo.spades import run_spades_se, run_spades_pe
from kakapo.translation_tables import TranslationTable


def run_spades(se_fastq_files, pe_fastq_files, dir_spades_assemblies,
               spades, dir_temp, ss, threads, ram, linfo=print):

    if len(se_fastq_files) > 0 or len(pe_fastq_files) > 0:
        if spades is None:
            linfo('spades is not available. ' +
                  'Cannot continue. Exiting.')
            exit(0)

    for se in se_fastq_files:
        dir_results = opj(dir_spades_assemblies, se + '__' + ss)
        fq_path = se_fastq_files[se]['vsearch_results_path' + '__' + ss]
        se_fastq_files[se]['spades_assembly' + '__' + ss] = None

        if ope(dir_results):
            linfo('SPAdes results for sample ' + se + ' already exist [' + ss + ']')
        else:
            make_dir(dir_results)
            linfo('Running SPAdes on: ' + se + ' [' + ss + ']')
            run_spades_se(spades,
                          out_dir=dir_results,
                          input_file=fq_path,
                          threads=threads,
                          memory=ram,
                          rna=True)

        assmbl_path = opj(dir_results, 'transcripts.fasta')
        if ope(assmbl_path):
            count = len(read_fasta(assmbl_path))
            tr_str = ' transcripts'
            if count == 1:
                tr_str = ' transcript'
            linfo('SPAdes produced ' + str(count) + tr_str + ' [' + ss + ']')
            se_fastq_files[se]['spades_assembly' + '__' + ss] = assmbl_path
        else:
            linfo('SPAdes produced no transcripts [' + ss + ']')

    for pe in pe_fastq_files:
        dir_results = opj(dir_spades_assemblies, pe + '__' + ss)
        fq_paths = pe_fastq_files[pe]['vsearch_results_path' + '__' + ss]
        pe_fastq_files[pe]['spades_assembly' + '__' + ss] = None

        if ope(dir_results):
            linfo('SPAdes results for sample ' + pe + ' already exist [' + ss + ']')
        else:
            make_dir(dir_results)
            linfo('Running SPAdes on: ' + pe + ' [' + ss + ']')

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
            count = len(read_fasta(assmbl_path))
            tr_str = ' transcripts'
            if count == 1:
                tr_str = ' transcript'
            linfo('SPAdes produced ' + str(count) + tr_str + ' [' + ss + ']')
            pe_fastq_files[pe]['spades_assembly' + '__' + ss] = assmbl_path
        else:
            linfo('SPAdes produced no transcripts [' + ss + ']')


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
        assemblies.append(a)

    assemblies = sorted(assemblies, key=itemgetter('name'), reverse=False)

    return assemblies


def makeblastdb_assemblies(assemblies, dir_prj_blast_assmbl, makeblastdb,
                           linfo=print):
    if len(assemblies) > 0:
        linfo('Building BLAST databases for assemblies')
        if makeblastdb is None:
            linfo('makeblastdb is not available. ' +
                  'Cannot continue. Exiting.')
            exit(0)
    for a in assemblies:
        assmbl_name = a['name']

        assmbl_blast_db_dir = opj(dir_prj_blast_assmbl, assmbl_name)
        assmbl_blast_db_file = opj(assmbl_blast_db_dir, assmbl_name)

        a['blast_db_path'] = assmbl_blast_db_file

        if ope(assmbl_blast_db_dir):
            linfo('BLAST database for ' + assmbl_name + ' already exists')
        else:
            linfo(assmbl_name)
            make_dir(assmbl_blast_db_dir)
            make_blast_db(exec_file=makeblastdb,
                          in_file=a['path'],
                          out_file=assmbl_blast_db_file,
                          title=assmbl_name)
