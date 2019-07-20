# -*- coding: utf-8 -*-

"""Trimmomatic"""

from kakapo.shell import call


def cluster_fast(vsearch, in_file, out_file):  # noqa

    cmd = [vsearch,
           '--cluster_fast', in_file,
           '--centroids', out_file,
           '--fasta_width', '0',
           '--id', '0.9']

    call(cmd)

def vsearch(vsearch, q_file, db_file, out_file, minlen):  # noqa

    cmd = [vsearch,
           '--usearch_global', q_file,
           '--db', db_file,
           '--userout', out_file,
           '--userfields', 'target',
           '--maxseqlength', '1000',
           '--minseqlength', str(minlen),
           '--threads', '0',
           '--strand', 'both',
           '--maxaccepts', '0',
           '--maxrejects', '0',
           '--iddef', '2',
           '--maxsubs', '3',
           '--maxgaps', '1',
           '--target_cov', '0.33',
           '--id', '0.85']

    call(cmd)
