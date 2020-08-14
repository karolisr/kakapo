"""vsearch."""

from kakapo.utils.subp import run


def run_cluster_fast(vsearch, ident, in_file, out_file):

    cmd = [vsearch,
           '--cluster_fast', in_file,
           '--centroids', out_file,
           '--fasta_width', '0',
           '--id', str(ident)]

    run(cmd, do_not_raise=True)


def run_vsearch(vsearch, ident, q_file, db_file, out_file, minlen):

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
           '--id', str(ident)]

    run(cmd, do_not_raise=True)
