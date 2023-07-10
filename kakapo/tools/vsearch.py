"""vsearch."""

from kakapo.utils.subp import run


def run_cluster_fast(vsearch, ident, in_file, out_file_centroids,
                     out_file_prefix_clusters=None,
                     clusterout_id=False,
                     iddef=None):

    cmd = [vsearch,
           '--cluster_fast', in_file,
           '--centroids', out_file_centroids,
           '--fasta_width', '0',
           '--id', str(ident)]

    if clusterout_id is True:
        cmd.append('--clusterout_id')

    if iddef is not None:
        cmd += ['--iddef', str(iddef)]

    if out_file_prefix_clusters is not None:
        cmd += ['--clusters', out_file_prefix_clusters]

    run(cmd, do_not_raise=False)


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
