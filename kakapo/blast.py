# -*- coding: utf-8 -*-

"""BLAST"""

import csv
from io import StringIO
import os.path

from operator import itemgetter
from collections import OrderedDict
from collections import Counter

from kakapo.shell import call
from kakapo.py_v_diffs import unicode

BLST_RES_COLS_1 = ['sseqid']
BLST_RES_COLS_2 = ['sseqid', 'evalue', 'sframe', 'sstart', 'send', 'qseqid']


def make_blast_db(exec_file, in_file, out_file, title, dbtype='nucl'):
    """Wrap makeblastdb"""
    cmd = [exec_file,
           '-in', in_file,
           '-out', out_file,
           '-title', title,
           '-dbtype', dbtype]

    call(cmd)


def run_blast(exec_file, task, threads, db_path, queries_file, out_file,
              evalue, qcov_hsp_perc, culling_limit, max_target_seqs,
              db_genetic_code, out_cols=BLST_RES_COLS_1):
    """Wrap blastn and tblastn"""
    exec_name = os.path.basename(exec_file)
    if exec_name in ['tblastn', ]:
        db_genetic_code = ['-db_gencode', db_genetic_code]
    else:
        db_genetic_code = []

    cmd = [exec_file,
           '-task', task,
           '-num_threads', str(threads),
           '-db', db_path,
           '-query', queries_file,
           '-out', out_file,
           '-evalue', evalue,
           '-outfmt', '6 delim=\t ' + ' '.join(out_cols),
           '-qcov_hsp_perc', qcov_hsp_perc,
           '-culling_limit', culling_limit,
           '-max_target_seqs', max_target_seqs]

    cmd = cmd + db_genetic_code
    # print(' '.join(cmd))
    call(cmd)


def _prs_blst_res(blast_results, col_names):  # noqa

    rdr = csv.DictReader(StringIO(blast_results),
                         fieldnames=col_names, delimiter='\t')

    lst = [row for row in rdr]
    return lst


def _prs_blst_res_f(blast_results_file, col_names):  # noqa
    with open(blast_results_file, 'r') as f:
        blast_results = unicode(f.read())
    parsed = _prs_blst_res(blast_results, col_names)
    return parsed


def parse_blast_results_file(blast_results_file, col_names=BLST_RES_COLS_2):  # noqa
    parsed_raw = _prs_blst_res_f(blast_results_file, col_names)
    parsed = dict()

    for rec in parsed_raw:
        subject = rec['sseqid']
        frame = int(rec['sframe'])
        start = int(rec['sstart'])
        end = int(rec['send'])

        if frame > 0:
            start = start - 1
        else:
            end = end - 1

        if subject not in parsed:
            parsed[subject] = []

        rec_dict = OrderedDict()

        rec_dict['frame'] = frame
        rec_dict['start'] = start
        rec_dict['end'] = end

        parsed[subject].append(rec_dict)

    return parsed


def collate_blast_results(parsed_blast_results):  # noqa

    coll = []

    for target in parsed_blast_results:

        hits = parsed_blast_results[target]

        t_frames = []
        t_starts = []
        t_ends = []

        for hit in hits:
            hit_frame = hit['frame']
            hit_start = hit['start']
            hit_end = hit['end']

            t_frames.append(hit_frame)
            t_starts.append(hit_start)
            t_ends.append(hit_end)

        t_frame_counts = Counter(t_frames)
        t_chosen_frame = t_frame_counts.most_common(1)[0][0]
        idx = [i for i, x in enumerate(t_frames) if x == t_chosen_frame]

        t_frames = [t_frames[i] for i in idx]
        t_starts = [t_starts[i] for i in idx]
        t_ends = [t_ends[i] for i in idx]

        t_frame = t_frames[0]
        t_start = None
        t_end = None

        if t_frame > 0:
            t_start = min(t_starts)
            t_end = max(t_ends)
        else:
            t_start = max(t_starts)
            t_end = min(t_ends)

        t_coll = OrderedDict()
        t_coll['sseqid'] = target
        t_coll['frame'] = t_frame
        t_coll['start'] = t_start
        t_coll['end'] = t_end

        coll.append(t_coll)

    coll_sorted = sorted(coll, key=itemgetter('sseqid'), reverse=False)

    return coll_sorted


# tblastn -help
# '-outfmt'
# Options 6, 7 and 10 can be additionally configured to produce
#    a custom format specified by space delimited format specifiers,
#    or by a token specified by the delim keyword.
#     E.g.: "17 delim=@ qacc sacc score".
#    The delim keyword must appear after the numeric output format
#    specification.
#    The supported format specifiers are:
#         qseqid means Query Seq-id
#            qgi means Query GI
#           qacc means Query accesion
#        qaccver means Query accesion.version
#           qlen means Query sequence length
#         sseqid means Subject Seq-id
#      sallseqid means All subject Seq-id(s), separated by a ';'
#            sgi means Subject GI
#         sallgi means All subject GIs
#           sacc means Subject accession
#        saccver means Subject accession.version
#        sallacc means All subject accessions
#           slen means Subject sequence length
#         qstart means Start of alignment in query
#           qend means End of alignment in query
#         sstart means Start of alignment in subject
#           send means End of alignment in subject
#           qseq means Aligned part of query sequence
#           sseq means Aligned part of subject sequence
#         evalue means Expect value
#       bitscore means Bit score
#          score means Raw score
#         length means Alignment length
#         pident means Percentage of identical matches
#         nident means Number of identical matches
#       mismatch means Number of mismatches
#       positive means Number of positive-scoring matches
#        gapopen means Number of gap openings
#           gaps means Total number of gaps
#           ppos means Percentage of positive-scoring matches
#         frames means Query and subject frames separated by a '/'
#         qframe means Query frame
#         sframe means Subject frame
#           btop means Blast traceback operations (BTOP)
#         staxid means Subject Taxonomy ID
#       ssciname means Subject Scientific Name
#       scomname means Subject Common Name
#     sblastname means Subject Blast Name
#      sskingdom means Subject Super Kingdom
#        staxids means unique Subject Taxonomy ID(s), separated by a ';'
#              (in numerical order)
#      sscinames means unique Subject Scientific Name(s), separated by a ';'
#      scomnames means unique Subject Common Name(s), separated by a ';'
#     sblastnames means unique Subject Blast Name(s), separated by a ';'
#              (in alphabetical order)
#     sskingdoms means unique Subject Super Kingdom(s), separated by a ';'
#              (in alphabetical order)
#         stitle means Subject Title
#     salltitles means All Subject Title(s), separated by a '<>'
#        sstrand means Subject Strand
#          qcovs means Query Coverage Per Subject
#        qcovhsp means Query Coverage Per HSP
#         qcovus means Query Coverage Per Unique Subject (blastn only)
#    When not provided, the default value is:
#    'qaccver saccver pident length mismatch gapopen qstart qend sstart send
#    evalue bitscore', which is equivalent to the keyword 'std'
#    Default = `0'
