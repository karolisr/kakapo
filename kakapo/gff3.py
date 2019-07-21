# -*- coding: utf-8 -*-

"""GFF3"""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import with_statement

import json

def gff_from_kakapo_ips5_json(json_dict):  # noqa

    gff_text = '##gff-version 3\n'

    for rec_key in json_dict:

        rec = json_dict[rec_key][0]
        ann = rec['kakapo_annotations']

        orf_begin = ann['orf_begin']
        orf_end = ann['orf_end']
        # orf_frame = ann['orf_frame']
        blast_hit_begin = ann['blast_hit_begin']
        blast_hit_end = ann['blast_hit_end']

        gff_seqid = rec_key
        gff_source = '.'
        gff_type = 'SO:0000236'
        gff_start = str(orf_begin + 1)
        gff_end = str(orf_end)
        gff_score = '.'
        gff_strand = '+'
        gff_phase = str(0)
        gff_attributes = ''

        gff_cols = [gff_seqid, gff_source, gff_type, gff_start, gff_end,
                    gff_score, gff_strand, gff_phase, gff_attributes]

        gff_row = '\t'.join(gff_cols)
        gff_text = gff_text + gff_row + '\n'

        gff_type = 'BLAST Hit'
        gff_start = str(blast_hit_begin + 1)
        gff_end = str(blast_hit_end)

        gff_cols = [gff_seqid, gff_source, gff_type, gff_start, gff_end,
                    gff_score, gff_strand, gff_phase, gff_attributes]

        gff_row = '\t'.join(gff_cols)
        gff_text = gff_text + gff_row + '\n'

    return gff_text


def gff_from_kakapo_ips5_json_file(json_path, gff_path=None):  # noqa

    with open(json_path, 'r') as f:
        json_dict = json.load(f)

    gff_text = gff_from_kakapo_ips5_json(json_dict)

    if gff_path is not None:
        with open(gff_path, 'w') as f:
            f.write(gff_text)

    return gff_text


if __name__ == '__main__':

    from os.path import join as opj
    from os.path import basename
    from os.path import splitext

    base_path = opj('tests', 'data', 'inter_pro_scan_result_samples')
    json_path = opj(base_path, 'sample_02.json')
    raw_fa_path = opj(base_path, 'sample_02_transcripts_nt.fasta')
    gff_path = opj(base_path, splitext(basename(raw_fa_path))[0] + '.gff')
    gff_from_kakapo_ips5_json_file(json_path, gff_path)
