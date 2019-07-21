# -*- coding: utf-8 -*-

"""GFF3"""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import with_statement

import json
import re

def gff_from_kakapo_ips5_json(json_dict):  # noqa

    orf_regex = r'.*__ORF:(\d+):(\d+):(\d)'

    gff_text = '##gff-version 3\n'

    for rec_key in json_dict:
        orf_info_raw = re.match(orf_regex, rec_key).groups()
        orf_info = [int(x) for x in orf_info_raw]

        beg = orf_info[0]
        end = orf_info[1]
        # frm = orf_info[2]

        gff_seqid = rec_key
        gff_source = '.'
        gff_type = 'SO:0000236'
        gff_start = str(beg + 1)
        gff_end = str(end)
        gff_score = '.'
        gff_strand = '+'
        gff_phase = str(0)
        gff_attributes = 'ID=kakapo predicted'

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
