# -*- coding: utf-8 -*-

"""GFF3"""

# https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
# https://useast.ensembl.org/info/website/upload/gff3.html
# http://gmod.org/wiki/GFF3

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import with_statement

import json
from collections import OrderedDict
# from kakapo.helpers import debug_print as pp


def gff_template():  # noqa
    t = OrderedDict()

    t['seqid'] = '.'
    t['source'] = '.'
    t['type'] = '.'
    t['start'] = '.'
    t['end'] = '.'
    t['score'] = '.'
    t['strand'] = '.'
    t['phase'] = '.'
    t['attributes'] = '.'

    return t


def gff_text(gff_dict):  # noqa
    return '\t'.join(list(gff_dict.values())) + '\n'


def gff_from_kakapo_json(json_dict):  # noqa
    gff = ''

    for rec_key in json_dict:

        rec = json_dict[rec_key][0]
        ann = rec['kakapo_annotations']

        # ORF annotation -----------------------------------------------------
        # orf_frame = ann['orf_frame']
        orf_begin = ann['orf_begin']
        orf_end = ann['orf_end']
        entry_orf = gff_template()
        entry_orf['seqid'] = rec_key
        entry_orf['source'] = 'kakapo'
        entry_orf['type'] = 'ORF'
        entry_orf['start'] = str(orf_begin + 1)
        entry_orf['end'] = str(orf_end)
        entry_orf['score'] = '.'
        entry_orf['strand'] = '+'
        entry_orf['phase'] = str(0)
        entry_orf['attributes'] = 'name=Kakapo Predicted ORF'
        gff = gff + gff_text(entry_orf)

        # BLAST Hit annotation -----------------------------------------------
        entry_blast_hit = entry_orf
        blast_hit_begin = ann['blast_hit_begin']
        blast_hit_end = ann['blast_hit_end']
        entry_blast_hit['type'] = 'BLAST Hit'
        entry_blast_hit['start'] = str(blast_hit_begin + 1)
        entry_blast_hit['end'] = str(blast_hit_end)
        entry_blast_hit['attributes'] = 'name=Kakapo Merged BLAST Hits'
        gff = gff + gff_text(entry_orf)
        # --------------------------------------------------------------------

    return gff


def gff_from_kakapo_ips5_json(json_dict):  # noqa
    gff = ''

    for rec_key in json_dict:
        rec = json_dict[rec_key][0]   # matches md5 sequence xref
        ann = rec['kakapo_annotations']
        orf_begin = ann['orf_begin']
        rec_matches = rec['matches']  # list of dicts
        for m in rec_matches:
            locations = m['locations'][0]
            beg = str(locations['start'] * 3 - 2 + orf_begin)
            end = str(locations['end'] * 3 + orf_begin)
            sig = m['signature']
            lib_info = sig['signatureLibraryRelease']
            accession = sig['accession']
            description = sig['description']
            lib_name = lib_info['library']
            # lib_vers = lib_info['version']

            attributes = ''

            if lib_name == 'PFAM':
                attributes = ('name=' + description + ';' +
                              'accession=' + accession + ';')
                score = str(locations['score'])
                g = gff_template()
                g['seqid'] = rec_key
                g['source'] = lib_name
                g['type'] = 'PROTEIN_FAMILY'
                g['start'] = beg
                g['end'] = end
                g['score'] = score
                g['strand'] = '+'
                g['phase'] = str(0)
                g['attributes'] = attributes

                gff = gff + gff_text(g)

    return gff


def gff_from_kakapo_ips5_json_file(json_path, gff_path=None):  # noqa

    with open(json_path, 'r') as f:
        json_dict = json.load(f)

    # gff = '##gff-version 3\n'

    gff = (gff_from_kakapo_json(json_dict) +
           gff_from_kakapo_ips5_json(json_dict))

    if gff_path is not None:
        with open(gff_path, 'w') as f:
            f.write(gff)

    return gff


# if __name__ == '__main__':

#     from os.path import join as opj
#     from os.path import basename
#     from os.path import splitext

#     print()

#     base_path = opj('tests', 'data', 'inter_pro_scan_result_samples')
#     json_path = opj(base_path, 'sample_02.json')
#     raw_fa_path = opj(base_path, 'sample_02_transcripts_nt.fasta')
#     gff_path = opj(base_path, splitext(basename(raw_fa_path))[0] + '.gff')
#     gff = gff_from_kakapo_ips5_json_file(json_path, gff_path)

#     print(gff)


# seqid - name of the chromosome or scaffold; chromosome names can be given
#     with or without the 'chr' prefix. Important note: the seq ID must be one
#     used within Ensembl, i.e. a standard chromosome name or an Ensembl
#     identifier such as a scaffold ID, without any additional content such as
#     species or assembly. See the example GFF output below.

# source - name of the program that generated this feature, or the data source
#     (database or project name)

# type - type of feature. Must be a term or accession from the SOFA sequence
#     ontology

# start - Start position of the feature, with sequence numbering starting at 1.

# end - End position of the feature, with sequence numbering starting at 1.

# score - A floating point value.

# strand - defined as + (forward) or - (reverse).

# phase - One of '0', '1' or '2'. '0' indicates that the first base of the
#     feature is the first base of a codon, '1' that the second base is the
#     first base of a codon, and so on..

# attributes - A semicolon-separated list of tag-value pairs, providing
#     additional information about each feature. Some of these tags are
#     predefined, e.g. ID, Name, Alias, Parent - see the GFF documentation
#     for more details.
