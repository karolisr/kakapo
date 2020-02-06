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
import re
from collections import OrderedDict


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


def gff_orf_blast_hit(ss, json_dict):  # noqa
    gff = ''

    for rec_key in json_dict:

        rec = json_dict[rec_key][0]
        ann = rec['kakapo_annotations__' + ss]

        # BLAST Hit annotation -----------------------------------------------
        if 'blast_hit_begin' not in ann or \
           'blast_hit_end' not in ann or \
           'frame' not in ann:
            continue

        frame = ann['frame']
        entry_blast_hit = gff_template()
        blast_hit_begin = ann['blast_hit_begin']
        blast_hit_end = ann['blast_hit_end']
        entry_blast_hit['seqid'] = rec_key
        entry_blast_hit['source'] = 'KAKAPO'
        entry_blast_hit['type'] = 'BLAST'
        entry_blast_hit['start'] = str(blast_hit_begin + 1)
        entry_blast_hit['end'] = str(blast_hit_end)
        entry_blast_hit['score'] = str(ann['evalue'])
        entry_blast_hit['strand'] = '+'
        entry_blast_hit['phase'] = str(0)
        query_name = ann['query_name'].replace('_', ' ')
        name = 'name=' + ann['query_name'] + '; Hit frame ' + str(frame)
        entry_blast_hit['attributes'] = name + ';note=Merged tblastn hits;'
        gff = gff + gff_text(entry_blast_hit)

        # ORF annotation -----------------------------------------------------
        if 'orf_begin' not in ann or \
           'orf_end' not in ann or \
           'frame' not in ann:
            continue

        orf_begin = ann['orf_begin']
        orf_end = ann['orf_end']
        entry_orf = entry_blast_hit
        entry_orf['type'] = 'ORF'
        entry_orf['start'] = str(orf_begin + 1)
        entry_orf['end'] = str(orf_end)
        entry_orf['score'] = str(ann['orf_grade'])
        name = 'name=ORF frame ' + str(frame)
        entry_orf['attributes'] = name
        gff = gff + gff_text(entry_orf)

        # --------------------------------------------------------------------

    return gff


def gff_orf_bad_blast_hit(ss, json_dict):  # noqa
    gff = ''

    for rec_key in json_dict:

        rec = json_dict[rec_key][0]
        ann = rec['kakapo_annotations__' + ss]

        # ORF bad annotation -------------------------------------------------
        if 'orfs_bad' not in ann:
            continue

        orfs_bad = ann['orfs_bad']

        for i, bad in enumerate(orfs_bad):

            orf_begin = bad['orf_begin']
            orf_end = bad['orf_end']
            orf_frame = bad['orf_frame']
            entry_orf = gff_template()
            entry_orf['seqid'] = rec_key
            entry_orf['source'] = 'KAKAPO'
            entry_orf['type'] = 'ORF BAD'
            entry_orf['start'] = str(orf_begin + 1)
            entry_orf['end'] = str(orf_end)
            entry_orf['score'] = str(bad['orf_grade'])
            entry_orf['strand'] = '+'
            entry_orf['phase'] = str(0)
            name = 'name=ORF BAD {} frame {}'.format(i, orf_frame)
            entry_orf['attributes'] = name
            gff = gff + gff_text(entry_orf)

        # --------------------------------------------------------------------

    return gff


def gff_pfam(ss, json_dict):  # noqa
    gff = ''

    for rec_key in json_dict:
        rec = json_dict[rec_key][0]
        ann = rec['kakapo_annotations__' + ss]

        if 'orf_begin' not in ann or \
           'matches' not in rec:
            continue

        orf_begin = ann['orf_begin']
        rec_matches = rec['matches']
        for m in rec_matches:
            locations = m['locations'][0]
            beg = str(locations['start'] * 3 - 2 + orf_begin)
            end = str(locations['end'] * 3 + orf_begin)
            sig = m['signature']
            lib_info = sig['signatureLibraryRelease']
            accession = sig['accession']
            description = sig['description']
            lib_name = lib_info['library']

            attributes = ''

            if lib_name == 'PFAM':

                fam_in_dec = re.findall(r'(^.*?)\s*family$', description, re.I)

                if len(fam_in_dec) != 0:
                    name = fam_in_dec[0]
                else:
                    name = description

                name = name.title()

                attributes = ('name=' + name + ': ' + accession + ';' +
                              'accession=' + accession + ';')

                ipr_acc = None
                entry = sig['entry']
                if entry is not None:
                    ipr_acc = entry['accession']
                if ipr_acc is not None:
                    attributes = attributes + 'interpro=' + ipr_acc + ';'

                score = str(locations['score'])
                g = gff_template()
                g['seqid'] = rec_key
                g['source'] = lib_name
                g['type'] = 'PfPfam'
                g['start'] = beg
                g['end'] = end
                g['score'] = score
                g['strand'] = '+'
                g['phase'] = str(0)
                g['attributes'] = attributes

                gff = gff + gff_text(g)

    return gff


def gff_phobius(ss, json_dict):  # noqa
    gff = ''

    phobius_translations = {'SIGNAL_PEPTIDE': 'SigPept',
                            'SIGNAL_PEPTIDE_N_REGION': 'SigPept1_N',
                            'SIGNAL_PEPTIDE_H_REGION': 'SigPept2_H',
                            'SIGNAL_PEPTIDE_C_REGION': 'SigPept3_C',
                            'NON_CYTOPLASMIC_DOMAIN': 'PhobNonCyt',
                            'TRANSMEMBRANE': 'PhobTrans',
                            'CYTOPLASMIC_DOMAIN': 'PhobCyt'}

    for rec_key in json_dict:
        rec = json_dict[rec_key][0]
        ann = rec['kakapo_annotations__' + ss]

        if 'orf_begin' not in ann or \
           'matches' not in rec:
            continue

        orf_begin = ann['orf_begin']
        rec_matches = rec['matches']
        for m in rec_matches:
            locations = m['locations'][0]
            beg = str(locations['start'] * 3 - 2 + orf_begin)
            end = str(locations['end'] * 3 + orf_begin)
            sig = m['signature']
            lib_info = sig['signatureLibraryRelease']
            accession = sig['accession']
            description = sig['description']
            lib_name = lib_info['library']

            if lib_name == 'PHOBIUS':
                name = sig['name']
                attributes = ('name=' + name + ';'
                              'description=' + description + ';')

                accession = phobius_translations[accession]

                g = gff_template()
                g['seqid'] = rec_key
                g['source'] = lib_name
                g['type'] = accession
                g['start'] = beg
                g['end'] = end
                g['score'] = '.'
                g['strand'] = '+'
                g['phase'] = str(0)
                g['attributes'] = attributes

                gff = gff + gff_text(g)

    return gff


def gff_panther(ss, json_dict):  # noqa
    gff = ''

    for rec_key in json_dict:
        rec = json_dict[rec_key][0]
        ann = rec['kakapo_annotations__' + ss]

        if 'orf_begin' not in ann or \
           'matches' not in rec:
            continue

        orf_begin = ann['orf_begin']
        rec_matches = rec['matches']

        for m in rec_matches:
            locations = m['locations'][0]
            beg = str(locations['start'] * 3 - 2 + orf_begin)
            end = str(locations['end'] * 3 + orf_begin)
            sig = m['signature']
            lib_info = sig['signatureLibraryRelease']
            accession = sig['accession']
            lib_name = lib_info['library']

            if lib_name == 'PANTHER':

                if ':' in accession:
                    continue

                name = sig['name'].title()
                if name == 'Family Not Named':
                    name = ''
                else:
                    name = ': ' + name
                attributes = ('name=' + accession + name + ';' +
                              'accession=' + accession + ';')

                ipr_acc = None
                entry = sig['entry']
                if entry is not None:
                    ipr_acc = entry['accession']
                if ipr_acc is not None:
                    attributes = attributes + 'interpro=' + ipr_acc + ';'

                g = gff_template()
                g['seqid'] = rec_key
                g['source'] = lib_name
                g['type'] = 'PfPanther'
                g['start'] = beg
                g['end'] = end
                g['score'] = '.'
                g['strand'] = '+'
                g['phase'] = str(0)
                g['attributes'] = attributes

                gff = gff + gff_text(g)

    return gff


def gff_from_kakapo_ips5_json_file(ss, json_path, gff_path=None):  # noqa

    with open(json_path, 'r') as f:
        json_dict = json.load(f)

    # gff = '##gff-version 3\n'

    gff = (gff_orf_blast_hit(ss, json_dict) +
           gff_orf_bad_blast_hit(ss, json_dict) +
           gff_pfam(ss, json_dict) +
           gff_phobius(ss, json_dict) +
           gff_panther(ss, json_dict))

    if gff_path is not None:
        with open(gff_path, 'w') as f:
            f.write(gff)

    return gff


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
