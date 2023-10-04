"""Sequence annotation parsers."""

import json
from collections import OrderedDict
from copy import deepcopy


def parse_ips5_json_file(json_path):
    with open(json_path, 'r') as f:
        json_dict = json.load(f)

    ann_types_supported = ('PFAM', 'PANTHER', 'PHOBIUS')

    phobius_translations = {'SIGNAL_PEPTIDE': 'SigPept',
                            'SIGNAL_PEPTIDE_N_REGION': 'SigPept1_N',
                            'SIGNAL_PEPTIDE_H_REGION': 'SigPept2_H',
                            'SIGNAL_PEPTIDE_C_REGION': 'SigPept3_C',
                            'NON_CYTOPLASMIC_DOMAIN': 'PhobNonCytoPlasm',
                            'TRANSMEMBRANE': 'PhobTransMemb',
                            'CYTOPLASMIC_DOMAIN': 'PhobCytoPlasm'}

    all_annotations = OrderedDict()

    for orf in json_dict:
        orf_annotations = list()
        # Assumption json_dict[orf] is always a list of length 1
        ips_data = json_dict[orf][0]
        ips_annotations = ips_data['matches']
        for ann in ips_annotations:
            # Assumption ann[locations] is always a list of length 1
            locations = ann['locations'][0]
            signature = ann['signature']
            name = signature['name']
            accession = signature['accession']
            description = signature['description']
            entry = signature['entry']
            ann_type = signature['signatureLibraryRelease']['library']
            # ann_type_version = signature['signatureLibraryRelease']['version']

            if ann_type in ann_types_supported:

                temp = orf.split('__ORF')
                orf_id = 'ORF' + temp[1]

                attributes = 'name=' + orf_id + ', '
                ann_type_new = ann_type

                if ann_type == 'PANTHER':
                    if ':' in accession:
                        continue

                    if name.lower() == 'family not named':
                        name = ''
                    else:
                        name = ': ' + name

                    attributes += (accession + name + ';'
                                   + 'accession=' + accession + ';')

                if ann_type == 'PFAM':
                    attributes += (accession + ': '
                                   + description + ';'
                                   + 'accession=' + accession + ';')

                if ann_type == 'PHOBIUS':
                    attributes += (name + ';'
                                   'description=' + description + ';')

                    ann_type_new = phobius_translations[accession]

                ipr_acc = None
                if entry is not None:
                    ipr_acc = entry['accession']
                if ipr_acc is not None:
                    attributes += 'interpro=' + ipr_acc + ';'

                beg = int(locations['start'])
                end = int(locations['end'])

                annotation = {'type': ann_type_new,
                              'beg': beg,
                              'end': end,
                              'attributes': attributes}

                orf_annotations.append(annotation)

        all_annotations[orf] = orf_annotations

    return all_annotations


def parse_kakapo_json_file(json_path):
    with open(json_path, 'r') as f:
        json_dict = json.load(f)

    all_annotations = OrderedDict()
    for transcript in json_dict:
        rec = json_dict[transcript][0]
        ann_key = tuple(rec.keys())[0]
        ann = rec[ann_key]
        all_annotations[transcript] = ann

    return all_annotations


def merge_kakapo_and_ips_annotations(json_path_kakapo, json_path_ips):
    ann_kakapo = parse_kakapo_json_file(json_path_kakapo)
    ann_ips = parse_ips5_json_file(json_path_ips)

    ann_all = deepcopy(ann_kakapo)

    for ips_key in ann_ips:
        temp = ips_key.split('__ORF')
        transcript = temp[0]
        orf = 'ORF' + temp[1]

        kakapo_good_orfs = ann_all[transcript]['orfs_good']

        orf_ips_ann = ann_ips[ips_key]

        for ips_ann in orf_ips_ann:

            beg = ips_ann['beg']
            end = ips_ann['end']

            ips_ann['beg'] = beg * 3 - 2 + kakapo_good_orfs[orf]['orf_begin']
            ips_ann['end'] = end * 3 + kakapo_good_orfs[orf]['orf_begin']

        kakapo_good_orfs[orf]['orf_ips_ann'] = orf_ips_ann

    return ann_all
