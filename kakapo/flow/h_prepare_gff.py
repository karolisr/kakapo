# -*- coding: utf-8 -*-

"""Kakapo workflow: Produce GFF3 files."""

import json

from os import remove as osremove
from os.path import exists as ope
from os.path import join as opj

from kakapo.gff3 import gff_from_kakapo_ips5_json_file
from kakapo.helpers import combine_text_files


def gff_from_json(ss, assemblies, dir_prj_ips, dir_prj_transcripts_combined,
                  prj_name, linfo=print):
    if len(assemblies) > 0:
        linfo('Producing GFF3 files [' + ss + ']')

    all_fas_paths = []
    all_gff_paths = []

    combined_fas_path = opj(dir_prj_transcripts_combined, prj_name + '__' + ss + '.fasta')
    combined_gff_path = opj(dir_prj_transcripts_combined, prj_name + '__' + ss + '.gff')

    for a in assemblies:

        if 'transcripts_nt_fasta_file__' + ss not in a:
            continue

        assmbl_name = a['name']
        transcripts_nt_path = a['transcripts_nt_fasta_file__' + ss]

        kakapo_json_path = opj(dir_prj_ips, assmbl_name + '_ann_kakapo__' + ss + '.json')
        ips_json_path = opj(dir_prj_ips, assmbl_name + '_ann_ips__' + ss + '.json')
        json_path = opj(dir_prj_ips, assmbl_name + '_ann__' + ss + '.json')

        gff_path = transcripts_nt_path.replace('.fasta', '.gff')

        ips_json_dict = {}
        kakapo_json_dict = {}

        if ope(ips_json_path):
            with open(ips_json_path, 'r') as f:
                ips_json_dict = json.load(f)

        if ope(kakapo_json_path):
            with open(kakapo_json_path, 'r') as f:
                kakapo_json_dict = json.load(f)

        json_dict = kakapo_json_dict.copy()
        json_dict.update(ips_json_dict)

        with open(json_path, 'w') as f:
            json.dump(json_dict, f, sort_keys=True, indent=4)

        if ope(json_path):
            linfo(assmbl_name)
            gff_from_kakapo_ips5_json_file(ss, json_path, gff_path)

            osremove(json_path)

            all_gff_paths.append(gff_path)
            all_fas_paths.append(transcripts_nt_path)

    combine_text_files(all_fas_paths, combined_fas_path)
    combine_text_files(all_gff_paths, combined_gff_path)
