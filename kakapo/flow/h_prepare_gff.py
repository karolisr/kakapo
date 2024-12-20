"""Kakapo workflow: Produce GFF3 files."""

from os.path import exists as ope
from os.path import join as opj

from kakapo.tools.gff3 import gff_from_json_dict
from kakapo.utils.misc import combine_text_files
from kakapo.tools.seq_annotations import parse_kakapo_json_file
from kakapo.tools.seq_annotations import merge_kakapo_and_ips_annotations


def gff_from_json(ss, assemblies, dir_prj_ips, dir_prj_transcripts_combined,
                  prj_name):

    all_fas_paths = []
    all_gff_paths = []

    combined_fas_path = opj(dir_prj_transcripts_combined, prj_name
                            + '__' + ss + '.fasta')
    combined_gff_path = opj(dir_prj_transcripts_combined, prj_name
                            + '__' + ss + '.gff')

    for a in assemblies:

        if 'transcripts_nt_fasta_file__' + ss not in a:
            continue

        assmbl_name = a['name']
        transcripts_nt_path = a['transcripts_nt_fasta_file__' + ss]

        kakapo_json_path = opj(dir_prj_ips, assmbl_name + '_ann_kakapo__'
                               + ss + '.json')
        ips_json_path = opj(dir_prj_ips, assmbl_name + '_ann_ips__'
                            + ss + '.json')

        gff_path = transcripts_nt_path.replace('.fasta', '.gff')

        json_dict = {}

        if ope(ips_json_path) and ope(kakapo_json_path):
            json_dict = merge_kakapo_and_ips_annotations(kakapo_json_path,
                                                         ips_json_path)

        elif ope(kakapo_json_path):
            json_dict = parse_kakapo_json_file(kakapo_json_path)

        # Log.inf('Producing GFF3 files for ' + ss + ' in ' + assmbl_name + '.')
        gff_from_json_dict(json_dict, gff_path)

        all_gff_paths.append(gff_path)
        all_fas_paths.append(transcripts_nt_path)

    combine_text_files(all_fas_paths, combined_fas_path)
    combine_text_files(all_gff_paths, combined_gff_path)
