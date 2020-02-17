# -*- coding: utf-8 -*-

"""Kakapo workflow: Run InterProScan."""

import json
import pickle

from collections import OrderedDict
from os import remove as osremove
from os.path import exists as ope
from os.path import join as opj
from time import sleep

from kakapo.bioio import read_fasta
from kakapo.config import PICKLE_PROTOCOL
from kakapo.ebi_iprscan5 import job_runner
from kakapo.ebi_iprscan5 import result_json
from kakapo.helpers import split_seq_defn_for_printing as split_seq_defn


def run_inter_pro_scan(ss, assemblies, email, dir_prj_ips, dir_cache_prj,
                       linfo=print):
    delay = 0.25

    if len(assemblies) > 0:
        linfo('Running InterProScan on translated transcripts [' + ss + ']')

    for a in assemblies:

        if 'transcripts_aa_orf_fasta_file__' + ss not in a:
            continue

        aa_file = a['transcripts_aa_orf_fasta_file__' + ss]

        if aa_file is None:
            continue

        assmbl_name = a['name']

        json_dump_file_path = opj(dir_prj_ips, assmbl_name + '_ann_ips__' + ss + '.json')

        if ope(json_dump_file_path):
            linfo('InterProScan results for assembly ' + assmbl_name + ', ' +
                  ss + ' have already been downloaded [' + ss + ']')
            continue

        seqs = read_fasta(aa_file)
        seqs = OrderedDict(sorted(seqs.items(),
                                  key=lambda x: x[0].split(' ')[1],
                                  reverse=True))

        _ = opj(dir_cache_prj, assmbl_name + '_ips_jobs__' + ss)

        if ope(_):
            with open(_, 'rb') as f:
                jobs = pickle.load(f)

        else:
            jobs = job_runner(email=email, dir_cache=dir_cache_prj,
                              seqs=seqs, logger=linfo)

            with open(_, 'wb') as f:
                pickle.dump(jobs, f, protocol=PICKLE_PROTOCOL)

        print()
        linfo('Downloading InterProScan results for transcripts in ' +
              assmbl_name + ' [' + ss + ']')
        print()

        all_ips_results = {}

        # Nicer printing
        max_title_a_len = 2 + max([len(split_seq_defn(x)[0]) for x in list(jobs['finished'].keys())])
        max_title_b_len = 2 + max([len(split_seq_defn(x)[1]) for x in list(jobs['finished'].keys())])

        for i, job in enumerate(jobs['finished']):

            job_id = jobs['finished'][job]

            titles_ab = split_seq_defn(job)
            title_a = titles_ab[0]
            title_b = titles_ab[1]

            progress = round(((i + 1) / len(jobs['finished'])) * 100)
            progress_str = '{:3d}'.format(progress) + '%'

            msg = (' ' * 12 +
                   title_a.ljust(max_title_a_len) +
                   title_b.ljust(max_title_b_len) +
                   progress_str + ' ' + job_id)

            linfo(msg)

            sleep(delay)

            ips_json = result_json(job_id)
            # ips_version = ips_json['interproscan-version']
            ips_json = ips_json['results']

            # These fields are set to 'EMBOSS_001' by default
            # Delete them
            del ips_json[0]['xref']

            job_no_def = job.split(' ')[0]

            # kakapo annotations
            ips_json[0]['kakapo_annotations__' + ss] = a['annotations__' + ss][job_no_def]

            all_ips_results[job_no_def] = ips_json

        print()

        with open(json_dump_file_path, 'w') as f:
            json.dump(all_ips_results, f, sort_keys=True, indent=4)

        # Removes cached jobs file.
        osremove(_)
