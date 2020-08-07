"""Kakapo workflow: Run InterProScan."""

import json
import pickle

from collections import OrderedDict
from os import remove as osremove
from os.path import exists as ope
from os.path import join as opj
from time import sleep

from kakapo.tools.bioio import read_fasta
from kakapo.tools.bioio import seq_records_to_dict
from kakapo.tools.config import PICKLE_PROTOCOL
from kakapo.tools.ebi_iprscan5 import job_runner
from kakapo.tools.ebi_iprscan5 import result_json
from kakapo.tools.seq import SEQ_TYPE_AA

from kakapo.utils.misc import split_seq_defn_for_printing as split_seq_defn


def run_inter_pro_scan(ss, assemblies, email, dir_prj_ips, dir_cache_prj,
                       parallel_run_count, max_title_a_len, max_run_id_len,
                       linfo=print):

    delay = 0.25

    for a in assemblies:

        if 'transcripts_aa_orf_fasta_file__' + ss not in a:
            continue

        aa_file = a['transcripts_aa_orf_fasta_file__' + ss]

        if aa_file is None:
            continue

        assmbl_name = a['name']

        json_dump_file_path = opj(dir_prj_ips, assmbl_name + '_ann_ips__' +
                                  ss + '.json')

        if ope(json_dump_file_path):
            linfo('InterProScan results for assembly ' + assmbl_name +
                  ', search strategy ' + ss + ' have already been downloaded.')
            continue
        else:
            print()
            linfo('Running InterProScan on translated ' + ss +
                  ' from ' + assmbl_name + '.')
            print()

        seqs = seq_records_to_dict(read_fasta(aa_file, SEQ_TYPE_AA))

        # Filter all ORFs except the first one.
        for seq_def in tuple(seqs.keys()):
            seq_def_prefix = seq_def.split(' ')[0]
            if not seq_def_prefix.endswith('ORF001'):
                del seqs[seq_def]

        seqs = OrderedDict(sorted(seqs.items(),
                                  key=lambda x: x[0].split(' ')[1],
                                  reverse=True))

        run_id = ss + '_' + assmbl_name

        _ = opj(dir_cache_prj, 'ips5_cache_done_' + run_id)

        if ope(_):
            with open(_, 'rb') as f:
                jobs = pickle.load(f)

        else:
            jobs = job_runner(email=email, dir_cache=dir_cache_prj,
                              seqs=seqs, run_id=run_id,
                              parallel_run_count=parallel_run_count,
                              max_title_a_len=max_title_a_len,
                              max_run_id_len=max_run_id_len)

            with open(_, 'wb') as f:
                pickle.dump(jobs, f, protocol=PICKLE_PROTOCOL)

        print()
        linfo('Downloading InterProScan results for ' + ss +
              ' in ' + assmbl_name + '.')
        print()

        all_ips_results = {}

        # Nicer printing
        for i, job in enumerate(jobs['finished']):

            job_id = jobs['finished'][job]

            titles_ab = split_seq_defn(job)
            title_a = titles_ab[0]

            progress = round(((i + 1) / len(jobs['finished'])) * 100)
            progress_str = '{:3d}'.format(progress) + '%'

            msg = (' ' * 12 +
                   title_a.ljust(max_title_a_len) +
                   run_id.ljust(max_run_id_len) +
                   progress_str.rjust(4) + ' ' + job_id)

            print(msg)

            sleep(delay)

            ips_json = result_json(job_id)
            if ips_json is None:
                continue
            # ips_version = ips_json['interproscan-version']
            ips_json = ips_json['results']

            # These fields are set to 'EMBOSS_001' by default
            # Delete them
            del ips_json[0]['xref']

            job_no_def = job.split(' ')[0]

            all_ips_results[job_no_def] = ips_json

        with open(json_dump_file_path, 'w') as f:
            json.dump(all_ips_results, f, sort_keys=True, indent=4)

        # Removes cached jobs file.
        osremove(_)
