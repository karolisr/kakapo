"""
EMBL-EBI InterProScan 5.

**Documentation** for Open API Interface at:

https://www.ebi.ac.uk/seqdb/confluence/display/JDSAT/InterProScan+5+Help+and+Documentation
"""

import io
import pickle
import sys
import tarfile
from collections import OrderedDict
from copy import deepcopy
from os import remove
from os.path import exists as ope
from os.path import join as opj
from time import sleep
from xml.etree import ElementTree

from kakapo.tools.config import PICKLE_PROTOCOL
from kakapo.utils.misc import split_seq_defn_for_printing as split_seq_defn
from kakapo.utils.http import get
from kakapo.utils.http import post

IPS_URL = 'https://www.ebi.ac.uk/Tools/services/rest/iprscan5'


def submit_job(email, title, sequence):
    data = {'email': email, 'title': title, 'sequence': sequence}
    url = IPS_URL + '/run'
    response_format = 'text'
    r = post(url, data, response_format)
    if r.status_code != 200:
        job_id = None
    else:
        job_id = r.text
    return job_id


def status(job_id):
    # Possible response values:
    #
    #   RUNNING: the job is currently being processed.
    #   FINISHED: job has finished, and the results can then be retrieved.
    #   ERROR: an error occurred attempting to get the job status.
    #   FAILURE: the job failed.
    #   NOT_FOUND: the job cannot be found.

    url = IPS_URL + '/status/' + job_id
    r = get(url=url, params=None, response_format='text')
    job_status = r.text
    return job_status


def result_types(job_id):
    url = IPS_URL + '/resulttypes/' + job_id
    r = get(url=url, params=None, response_format='xml')
    root = ElementTree.fromstring(r.text)
    types = root.findall('type')
    r = {x.find('identifier').text: {
        'Accept': x.find('mediaType').text} for x in types}
    return r


def _result(job_id, result_type):
    available_result_types = result_types(job_id)

    if result_type in available_result_types:
        response_format = available_result_types[result_type]
    else:
        # raise Exception
        # ToDo: Just skips if, for some reason, the server
        #       did not report an error but also did not
        #       return anything.

        # I noticed that in these cases waiting a little and retrying
        # sometimes works. Should implement wait and retry in this
        # function.
        return None

    url = IPS_URL + '/result/' + job_id + '/' + result_type
    response = get(url=url, params=None, response_format=response_format)
    # Returns unparsed response!
    return response


def result_json(job_id):
    r = _result(job_id, result_type='json')
    # ToDo: Just skips if, for some reason, the server
    #       did not report an error but also did not
    #       return anything.
    if r is None:
        return None
    ##################################################
    return r.json()


def result_sequence(job_id):
    r = _result(job_id, result_type='sequence')
    # ToDo: Just skips if, for some reason, the server
    #       did not report an error but also did not
    #       return anything.
    if r is None:
        return None
    ##################################################
    return r.text


def result_html(job_id, out_dir):
    r = _result(job_id, result_type='htmltarball')
    # ToDo: Just skips if, for some reason, the server
    #       did not report an error but also did not
    #       return anything.
    if r is None:
        return None
    ##################################################

    tar_ref = tarfile.open(fileobj=io.BytesIO(r.content), mode='r:gz')
    tar_ref.extractall(out_dir)
    tar_ref.close()

    return None


def result_gff(job_id, out_file=None):
    r = _result(job_id, result_type='gff')
    # ToDo: Just skips if, for some reason, the server
    #       did not report an error but also did not
    #       return anything.
    if r is None:
        return None
    ##################################################
    gff_text = r.text

    if out_file is not None:
        with open(out_file, 'w') as f:
            f.write(gff_text)

    return gff_text


def job_runner(email, dir_cache, seqs=None, run_id='', parallel_run_count=1,
               max_title_a_len=0, max_run_id_len=0):
    """Run InterProScan 5."""
    max_jobs = int(30 / parallel_run_count)
    delay = 0.333
    cfp = opj(dir_cache, 'ips5_cache_running_' + run_id)

    def load():
        with open(cfp, 'rb') as f:
            c = pickle.load(f)
        return c

    def dump(c):
        with open(cfp, 'wb') as f:
            pickle.dump(c, f, protocol=PICKLE_PROTOCOL)

    seqs_orig = deepcopy(seqs)
    seqs = deepcopy(seqs)

    jobs = None

    if ope(cfp):
        jobs = load()
    else:
        if seqs is not None:
            jobs = {'queue': seqs,
                    'running': OrderedDict(),
                    'finished': OrderedDict(),
                    'error': OrderedDict(),
                    'failure': OrderedDict(),
                    'not_found': OrderedDict(),
                    'other': OrderedDict()}

            dump(jobs)

        else:
            print('No sequences provided.')
            sys.exit(0)

    ##########################################################################

    queue = jobs['queue']
    running = jobs['running']
    finished = jobs['finished']
    # error = jobs['error']
    failure = jobs['failure']
    # not_found = jobs['not_found']
    other = jobs['other']

    retry_list = list()

    queue_size = len(seqs_orig)

    busy = True
    submit_message = True
    while busy:
        dump(jobs)
        sleep(delay)

        if len(queue) > 0:
            if len(running) < max_jobs:
                if submit_message is True:
                    pass
                    # print()
                    # print('Submitting jobs: ' + run_id)
                    # print()
                job_status = 'SUBMITTED '
                title = list(queue.keys()).pop()
                sequence = queue.pop(title)
                job_id = submit_job(email, title, sequence)
                if job_id is None:
                    job_status = 'WILL_RETRY'
                    queue[title] = sequence
                    job_id = ''
                else:
                    running[title] = job_id
                titles_ab = split_seq_defn(title)
                title_a = titles_ab[0]

                msg = (job_status + '  ' +
                       title_a.ljust(max_title_a_len) +
                       run_id.ljust(max_run_id_len) +
                       ' ' * 5 + job_id)

                print(msg)

        if len(running) < max_jobs and len(queue) > 0:
            submit_message = False
            continue

        if len(running) > 0:
            submit_message = True
            # print()
            # print('Looking for finished jobs: ' + run_id)
            # print()
            finished_jobs = False
            sleep(delay + 7)
            job_statuses = {}

            for title in running:
                sleep(delay)
                job_id = running[title]
                job_status = status(job_id)
                job_statuses[title] = {'job_id': job_id,
                                       'job_status': job_status}
                titles_ab = split_seq_defn(title)
                title_a = titles_ab[0]
                # ToDo: Refactor
                if job_status == 'RUNNING':
                    print(' ' * 10 + '- ' + title_a.ljust(max_title_a_len) +
                          run_id.ljust(max_run_id_len) + ' ' * 5 + job_id)
                else:
                    print(' ' * 10 + '+ ' + title_a.ljust(max_title_a_len) +
                          run_id.ljust(max_run_id_len) + ' ' * 5 + job_id)
                    finished_jobs = True

            if finished_jobs is True:
                print()
            else:
                continue

            for title in job_statuses:

                job_id = job_statuses[title]['job_id']
                job_status = job_statuses[title]['job_status']

                if job_status == 'RUNNING':
                    pass

                elif job_status == 'FINISHED':
                    job_id = running.pop(title)
                    if 'error' in result_types(job_id):
                        job_status = 'FAILURE'
                        failure[title] = job_id
                        if retry_list.count(title) < 3:
                            job_status = 'WILL_RETRY'
                            queue[title] = seqs_orig[title]
                            retry_list.append(title)
                    else:
                        finished[title] = job_id

                elif job_status == 'FAILURE':
                    job_id = running.pop(title)
                    failure[title] = job_id

                elif job_status == 'ERROR':
                    continue
                    # job_id = running.pop(title)
                    # error[title] = job_id

                elif job_status == 'NOT_FOUND':
                    job_status = 'WILL_RETRY'
                    job_id = running.pop(title)
                    queue[title] = seqs_orig[title]
                    # not_found[title] = job_id

                else:
                    job_id = running.pop(title)
                    other[title] = job_id

                if job_status != 'RUNNING':
                    progress = round((len(finished) / queue_size) * 100)
                    progress_str = '{:3d}'.format(progress) + '%'
                    titles_ab = split_seq_defn(title)
                    title_a = titles_ab[0]
                    job_status_msg = job_status.ljust(10)
                    msg = (job_status_msg + '  ' +
                           title_a.ljust(max_title_a_len) +
                           run_id.ljust(max_run_id_len) +
                           progress_str.rjust(4) + ' ' + job_id)

                    print(msg)

        if len(running) == 0 and len(queue) == 0:
            busy = False

    if ope(cfp):
        remove(cfp)

    return jobs
