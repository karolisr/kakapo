# -*- coding: utf-8 -*-

"""
EMBL-EBI InterProScan 5

Documentation: https://www.ebi.ac.uk/seqdb/confluence/display/JDSAT/InterProScan+5+Help+and+Documentation#InterProScan5HelpandDocumentation-OpenAPIInterface  # noqa
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

from kakapo.http_k import get
from kakapo.http_k import post

IPS_URL = 'https://www.ebi.ac.uk/Tools/services/rest/iprscan5'


def submit_job(email, title, sequence):  # noqa
    data = {'email': email, 'title': title, 'sequence': sequence}
    url = IPS_URL + '/run'
    response_format = 'text'
    r = post(url, data, response_format)
    job_id = r.text
    return job_id


def status(job_id):  # noqa

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


def result_types(job_id):  # noqa
    url = IPS_URL + '/resulttypes/' + job_id
    r = get(url=url, params=None, response_format='xml')
    root = ElementTree.fromstring(r.text)
    types = root.findall('type')
    r = {x.find('identifier').text: {
         'Accept': x.find('mediaType').text} for x in types}
    return r


def _result(job_id, result_type):

    valid_result_types = result_types(job_id)

    if result_type in valid_result_types:
        response_format = valid_result_types[result_type]
    else:
        raise Exception

    url = IPS_URL + '/result/' + job_id + '/' + result_type
    response = get(url=url, params=None, response_format=response_format)
    # Returns unparsed response!
    return response


def result_json(job_id):  # noqa
    r = _result(job_id, result_type='json')
    return r.json()


def result_sequence(job_id):  # noqa
    r = _result(job_id, result_type='sequence')
    return r.text


def result_html(job_id, out_dir):  # noqa
    r = _result(job_id, result_type='htmltarball')

    tar_ref = tarfile.open(fileobj=io.BytesIO(r.content), mode='r:gz')
    tar_ref.extractall(out_dir)
    tar_ref.close()

    return None


def result_gff(job_id, out_file=None):  # noqa
    r = _result(job_id, result_type='gff')
    gff_text = r.text

    if out_file is not None:
        with open(out_file, 'w') as f:
            f.write(gff_text)

    return gff_text


def job_runner(email, dir_cache, seqs=None, logger=print):
    """Run InterProScan 5"""
    ##########################################################################
    MAX_JOBS = 25
    DELAY = 5
    CFP = opj(dir_cache, 'ips5_cache')
    ##########################################################################

    ##########################################################################
    def load():
        with open(CFP, 'rb') as f:
            c = pickle.load(f)
        return c

    def dump(c):
        with open(CFP, 'wb') as f:
            pickle.dump(c, f, protocol=2)
    ##########################################################################

    seqs = deepcopy(seqs)

    jobs = None

    if ope(CFP):
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
            # print('No sequences provided.')
            sys.exit(0)

    ##########################################################################

    queue = jobs['queue']
    running = jobs['running']
    finished = jobs['finished']
    error = jobs['error']
    failure = jobs['failure']
    not_found = jobs['not_found']
    other = jobs['other']

    ###
    # Nicer printing
    max_title_len = max([len(x) for x in list(queue.keys())])
    ###

    busy = True
    while busy:

        dump(jobs)

        sleep(DELAY)

        if len(queue) > 0:
            if len(running) < MAX_JOBS:
                title = list(queue.keys()).pop()
                sequence = queue.pop(title)
                job_id = submit_job(email, title, sequence)
                running[title] = job_id
                title_len = len(title)
                title_len_diff = max_title_len - title_len
                msg = '' + title + title_len_diff * ' ' + \
                      '\t' + job_id + '\t' + 'SUBMITTED'
                logger(msg)

        if len(running) > 0:

            job_statuses = {}

            for title in running:
                job_id = running[title]
                job_status = status(job_id)

                job_statuses[title] = {'job_id': job_id,
                                       'job_status': job_status}

            for title in job_statuses:

                job_id = job_statuses[title]['job_id']
                job_status = job_statuses[title]['job_status']

                if job_status == 'RUNNING':
                    pass

                elif job_status == 'FINISHED':
                    job_id = running.pop(title)
                    finished[title] = job_id

                elif job_status == 'ERROR':
                    job_id = running.pop(title)
                    error[title] = job_id

                elif job_status == 'FAILURE':
                    job_id = running.pop(title)
                    failure[title] = job_id

                elif job_status == 'NOT_FOUND':
                    job_id = running.pop(title)
                    not_found[title] = job_id

                else:
                    job_id = running.pop(title)
                    other[title] = job_id

                if job_status != 'RUNNING':
                    title_len = len(title)
                    title_len_diff = max_title_len - title_len
                    msg = '' + title + title_len_diff * ' ' + \
                          '\t' + job_id + '\t' + job_status
                    logger(msg)

        if len(running) == 0 and len(queue) == 0:
            busy = False

    if ope(CFP):
        remove(CFP)

    return jobs
