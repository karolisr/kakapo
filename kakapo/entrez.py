# -*- coding: utf-8 -*-
"""
Wraps with NCBI's Entrez Programming Utilities (E-utilities).

More information on E-utilities at:
    http://www.ncbi.nlm.nih.gov/books/NBK25497

Database names and unique identifiers returned can be found here:
    http://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly

"""

import re
import os

from io import StringIO
from math import ceil
from time import sleep
from xml.etree import ElementTree

from kakapo.http_k import get
from kakapo.http_k import post

from kakapo.parsers import parse_efetch_sra_csv_text
from kakapo.parsers import parse_gbseq_xml_text
from kakapo.parsers import parse_esummary_xml_text
from kakapo.bioio import read_fasta

ENTREZ_BASE_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
DELAY = 0.5


def _check_for_api_key():

    global_variables = globals()

    if 'ENTREZ_KEY' in global_variables:
        ncbi_api_key = global_variables['ENTREZ_KEY']
    elif 'ENTREZ_KEY' in os.environ:
        ncbi_api_key = os.environ['ENTREZ_KEY']
    else:
        print('Warning: ENTREZ_KEY is not defined.')
        ncbi_api_key = None

    return ncbi_api_key


term = '"RefSeq"[Keyword] AND "arabidopsis"[Primary Organism] AND "chloroplast"[filter]'
db = 'nuccore'

def esearch(term, db, api_key=None, ret_type='uilist'):  # noqa

    # rettype='uilist'
    # rettype='acc'

    if api_key is None:
        api_key = _check_for_api_key()

    eutil = 'esearch.fcgi'
    url = ENTREZ_BASE_URL + eutil

    params = {'db': db, 'term': term, 'idtype': 'acc', 'rettype': ret_type,
              'retmode': 'json', 'usehistory': 'y', 'api_key': api_key,
              'retmax': 0}

    response = get(url, params, 'json')
    parsed = response.json()
    data = parsed['esearchresult']

    counts = [int(data['count'])]
    query_keys = [int(data['querykey'])]
    web_env = data['webenv']

    return_dict = {
        'db': db,
        'counts': counts,
        'query_keys': query_keys,
        'web_env': web_env}

    return return_dict


def epost(ids, db, api_key=None):  # noqa

    if api_key is None:
        api_key = _check_for_api_key()

    eutil = 'epost.fcgi'
    url = ENTREZ_BASE_URL + eutil

    db = db
    ids_count = len(ids)
    counts = list()
    query_keys = list()
    web_env = None
    max_count = 100

    for strt in range(0, ids_count, max_count):

        # print('epost ->', strt + 1, ids_count)

        end = strt + min(ids_count - strt, max_count)
        counts.append(end - strt)
        ids_temp = ids[strt:end]

        data = {'db': db, 'id': ','.join(ids_temp), 'api_key': api_key,
                'WebEnv': web_env, 'usehistory': 'y'}

        response = post(url, data, 'xml')

        root = ElementTree.fromstring(response.text)

        for child in root:
            if child.tag == 'QueryKey':
                query_keys.append(int(child.text))
            if child.tag == 'WebEnv':
                web_env = child.text

        sleep(DELAY)

    return_dict = {
        'db': db,
        'counts': counts,
        'query_keys': query_keys,
        'web_env': web_env}

    return return_dict


def efetch(data, parser, ret_type=None, ret_mode='xml', api_key=None):  # noqa

    # TODO: Print progress.

    if api_key is None:
        api_key = _check_for_api_key()

    eutil = 'efetch.fcgi'
    url = ENTREZ_BASE_URL + eutil

    db = data['db']
    counts = data['counts']
    query_keys = data['query_keys']
    web_env = data['web_env']

    ret_max = 10
    ret_start = 0

    return_list = []

    for i in range(0, len(query_keys)):
        query_key = query_keys[i]
        count = counts[i]

        for ret_start in range(0, count, ret_max):

            # print('efetch ->', i + 1, len(query_keys), ret_start, count)

            if ret_start > 0:
                sleep(DELAY)

            params = {'db': db, 'query_key': query_key, 'WebEnv': web_env,
                      'retstart': ret_start, 'retmax': ret_max,
                      'rettype': ret_type, 'retmode': ret_mode,
                      'usehistory': 'y', 'api_key': api_key,
                      'idtype': 'acc'}

            response = get(url, params, ret_mode)
            parsed = parser(response.text)

            for item in parsed:
                # print(item)
                return_list.append(item)

    return return_list


def esummary(accs, db, api_key=None):  # noqa

    if api_key is None:
        api_key = _check_for_api_key()

    eutil = 'esummary.fcgi'
    url = ENTREZ_BASE_URL + eutil

    data = epost(accs, db, api_key=None)

    db = data['db']
    counts = data['counts']
    query_keys = data['query_keys']
    web_env = data['web_env']

    ret_max = 500
    ret_start = 0

    return_list = []

    for i in range(0, len(query_keys)):
        query_key = query_keys[i]
        count = counts[i]

        for ret_start in range(0, count, ret_max):

            if ret_start > 0:
                sleep(DELAY)

            params = {'db': db, 'query_key': query_key, 'WebEnv': web_env,
                      'retstart': ret_start, 'retmax': ret_max,
                      'retmode': 'json',
                      'usehistory': 'y', 'api_key': api_key,
                      'idtype': 'acc'}

            response = get(url, params, 'json')
            parsed = response.json()['result']
            keys = parsed['uids']

            for k in keys:
                return_list.append(parsed[k])

    return return_list


def taxids_for_accs(accs, db):  # noqa
    summ = esummary(accs, db)
    ret_dict = dict()
    for x in summ:
        acc = x['accessionversion']
        taxid = x['taxid']
        ret_dict[acc] = taxid

    return ret_dict


def accessions(data):  # noqa
    ids = efetch(data=data, ret_type='uilist', ret_mode='text',
                 parser=lambda x: x.strip().split('\n'))
    return ids


def dnld_seqs(accs, db):  # noqa
    data = epost(accs, db)
    efetch_results = efetch(data, parse_gbseq_xml_text, 'gb')
    return efetch_results


def dnld_seqs_gb_format(accs, db):  # noqa
    data = epost(accs, db)
    efetch_results = efetch(data, lambda x: x.split('\n\n')[:-1], 'gb', 'text')
    return efetch_results


def _process_fasta_efetch_results(efetch_results):  # noqa
    ret_list = []
    for x in efetch_results:
        x = x.strip('>\n')
        x = '>' + x
        ret_list.append(x)
    ret_list = sorted(list(set(ret_list)))
    parsed_fasta = read_fasta(StringIO('\n'.join(ret_list)))
    return parsed_fasta


def dnld_seqs_fasta_format(accs, db):  # noqa
    data = epost(accs, db)
    efetch_results = efetch(data, lambda x: x.split('\n>'), 'fasta', 'text')
    return _process_fasta_efetch_results(efetch_results)


def dnld_cds_nt_fasta(accs):  # noqa
    data = epost(accs, 'nuccore')
    efetch_results = efetch(data, lambda x: x.split('\n>'),
                            'fasta_cds_na', 'text')
    return _process_fasta_efetch_results(efetch_results)


def cds_acc_for_prot_acc(prot_accessions):  # noqa
    prot_dict_list = dnld_seqs(prot_accessions, 'protein')
    ret_dict = dict()
    for rec in prot_dict_list:
        acc = rec['accession']
        ver = rec['version']
        accv = acc + '.' + ver
        dbsource = rec['db_source']
        dbsource = re.findall('.*accession\s+(.*)', dbsource)
        if len(dbsource) != 0:
            dbsource = dbsource[0]
        else:
            dbsource = None
        ret_dict[accv] = dbsource
    return ret_dict


def sra_run_info(acc_list):  # noqa
    assert type(acc_list) in (tuple, list, set)
    ret_list = []
    page_size = 75
    tot_acc = len(acc_list)
    temp_acc_list = acc_list
    pages = 1
    if tot_acc > page_size:
        pages = ceil(tot_acc / page_size)
    for p in range(0, pages):
        beg = page_size * p
        end = beg + page_size
        temp_acc_list = acc_list[beg:end]
        term = ' OR '.join(temp_acc_list)
        esearch_results = esearch(term, 'sra')
        if esearch_results is None:
            continue
        efetch_results = efetch(esearch_results, parse_efetch_sra_csv_text,
                                'runinfo', 'csv')
        ret_list = ret_list + efetch_results
    return ret_list
