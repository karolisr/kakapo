# -*- coding: utf-8 -*-
"""
Wraps with NCBI's Entrez Programming Utilities (E-utilities).

More information on E-utilities at:
    http://www.ncbi.nlm.nih.gov/books/NBK25497

Database names and unique identifiers returned can be found here:
    http://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly

"""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import with_statement

import re
import os

from time import sleep
from xml.etree import ElementTree
from math import ceil

from kakapo.http_k import get
from kakapo.http_k import post

from kakapo.parsers import parse_efetch_sra_csv_text
from kakapo.parsers import parse_gbseq_xml_text
from kakapo.parsers import parse_esummary_xml_text

ENTREZ_BASE_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
DELAY = 1


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


def esearch(term, db, api_key=None):  # noqa

    if api_key is None:
        api_key = _check_for_api_key()

    eutil = 'esearch.fcgi'
    url = ENTREZ_BASE_URL + eutil

    params = {'db': db, 'term': term, 'idtype': 'acc', 'rettype': 'count',
              'retmode': 'json', 'usehistory': 'y', 'api_key': api_key}

    response = get(url, params, 'json')
    parsed = response.json()
    esearch_result = parsed['esearchresult']
    if 'ERROR' in esearch_result:
        raise Exception(esearch_result['ERROR'])

    record_count = int(esearch_result['count'])

    ret_max = 5000
    query_key = None
    web_env = None
    id_set = set()

    for ret_start in range(0, record_count, ret_max):

        if ret_start > 0:
            sleep(DELAY)

        params = {'db': db, 'term': term, 'idtype': 'acc',
                  'retstart': ret_start, 'retmax': ret_max,
                  'rettype': 'uilist', 'retmode': 'json', 'usehistory': 'y',
                  'api_key': api_key}

        if query_key is not None:
            params['query_key'] = query_key

        if web_env is not None:
            params['WebEnv'] = web_env

        response = get(url, params, 'json')
        parsed = response.json()
        data = parsed['esearchresult']

        id_set = id_set | set(data['idlist'])
        query_key = data['querykey']
        web_env = data['webenv']

    id_tuple = tuple(sorted(id_set))

    return_dict = {
        'db': db,
        'record_count': record_count,
        'ids': id_tuple,
        'query_key': query_key,
        'web_env': web_env}

    return return_dict


def epost(ids, db, api_key=None):  # noqa

    if api_key is None:
        api_key = _check_for_api_key()

    eutil = 'epost.fcgi'
    url = ENTREZ_BASE_URL + eutil
    data = {'db': db, 'id': ','.join(ids), 'api_key': api_key}

    response = post(url, data, 'xml')
    root = ElementTree.fromstring(response.text)

    query_key = None
    web_env = None

    for child in root:
        if child.tag == 'QueryKey':
            query_key = child.text
        if child.tag == 'WebEnv':
            web_env = child.text

    record_count = len(ids)

    return_dict = {
        'db': db,
        'record_count': record_count,
        'query_key': query_key,
        'web_env': web_env}

    return return_dict


def efetch_data(data, parser, ret_type=None, ret_mode='xml', api_key=None):  # noqa

    if api_key is None:
        api_key = _check_for_api_key()

    eutil = 'efetch.fcgi'
    url = ENTREZ_BASE_URL + eutil

    db = data['db']
    record_count = data['record_count']
    query_key = data['query_key']
    web_env = data['web_env']

    ret_max = 500
    ret_start = 0

    return_list = []

    for ret_start in range(0, record_count, ret_max):

        if ret_start > 0:
            sleep(DELAY)

        params = {'db': db, 'query_key': query_key, 'WebEnv': web_env,
                  'retstart': ret_start, 'retmax': ret_max,
                  'rettype': ret_type, 'retmode': ret_mode, 'usehistory': 'y',
                  'api_key': api_key}

        response = get(url, params, ret_mode)
        parsed = parser(response.text)

        for item in parsed:
            return_list.append(item)

    return return_list


def efetch(db, params, ret_type=None, ret_mode='xml', api_key=None):  # noqa

    if api_key is None:
        api_key = _check_for_api_key()

    eutil = 'efetch.fcgi'
    url = ENTREZ_BASE_URL + eutil

    params_all = {'db': db, 'rettype': ret_type, 'retmode': ret_mode,
                  'api_key': api_key}

    params_all.update(params)

    response = get(url, params_all, ret_mode)

    return response.text


def esummary(data, api_key=None):  # noqa

    if api_key is None:
        api_key = _check_for_api_key()

    eutil = 'esummary.fcgi'

    db = data['db']
    record_count = data['record_count']
    query_key = data['query_key']
    web_env = data['web_env']

    ret_max = 500

    return_list = []

    for ret_start in range(0, record_count, ret_max):

        if ret_start > 0:
            sleep(DELAY)

        url = ENTREZ_BASE_URL + eutil

        params = {'db': db, 'query_key': query_key, 'WebEnv': web_env,
                  'retstart': ret_start, 'retmax': ret_max, 'retmode': 'json'}

        response = get(url, params, 'json')
        parsed = response.json()['result']
        keys = parsed['uids']

        for k in keys:
            return_list.append(parsed[k])

    return return_list


def esearch_epost(term, db):  # noqa
    if type(term) in (list, tuple):
        term = ' OR '.join(term)
    esearch_results = esearch(term, db)
    id_list = esearch_results['ids']
    if len(id_list) > 0:
        epost_results = epost(id_list, db)
    else:
        epost_results = None
    return epost_results


def summary(term, db):  # noqa
    epost_results = esearch_epost(term, db)
    esummary_results = esummary(epost_results, parse_esummary_xml_text)
    return esummary_results


def taxids_for_acc(accessions, db):  # noqa
    summ = summary(accessions, db)
    ret_dict = dict()
    for x in summ:
        acc = x['accessionversion']
        taxid = x['taxid']
        ret_dict[acc] = taxid

    return ret_dict

def dnld_seqs(term, db):  # noqa
    epost_results = esearch_epost(term, db)
    efetch_results = efetch_data(epost_results, parse_gbseq_xml_text, 'gb')
    return efetch_results


def dnld_seqs_gb_format(term, db):  # noqa
    epost_results = esearch_epost(term, db)
    efetch_results = efetch_data(epost_results, lambda x: x.split('\n\n')[:-1],
                                 'gb', 'text')
    return efetch_results


def dnld_seqs_fasta_format(term, db):  # noqa
    epost_results = esearch_epost(term, db)
    efetch_results = efetch_data(epost_results, lambda x: x.split('\n>'),
                                 'fasta', 'text')

    ret_list = []
    for x in efetch_results:
        x = x.strip('>\n')
        x = '>' + x
        ret_list.append(x)

    ret_list = sorted(list(set(ret_list)))

    return ret_list


def dnld_cds_nt_fasta(term):  # noqa
    epost_results = esearch_epost(term, 'nuccore')
    efetch_results = efetch_data(epost_results, lambda x: x.split('\n>'),
                                 'fasta_cds_na', 'text')

    ret_list = []
    for x in efetch_results:
        x = x.strip('>\n')
        x = '>' + x
        ret_list.append(x)

    ret_list = sorted(list(set(ret_list)))

    return ret_list


def cds_acc_for_prot_acc(prot_accessions):  # noqa
    ret_dict = dict()
    prot_dict_list = dnld_seqs(prot_accessions, 'protein')

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

        epost_results = esearch_epost(term, 'sra')
        if epost_results is None:
            continue

        efetch_results = efetch_data(epost_results, parse_efetch_sra_csv_text,
                                     'runinfo', 'csv')

        ret_list = ret_list + efetch_results

    return ret_list
