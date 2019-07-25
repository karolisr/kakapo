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

from time import sleep
from xmltodict import parse as parse_xml

from kakapo.http_k import get
from kakapo.http_k import post

from kakapo.parsers import parse_efetch_sra_xml_text
from kakapo.parsers import parse_gbseq_xml_text
from kakapo.parsers import parse_esummary_xml_text

ENTREZ_BASE_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
DELAY = 0


def esearch(db, term):
    """

    Wrap ESearch E-utility.

    :param db: Name of the Entrez database to search.
    :type db: str

    :param term: Search terms.
    :type term: str

    :returns: A dictionary with these keys: Database, Count, IdList, QueryKey,
        WebEnv
    :rtype: dict
    """
    eutil = 'esearch.fcgi'

    url = ENTREZ_BASE_URL + eutil

    params = {'db': db, 'term': term, 'rettype': 'count'}
    response = get(url, params, 'xml')
    total_count = int(parse_xml(response.text)['eSearchResult']['Count'])

    # Now download the uids
    retmax = 5000
    retstart = 0
    query_key = ''
    web_env = ''
    id_set = set()
    return_value = list()

    for retstart in range(0, total_count, retmax):

        params = {'db': db, 'query_key': query_key, 'WebEnv': web_env,
                  'retstart': str(retstart), 'retmax': str(retmax),
                  'usehistory': 'y', 'term': term}

        response = get(url, params, 'xml')

        data = parse_xml(response.text)['eSearchResult']

        if total_count == 1:
            id_set.add(data['IdList']['Id'])
        else:
            for uid in data['IdList']['Id']:
                id_set.add(uid)

        query_key = data['QueryKey']
        web_env = data['WebEnv']

    id_list = list(id_set)
    id_list.sort()

    # if count >= 99999:
    #     message = (
    #         'There are more than 99,999 unique identifiers: {c}.')
    #     message = message.format(c=count)
    #     raise Error(message)

    return_value = {
        'Database': db,
        'Count': total_count,
        'IdList': id_list,
        'QueryKey': query_key,
        'WebEnv': web_env}

    sleep(DELAY)

    return return_value


def epost(db, id_list):
    """

    Wrap EPost E-utility.

    :param db: Name of the Entrez database to search.
    :type db: str

    :param id_list: List of unique identifiers.
    :type id_list: list

    :returns: A dictionary with these keys: Database, Count, QueryKey, WebEnv
    :rtype: dict

    """
    eutil = 'epost.fcgi'

    id_list_string = ','.join(id_list)

    url = ENTREZ_BASE_URL + eutil

    data = {'db': db, 'id': id_list_string}

    response = post(url, data, 'xml')

    results = parse_xml(response.text)['ePostResult']

    query_key = results['QueryKey']
    web_env = results['WebEnv']

    count = len(id_list)

    return_value = {
        'Database': db,
        'Count': count,
        'QueryKey': query_key,
        'WebEnv': web_env}

    sleep(DELAY)

    return return_value


def efetch(data, parser, ret_type, retmode='xml'):
    """

    Wrap EFetch E-utility.

    :param data: A dictionary returned by :func:`esearch` or :func:`epost` with
        these keys: Database, Count, QueryKey, WebEnv
    :type data: dict

    :param parser: A function that will be called to interpret downloaded data.
        This function may be called several times as :func:`efetch` downloads
        downloads data in batches.
    :type parser: function

    :param ret_type: Retrieval type. This parameter specifies the record view
        returned, such as Abstract or MEDLINE from PubMed, or GenPept or FASTA
        from protein.
    :type ret_type: str

    :returns: A list of one or more items which will be of the type produced by
        the parser.
    :rtype: list
    """
    eutil = 'efetch.fcgi'

    db = data['Database']
    query_key = data['QueryKey']
    web_env = data['WebEnv']
    count = data['Count']

    retmax = 500
    retstart = 0

    return_value = []

    for retstart in range(0, count, retmax):

        url = ENTREZ_BASE_URL + eutil

        params = {'db': db, 'query_key': query_key, 'WebEnv': web_env,
                  'retstart': str(retstart), 'retmax': str(retmax),
                  'rettype': ret_type, 'retmode': retmode}

        response = get(url, params, retmode)

        parsed = parser(response.text)

        for item in parsed:
            return_value.append(item)

    # ret_mode Retrieval mode. This parameter specifies the data format
    #     of the records returned, such as plain text, HMTL or XML.

    # See the link below for the possible values of ret_type and ret_mode:
    #     http://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly  # noqa

    sleep(DELAY)

    return return_value


def esummary(data, parser):  # noqa

    eutil = 'esummary.fcgi'

    db = data['Database']
    query_key = data['QueryKey']
    web_env = data['WebEnv']
    count = data['Count']

    retmax = 500
    retstart = 0

    return_value = []

    for retstart in range(0, count, retmax):

        url = ENTREZ_BASE_URL + eutil

        params = {'db': db, 'query_key': query_key, 'WebEnv': web_env,
                  'retstart': str(retstart), 'retmax': str(retmax),
                  'rettype': 'docsum'}

        response = get(url, params, 'xml')

        parsed = parser(response.text)

        for item in parsed:
            return_value.append(item)

    sleep(DELAY)

    return return_value


def esearch_epost(term, db):  # noqa
    if type(term) in [list, tuple]:
        term = ' OR '.join(term)
    esearch_results = esearch(db, term)
    id_list = esearch_results['IdList']
    epost_results = epost(db, id_list)
    return epost_results


def summary(term, db):  # noqa
    epost_results = esearch_epost(term, db)
    esummary_results = esummary(epost_results, parse_esummary_xml_text)
    return esummary_results


def dnld_seqs(term, db):  # noqa
    epost_results = esearch_epost(term, db)
    efetch_results = efetch(epost_results, parse_gbseq_xml_text, 'gb')
    return efetch_results


def sra_info(term):  # noqa
    epost_results = esearch_epost(term, 'sra')
    efetch_results = efetch(epost_results, parse_efetch_sra_xml_text,
                            'runinfo')
    return efetch_results
