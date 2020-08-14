"""EMBL-EBI Search."""

import csv

from kakapo.utils.http import get

# EBI Search RESTful Web Services https://www.ebi.ac.uk/ebisearch/swagger.ebi

EBI_URL = 'https://www.ebi.ac.uk'
EBI_SEARCH_URL = EBI_URL + '/ebisearch/ws/rest'
PFAM_ENTRIES_URL = EBI_SEARCH_URL + '/PFAM_ENTRIES'
PFAM_SEQS_URL = EBI_SEARCH_URL + '/PFAM_SEQS'

EBI_SEARCH_TOTAL_RECORDS_HEADER_KEY = 'X-EBI-Search-Total-Results'
EBI_SEARCH_TOTAL_RECORDS_KEY = 'hitCount'
EBI_SEARCH_RECORD_OFFSET_KEY = 'start'
EBI_SEARCH_PAGE_SIZE_KEY = 'size'


def _ebi_search_parse_csv_text(text):
    text_split = text.split('\n')
    text_split.pop()
    parsed_csv = csv.DictReader(text_split)
    parsed_csv = {'entries': list(parsed_csv)}
    return parsed_csv


def _ebi_search_paged_get(url, params, page_size, data_key,
                          response_format='json',
                          record_count_in_header=False):

    params[EBI_SEARCH_RECORD_OFFSET_KEY] = 0
    params[EBI_SEARCH_PAGE_SIZE_KEY] = page_size

    response = get(url, params, response_format)

    res_parsed = None

    if response_format == 'csv':
        res_parsed = _ebi_search_parse_csv_text(response.text)

    elif response_format == 'json':
        res_parsed = response.json()

    data = res_parsed[data_key]

    if record_count_in_header:
        total_records = response.headers[EBI_SEARCH_TOTAL_RECORDS_HEADER_KEY]
        total_records = int(total_records)
    else:
        total_records = res_parsed[EBI_SEARCH_TOTAL_RECORDS_KEY]

    if total_records > page_size:
        pages = total_records // page_size
        for p in range(1, pages + 1):
            params[EBI_SEARCH_RECORD_OFFSET_KEY] = p * page_size
            response = get(url, params, response_format)
            if response_format == 'csv':
                res_parsed = _ebi_search_parse_csv_text(response.text)
            elif response_format == 'json':
                res_parsed = response.json()
            data = data + res_parsed[data_key]

    return data


def pfam_entry(query, fields=()):

    params = {'query': query,
              'fields': ','.join(fields)}

    response = _ebi_search_paged_get(url=PFAM_ENTRIES_URL, params=params,
                                     page_size=100, data_key='entries',
                                     response_format='json')

    return response


def pfam_seqs(query, fields=('NCBI_TAXID', 'UNIPROTKB')):

    params = {'query': query,
              'fields': ','.join(fields)}

    response = _ebi_search_paged_get(url=PFAM_SEQS_URL, params=params,
                                     page_size=100, data_key='entries',
                                     response_format='csv',
                                     record_count_in_header=True)

    return response


def pfam_tax_prot_ids(pfam_seqs_list):
    txids = [int(x['NCBI_TAXID']) for x in pfam_seqs_list]
    prids = [x['UNIPROTKB'] for x in pfam_seqs_list]
    return {'tax_ids': txids, 'uniprot_ids': prids}


def prot_ids_for_tax_ids(pfam_seqs_list, tax_ids):
    if type(tax_ids) not in (list, tuple, set):
        tax_ids = [int(tax_ids), ]

    uniprot_data = pfam_tax_prot_ids(pfam_seqs_list)
    uniprot_tax_ids = uniprot_data['tax_ids']
    uniprot_prot_ids = uniprot_data['uniprot_ids']

    idx = [i for i, value in enumerate(uniprot_tax_ids) if value in tax_ids]
    prot_ids = [uniprot_prot_ids[i] for i in idx]

    return prot_ids
