"""EMBL-EBI Proteins."""

from math import ceil

from kakapo.utils.http import get

# Proteins REST API https://www.ebi.ac.uk/proteins/api/doc -------------------

EBI_URL = 'https://www.ebi.ac.uk'

PROT_API_URL = EBI_URL + '/proteins/api'


def fasta_by_accession_list(acc_list):
    url = PROT_API_URL + '/proteins'

    max_recs = 90
    pages = int(ceil(len(acc_list) / max_recs))

    data = ''

    for p in range(pages):
        offset = p * max_recs
        last = offset + max_recs
        acc_param = ','.join(x for x in acc_list[offset:last])
        params = {'accession': acc_param}
        response = get(url=url, params=params, response_format='fasta')
        txt = ''
        if response is not None:
            txt = response.text
        data += txt

    return data
