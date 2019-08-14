# -*- coding: utf-8 -*-

"""
EMBL-EBI Proteins
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import with_statement

from math import ceil

from kakapo.http_k import get

# Proteins REST API https://www.ebi.ac.uk/proteins/api/doc ------------------

EBI_URL = 'https://www.ebi.ac.uk'

PROT_API_URL = EBI_URL + '/proteins/api'


def fasta_by_accession_list(acc_list): # noqa
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
        if response is None:
            response = ''
        data = data + response.text

    return data
