# -*- coding: utf-8 -*-

"""http"""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import with_statement

import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
from requests.exceptions import HTTPError

# Possible values for the Accept request-header field:
ACC_HEAD = {'json': {'Accept': 'application/json'},
            'csv': {'Accept': 'text/csv'},
            'fasta': {'Accept': 'text/x-fasta'},
            'plain_text': {'Accept': 'text/plain'},
            'xml': {'Accept': 'application/xml'}}


def _valid_response_formats():
    return tuple(ACC_HEAD.keys())


def retry_session(retries=5, backoff_factor=1,
                  status_forcelist=(500, 502, 504)):  # noqa

    session = requests.Session()

    retry = Retry(
        total=retries,
        read=retries,
        connect=retries,
        backoff_factor=backoff_factor,
        status_forcelist=status_forcelist)

    adapter = HTTPAdapter(max_retries=retry)

    session.mount('http://', adapter)
    session.mount('https://', adapter)

    return session


def get(url, params=None, response_format='json'):
    """
    Wrap requests.get

    :param url: url
    :type url: str

    :param params: parameter names and values
    :type params: dict

    :param response_format: 'json', etc.
    :type response_format: str

    :returns: parsed response
    :rtype: dict or list
    """
    if type(response_format) in (dict, ):
        headers = response_format
    else:
        assert response_format in _valid_response_formats()
        headers = ACC_HEAD[response_format]

    with retry_session() as session:
        response = session.get(url=url, params=params, headers=headers)

    try:
        response.raise_for_status()
    except HTTPError as e:
        pass
        # print(e)

    return response


def post(url, data, response_format):
    """
    Wrap requests.post
    """
    assert response_format in _valid_response_formats()
    headers = ACC_HEAD[response_format]

    with retry_session() as session:
        response = session.post(url=url, data=data, headers=headers)

    try:
        response.raise_for_status()
    except HTTPError as e:
        pass
        # print(e)

    return response


def download_file(url, local_path):  # noqa
    r = get(url)

    with open(local_path, 'wb') as f:
        f.write(r.content)
