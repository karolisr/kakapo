# -*- coding: utf-8 -*-

"""http."""

from sys import exit
from typing import Union

from requests import Response, Session
from requests.adapters import HTTPAdapter
from requests.exceptions import HTTPError
from urllib3.util import Retry

# Possible values for the Accept request-header field:
ACC_HEAD: dict[str, dict[str, str]] = {
    'csv': {'Accept': 'text/csv'},
    'fasta': {'Accept': 'text/x-fasta'},
    'json': {'Accept': 'application/json'},
    'text': {'Accept': 'text/plain'},
    'asn.1': {'Accept': 'text/plain'},
    'xml': {'Accept': 'application/xml'}}


def _valid_response_formats():
    return tuple(ACC_HEAD.keys())


def retry_session(retries=5, backoff_factor=1,
                  status_forcelist=(500, 502, 503, 504)):

    session = Session()

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


def get(url, params=None, response_format: Union[str, dict] = 'json'
        ) -> Response:
    """
    Wrap requests.get.

    :param url: url
    :type url: str

    :param params: parameter names and values
    :type params: dict

    :param response_format: 'json', etc.
    :type response_format: str

    :returns: Response
    :rtype: Response
    """

    headers: dict[str, str]

    if type(response_format) == dict:
        headers = response_format
    else:
        assert response_format in _valid_response_formats()
        headers = ACC_HEAD[response_format]

    with retry_session() as session:
        r: Response = session.get(
            url=url, params=params, headers=headers)

    try:
        r.raise_for_status()
    except HTTPError as e:
        print(e)
        exit(1)

    return r


def post(url, data, response_format) -> Response:
    """Wrap requests.post."""
    assert response_format in _valid_response_formats()
    headers = ACC_HEAD[response_format]

    with retry_session() as session:
        r: Response = session.post(url=url, data=data, headers=headers)

    try:
        r.raise_for_status()
    except HTTPError as e:
        print(e)
        exit(1)

    return r


def download_file(url, local_path, response_format: Union[str, dict] = 'json'):
    r: Response = get(url, response_format=response_format)

    with open(local_path, 'wb') as f:
        f.write(r.content)
