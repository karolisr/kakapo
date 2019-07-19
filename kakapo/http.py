# -*- coding: utf-8 -*-

"""http"""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import with_statement

import requests

# Possible values for the Accept request-header field:
ACC_HEAD = {'json': {'Accept': 'application/json'},
            'csv': {'Accept': 'text/csv'},
            'fasta': {'Accept': 'text/x-fasta'},
            'plain_text': {'Accept': 'text/plain'},
            'xml': {'Accept': 'application/xml'}}


def _valid_response_formats():
    return tuple(ACC_HEAD.keys())


def get(url, params, response_format):
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

    response = requests.get(url=url, params=params, headers=headers)

    if response.ok:
        pass
    else:
        return None

    return response


def post(url, data, response_format):
    """
    Wrap requests.post
    """
    assert response_format in _valid_response_formats()
    headers = ACC_HEAD[response_format]
    response = requests.post(url=url, data=data, headers=headers)
    return response
