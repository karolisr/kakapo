# -*- coding: utf-8 -*-

"""http"""

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

    # print('*' * 100)
    # print(response.headers)
    # print(response.status_code, response.reason)

    if response.ok:
        pass
        # if type(response_format) in (dict, ):
        #     response_parsed = response
        # else:
        #     if response_format == 'json':
        #         response_parsed = response.json()
        #     elif response_format == 'csv':
        #         response_parsed = response.text
        #     elif response_format == 'fasta':
        #         response_parsed = response.text
        #     elif response_format == 'plain_text':
        #         response_parsed = response.text
        #     elif response_format == 'xml':
        #         response_parsed = response.text
        #     else:
        #         # response_parsed = response.text
        #         raise NotImplementedError

    else:
        pass
        # response_parsed = None

    return response


def post(url, data, response_format):
    """
    Wrap requests.post
    """
    assert response_format in _valid_response_formats()
    headers = ACC_HEAD[response_format]
    response = requests.post(url=url, data=data, headers=headers)

    print('*' * 100)
    print(response.headers)
    print(response.status_code, response.reason)

    return response
