# -*- coding: utf-8 -*-
"""
SRA Hacks
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import with_statement

import re

from time import sleep

from kakapo.http_k import get


def confirm_single_ended(sra_run_acc):  # noqa
    sleep(3)
    url = 'https://trace.ncbi.nlm.nih.gov/Traces/sra'
    params = {'run_spot': sra_run_acc, 'page_size': 1}
    response_format = 'plain_text'
    r = get(url, params, response_format)
    read_mask_bio = re.findall(r'read_mask_bio:(\d)', r.text, re.MULTILINE)
    if len(read_mask_bio) != 0:
        read_mask_bio = int(read_mask_bio[0])
    confirm = False
    # This is equal to 1 for Illumina, but may not be true for other platforms.
    if read_mask_bio == 1:
        confirm = True
    return confirm
