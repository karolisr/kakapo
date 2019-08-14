# -*- coding: utf-8 -*-

"""Logging"""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import with_statement

import logging
import sys

from kakapo.config import CONYELL, CONSDFL


def prepare_logger(console=True, stream=None, file=None):  # noqa

    handlers = logging.getLogger().handlers

    for h in handlers:
        logging.getLogger().removeHandler(h)

    format_console = CONYELL + '%(asctime)s' + CONSDFL + ' - %(message)s'
    format_file = '%(asctime)s - %(message)s'
    date_time_format = '%Y/%m/%d %H:%M:%S'

    formatter_console = logging.Formatter(fmt=format_console,
                                          datefmt=date_time_format)

    formatter_file = logging.Formatter(fmt=format_file,
                                       datefmt=date_time_format)

    console_log_handler = None
    stream_log_handler = None
    file_log_handler = None

    if console is True:
        console_log_handler = logging.StreamHandler(sys.stdout)
        console_log_handler.setFormatter(formatter_console)
        logging.getLogger().addHandler(console_log_handler)

    if stream is not None:
        stream_log_handler = logging.StreamHandler(stream)
        stream_log_handler.setFormatter(formatter_file)
        logging.getLogger().addHandler(stream_log_handler)

    if file is not None:
        file_log_handler = logging.FileHandler(file)
        file_log_handler.setFormatter(formatter_file)
        logging.getLogger().addHandler(file_log_handler)

    logging.getLogger().setLevel(logging.INFO)

    return logging, stream_log_handler
