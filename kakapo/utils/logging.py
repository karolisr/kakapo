"""Kakapo logging class."""

import logging
import sys
from datetime import datetime

from click import echo
from click import style


from kakapo.tools.config import CONYELL, CONSDFL


def prepare_logger(console=True, stream=None, file=None):

    handlers = logging.getLogger().handlers

    for h in handlers:
        logging.getLogger().removeHandler(h)

    format_console = CONYELL + '%(asctime)s' + CONSDFL + ' %(message)s'
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


class Log:
    _console = True
    _write = False
    _file = None

    @classmethod
    def console(cls):
        return cls._console

    @classmethod
    def set_console(cls, value):
        assert type(value) == bool
        cls._console = value

    @classmethod
    def write(cls):
        return cls._write

    @classmethod
    def set_write(cls, value):
        assert type(value) == bool
        cls._write = value

    @classmethod
    def file(cls):
        return cls._file

    @classmethod
    def set_file(cls, value):
        assert type(value) == str
        cls._file = value

    @staticmethod
    def time_stamp():
        now = datetime.now()
        return now.strftime('%Y-%m-%d %H:%M:%S')

    @classmethod
    def log(cls, msg='', wrn='', err='', s=''):

        ts = cls.time_stamp()

        if msg != '':
            msg = str(msg)
        if wrn != '':
            wrn = str(wrn)
        if err != '':
            err = str(err)
        if s != '':
            s = str(s)

        ts_style = style(ts, fg='magenta', bold=False)

        s_style = s
        txt = '{ts} {msg}{wrn}{err} {s}'
        if msg + wrn + err == '':
            s_style = style(s, fg='blue', bold=True, underline=False)
            txt = '{ts} {msg}{wrn}{err}{s}'

        msg_style = style(msg, fg='green', bold=False)
        wrn_style = style(wrn, fg='yellow', bold=False)
        err_style = style(err, fg='red', bold=False)

        msg_c = txt.format(ts=ts_style, msg=msg_style, wrn=wrn_style,
                           err=err_style, s=s_style)

        msg_f = txt.format(ts=ts, msg=msg, wrn=wrn, err=err, s=s) + '\n'

        if cls._file is not None and cls._write is True:
            with open(cls._file, 'a+') as f:
                f.write(msg_f)

        if cls._console is True:
            echo(msg_c)

        return msg_f

    @classmethod
    def inf(cls, s=''):
        return cls.log(s=s)

    @classmethod
    def msg(cls, m, s=''):
        return cls.log(msg=m, s=s)

    @classmethod
    def wrn(cls, w, s=''):
        return cls.log(wrn=w, s=s)

    @classmethod
    def err(cls, e, s=''):
        return cls.log(err=e, s=s)


if __name__ == '__main__':

    Log.file = 'log.txt'
    Log.write = False
    Log.console = True

    print(Log.inf('Log.inf'), end='')
    print(Log.msg('Log.msg', 'str'), end='')
    print(Log.msg('Log.msg'), end='')
    print(Log.wrn('Log.wrn', 'str'), end='')
    print(Log.wrn('Log.wrn'), end='')
    print(Log.err('Log.err', 'str'), end='')
    print(Log.err('Log.err'), end='')
