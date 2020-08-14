"""Kakapo logging class."""

import io

from datetime import datetime
from io import StringIO
from os.path import abspath
from os.path import expanduser

from click import echo
from click import style


HANDLE_TYPES = (io.IOBase, StringIO)


class Log:
    _console = True
    _write = False
    _file = None
    _colors = False

    @classmethod
    def colors(cls):
        return cls._colors

    @classmethod
    def set_colors(cls, value):
        assert type(value) == bool
        cls._colors = value

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
        if type(value) == str:
            cls._file = abspath(expanduser(value))
        elif isinstance(value, HANDLE_TYPES):
            cls._file = value

    @staticmethod
    def time_stamp():
        now = datetime.now()
        return now.strftime('%Y-%m-%d %H:%M:%S')

    @classmethod
    def log(cls, msg_inf='', msg='', wrn='', err='', s=''):

        ts = cls.time_stamp()

        if msg_inf != '':
            msg_inf = str(msg_inf)
        if msg != '':
            msg = str(msg)
        if wrn != '':
            wrn = str(wrn)
        if err != '':
            err = str(err)
        if s != '':
            s = str(s)

        ts_color = 'reset'
        inf_color = 'reset'
        msg_inf_color = inf_color
        msg_color = 'reset'
        wrn_color = 'reset'
        err_color = 'reset'

        if cls.colors() is True:
            ts_color = 'bright_yellow'
            inf_color = 'bright_blue'
            msg_inf_color = inf_color
            msg_color = 'green'
            wrn_color = 'yellow'
            err_color = 'red'

        ts_style = style(ts, fg=ts_color, bold=False)

        s_style = style(s, bold=False, underline=False)
        txt = '{ts} {msg_inf}{msg}{wrn}{err} {s}'
        if msg_inf + msg + wrn + err == '':
            s_style = style(s, fg=inf_color, bold=False, underline=False)
            txt = '{ts} {msg_inf}{msg}{wrn}{err}{s}'

        msg_inf_style = style(msg_inf, fg=msg_inf_color, bold=False)
        msg_style = style(msg, fg=msg_color, bold=False)
        wrn_style = style(wrn, fg=wrn_color, bold=False)
        err_style = style(err, fg=err_color, bold=False)

        msg_c = txt.format(ts=ts_style, msg_inf=msg_inf_style, msg=msg_style,
                           wrn=wrn_style, err=err_style, s=s_style)

        msg_c = msg_c.replace('\n', '\n' + ' ' * (len(ts) + 1))

        msg_f = txt.format(ts=ts, msg_inf=msg_inf, msg=msg, wrn=wrn, err=err,
                           s=s)

        msg_f = msg_f.replace('\n', '\n' + ' ' * (len(ts) + 1)) + '\n'

        if cls._file is not None and cls._write is True:

            handle = False
            if isinstance(cls._file, HANDLE_TYPES):
                handle = True
                f = cls._file
            if handle is False:
                f = open(cls._file, 'a+')

            f.write(msg_f)

            if handle is False:
                f.close()

        if cls._console is True:
            echo(msg_c)

        return msg_f

    @classmethod
    def inf(cls, s=''):
        return cls.log(s=s)

    @classmethod
    def msg_inf(cls, m, s=''):
        return cls.log(msg_inf=m, s=s)

    @classmethod
    def msg(cls, m, s=''):
        return cls.log(msg=m, s=s)

    @classmethod
    def wrn(cls, w, s=''):
        return cls.log(wrn=w, s=s)

    @classmethod
    def err(cls, e, s=''):
        return cls.log(err=e, s=s)
