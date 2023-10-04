"""Kakapo logging class."""

from datetime import datetime
from io import IOBase, StringIO
from os.path import abspath, expanduser
from typing import Union

from click import echo, style

HANDLE_TYPES = (IOBase, StringIO)


class Log:
    _console = True
    _write = False
    _file: Union[str, None] = None
    _handle: Union[IOBase, StringIO, None] = None
    _colors = False
    _time_stamp = True

    @classmethod
    def colors(cls):
        return cls._colors

    @classmethod
    def set_colors(cls, value: bool):
        assert type(value) == bool
        cls._colors = value

    @classmethod
    def console(cls):
        return cls._console

    @classmethod
    def set_console(cls, value: bool):
        assert type(value) == bool
        cls._console = value

    @classmethod
    def write(cls):
        return cls._write

    @classmethod
    def set_write(cls, value: bool):
        assert type(value) == bool
        cls._write = value

    @classmethod
    def file(cls):
        return cls._file

    @classmethod
    def set_file(cls, value: Union[str, IOBase, StringIO]):
        if type(value) == str:
            cls._file = abspath(expanduser(value))
            cls._handle = None
        elif isinstance(value, HANDLE_TYPES):
            cls._handle = value
            cls._file = None

    @classmethod
    def time_stamp(cls):
        return cls._time_stamp

    @classmethod
    def set_time_stamp(cls, value: bool):
        assert type(value) == bool
        cls._time_stamp = value

    @classmethod
    def ts(cls) -> str:
        if cls.time_stamp() is True:
            now = datetime.now()
            return now.strftime('%Y-%m-%d %H:%M:%S') + ' '
        else:
            return ''

    @classmethod
    def log(cls,
            inf: str = '',
            msg: str = '',
            wrn: str = '',
            err: str = '',
            s: str = '',
            timestamp: bool = True) -> str:

        sts = cls.ts()
        lts = len(sts)

        if timestamp is True:
            ts = sts
        else:
            ts = ' ' * lts

        if inf != '':
            inf = str(inf)
        if msg != '':
            msg = str(msg)
        if wrn != '':
            wrn = str(wrn)
        if err != '':
            err = str(err)
        if s != '':
            s = ' ' + str(s)

        ts_color = 'reset'
        inf_color = 'reset'
        msg_color = 'reset'
        wrn_color = 'reset'
        err_color = 'reset'

        if cls.colors() is True:
            ts_color = 'bright_yellow'
            inf_color = 'bright_blue'
            msg_color = 'green'
            wrn_color = 'yellow'
            err_color = 'red'

        ts_style = style(ts, fg=ts_color)

        s_style = style(s, fg='reset', reset=False)
        txt = '{ts}{inf}{msg}{wrn}{err}{s}'

        inf_style = style(inf, fg=inf_color, bold=True)
        msg_style = style(msg, fg=msg_color)
        wrn_style = style(wrn, fg=wrn_color)
        err_style = style(err, fg=err_color)

        msg_c = txt.format(ts=ts_style, inf=inf_style, msg=msg_style,
                           wrn=wrn_style, err=err_style, s=s_style)

        msg_c = msg_c.replace('\n', '\n' + ' ' * lts)

        msg_f = txt.format(ts=ts, inf=inf, msg=msg, wrn=wrn, err=err, s=s)

        msg_f = msg_f.replace('\n', '\n' + ' ' * lts) + '\n'

        if cls._write is True:

            handle = False
            if cls._handle is not None:
                handle = True
                f = cls._handle
            else:
                assert cls._file is not None
                f = open(cls._file, 'a+')

            f.write(msg_f)

            if handle is False:
                f.close()

        if cls._console is True:
            echo(msg_c)

        return msg_f

    @classmethod
    def inf(cls, i: str) -> str:
        return cls.log(inf=i, timestamp=True)

    @classmethod
    def msg(cls, m: str, s: str) -> str:
        return cls.log(msg=m, s=s, timestamp=True)

    @classmethod
    def wrn(cls, w: str, s: str) -> str:
        return cls.log(wrn=w, s=s, timestamp=True)

    @classmethod
    def err(cls, e: str, s: str) -> str:
        return cls.log(err=e, s=s, timestamp=True)
