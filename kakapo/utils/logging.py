"""Kakapo logging class."""

import io

from datetime import datetime
from io import StringIO
from os.path import abspath
from os.path import expanduser

from multipledispatch import dispatch
from click import echo
from click import style


LOG_NS = dict()
HANDLE_TYPES = (io.IOBase, StringIO)


class Log:
    _console = True
    _write = False
    _file = None
    _colors = False
    _time_stamp = True

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

    @classmethod
    def time_stamp(cls):
        return cls._time_stamp

    @classmethod
    def set_time_stamp(cls, value):
        assert type(value) == bool
        cls._time_stamp = value

    @classmethod
    def ts(cls):
        if cls.time_stamp() is True:
            now = datetime.now()
            return now.strftime('%Y-%m-%d %H:%M:%S') + ' '
        else:
            return ''

    @classmethod
    def log(cls, inf='', msg='', wrn='', err='', s='', ts=True):

        sts = cls.ts()
        lts = len(sts)

        if ts is True:
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

    # inf --------------------------------------------------------------------
    @classmethod
    @dispatch(object, str, namespace=LOG_NS)
    def inf(cls, i):
        return cls.log(inf=i, ts=True)

    @classmethod
    @dispatch(object, str, bool, namespace=LOG_NS)
    def inf(cls, i, ts):
        return cls.log(inf=i, ts=ts)

    @classmethod
    @dispatch(object, str, str, namespace=LOG_NS)
    def inf(cls, i, s):
        return cls.log(inf=i, s=s, ts=True)

    @classmethod
    @dispatch(object, str, str, bool, namespace=LOG_NS)
    def inf(cls, i, s, ts):
        return cls.log(inf=i, s=s, ts=ts)

    # msg --------------------------------------------------------------------
    @classmethod
    @dispatch(object, str, namespace=LOG_NS)
    def msg(cls, m):
        return cls.log(msg=m, ts=True)

    @classmethod
    @dispatch(object, str, bool, namespace=LOG_NS)
    def msg(cls, m, ts):
        return cls.log(msg=m, ts=ts)

    @classmethod
    @dispatch(object, str, str, namespace=LOG_NS)
    def msg(cls, m, s):
        return cls.log(msg=m, s=s, ts=True)

    @classmethod
    @dispatch(object, str, str, bool, namespace=LOG_NS)
    def msg(cls, m, s, ts):
        return cls.log(msg=m, s=s, ts=ts)

    # wrn --------------------------------------------------------------------
    @classmethod
    @dispatch(object, str, namespace=LOG_NS)
    def wrn(cls, w):
        return cls.log(wrn=w, ts=True)

    @classmethod
    @dispatch(object, str, bool, namespace=LOG_NS)
    def wrn(cls, w, ts):
        return cls.log(wrn=w, ts=ts)

    @classmethod
    @dispatch(object, str, str, namespace=LOG_NS)
    def wrn(cls, w, s):
        return cls.log(wrn=w, s=s, ts=True)

    @classmethod
    @dispatch(object, str, str, bool, namespace=LOG_NS)
    def wrn(cls, w, s, ts):
        return cls.log(wrn=w, s=s, ts=ts)

    # err --------------------------------------------------------------------
    @classmethod
    @dispatch(object, str, namespace=LOG_NS)
    def err(cls, e):
        return cls.log(err=e, ts=True)

    @classmethod
    @dispatch(object, str, bool, namespace=LOG_NS)
    def err(cls, e, ts):
        return cls.log(err=e, ts=ts)

    @classmethod
    @dispatch(object, str, str, namespace=LOG_NS)
    def err(cls, e, s):
        return cls.log(err=e, s=s, ts=True)

    @classmethod
    @dispatch(object, str, str, bool, namespace=LOG_NS)
    def err(cls, e, s, ts):
        return cls.log(err=e, s=s, ts=ts)


# Test Log:
if __name__ == '__main__':
    log = Log()
    log.set_colors(True)

    # Time-stamp will be printed by default:
    log.inf('inf', 'inf')
    log.inf('inf')
    log.msg('msg', 'msg')
    log.msg('msg')
    log.wrn('wrn', 'wrn')
    log.wrn('wrn')
    log.err('err', 'err')
    log.err('err')

    # Time-stamp will not be printed:
    log.inf('inf', 'inf', False)
    log.inf('inf', False)
    log.msg('msg', 'msg', False)
    log.msg('msg', False)
    log.wrn('wrn', 'wrn', False)
    log.wrn('wrn', False)
    log.err('err', 'err', False)
    log.err('err', False)

    log.msg('msg1', '\nmsg2')
    log.msg('msg1\nmsg2')

    # Time-stamp will NEVER be printed:
    log.set_time_stamp(False)

    log.inf('inf', 'inf')
    log.inf('inf')
    log.msg('msg', 'msg')
    log.msg('msg')
    log.wrn('wrn', 'wrn')
    log.wrn('wrn')
    log.err('err', 'err')
    log.err('err')

    log.inf('inf', 'inf', False)
    log.inf('inf', False)
    log.msg('msg', 'msg', False)
    log.msg('msg', False)
    log.wrn('wrn', 'wrn', False)
    log.wrn('wrn', False)
    log.err('err', 'err', False)
    log.err('err', False)

    log.msg('msg1', '\nmsg2')
    log.msg('msg1\nmsg2')
