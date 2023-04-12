"""Spawn new processes."""

import re

from subprocess import run as subp_run


def run(cmd, in_txt=None, capture=True, cwd=None, do_not_raise=False,
        text=True):
    out = subp_run(cmd, input=in_txt, capture_output=capture,
                   cwd=cwd, text=text)
    if do_not_raise is False and out.returncode > 0:
        raise Exception(out.stderr)
    return out


def run_then_grep(cmd, regexp, do_not_raise=True):
    out = run(cmd, do_not_raise=do_not_raise)
    o = out.stdout
    e = out.stderr
    s = o
    if e != '':
        s = o + e
    matched = re.findall(regexp, s)
    return matched


def which(executable):
    return run(['which', executable]).stdout.strip()
