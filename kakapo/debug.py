# -*- coding: utf-8 -*-

"""
Debugging functions, etc.
"""


def debug_print(msg=''): # noqa
    from kakapo.config import DEBUG_MODE
    if DEBUG_MODE:
        import pprint
        PP = pprint.PrettyPrinter(indent=1, width=110, compact=True)
        PP.pprint(msg)
    else:
        pass
