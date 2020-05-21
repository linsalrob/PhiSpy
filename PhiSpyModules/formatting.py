"""
This is a helper module to just make colored text if we are *not* redirected

Colors come from https://stackoverflow.com/questions/287871/print-in-terminal-with-colors
Testing redirect comes from https://stackoverflow.com/questions/1512457/determining-if-stdout-for-a-python-process-is-redirected
"""

import os
import sys

from .errors import ColorNotFoundError


"""
if os.fstat(0) == os.fstat(1):
    we are not redirected
"""

"""
Colors that you can import and make the text look pretty


"""

__author__ = 'Rob Edwards'


class colors(object):
    color = {
        'HEADER': '\033[95m',
        'OKBLUE': '\033[94m',
        'OKGREEN': '\033[92m',
        'WARNING': '\033[93m',
        'FAIL': '\033[91m',
        'ENDC': '\033[0m',
        'BOLD': '\033[1m',
        'UNDERLINE': '\033[4m',
        'PINK': '\033[95m',
        'BLUE': '\033[94m',
        'GREEN': '\033[92m',
        'YELLOW': '\033[93m',
        'RED': '\033[91m',
        'WHITE': '\033[0m',
        }


def message(msg, c, stream):
    """
    Print a message to stderr using color
    :param msg: the message to print
    :param c: the colour to use
    :param stream: either stderr or stdout
    :return: nothing
    """

    c = c.upper()

    if c not in colors.color:
        raise ColorNotFoundError(f"Color {c} was not found")

    if stream.lower() == 'stderr':
        if os.fstat(0) == os.fstat(1):
            #  stderr is not redirected
            sys.stderr.write(f"{colors.color[c]}{msg}{colors.color['ENDC']}\n")
        else:
            sys.stderr.write(f"{msg}\n")
    elif stream.lower() == 'stdout':
        if os.fstat(0) == os.fstat(1):
            #  stderr is not redirected
            sys.stdout.write(f"{colors.color[c]}{msg}{colors.color['ENDC']}\n")
        else:
            sys.stdout.write(f"{msg}\n")
    else:
        raise IOError(f"There is no IO stream {stream}")