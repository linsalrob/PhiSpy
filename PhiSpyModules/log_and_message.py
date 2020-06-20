"""
This is  a critical component but it causes circular dependencies
because everything needs it, so we abstract it in its own class
"""

import os
import sys
import logging
from .formatting import Colors, ColorNotFoundError


def message(msg, c, stream):
    """
    Print a message to stderr using color
    :param msg: the message to print
    :param c: the colour to use
    :param stream: either stderr or stdout
    :return: nothing
    """

    c = c.upper()
    # we strip any newline off message in case it was added and then we add it here! Ensures always only one newline
    msg = msg.strip()
    if c not in Colors.color:
        raise ColorNotFoundError(f"Color {c} was not found")

    if stream.lower() == 'stderr':
        if os.fstat(0) == os.fstat(1):
            #  stderr is not redirected
            sys.stderr.write(f"{Colors.color[c]}{msg}{Colors.color['ENDC']}\n")
        else:
            sys.stderr.write(f"{msg}\n")
    elif stream.lower() == 'stdout':
        if os.fstat(0) == os.fstat(1):
            #  stderr is not redirected
            sys.stdout.write(f"{Colors.color[c]}{msg}{Colors.color['ENDC']}\n")
        else:
            sys.stdout.write(f"{msg}\n")
    else:
        raise IOError(f"There is no IO stream {stream}")


def log_and_message(msg, c="WHITE", stderr=False, stdout=False, quiet=False, loglevel="INFO"):
    """
    Write the message to both the log and an output stream. By default we will just log the message in the
    log.

    Note that we also adhere to the quiet option of self, and will only write to the log
    Set either stderr or stdout to true to write to those streams too (but you will need to reset quiet if appropriate)

    :param self: The arg parser object
    :param msg: the message to write
    :param c: the color to write the message
    :param stderr: write the message to stderr
    :param stdout: write the message to stdout
    :param loglevel: the logging level. See https://docs.python.org/3/library/logging.html#levels for a list
    :return:
    """

    if loglevel == 'CRITICAL':
        logging.getLogger('PhiSpy').critical(msg.strip())
    elif loglevel == 'ERROR':
        logging.getLogger('PhiSpy').error(msg.strip())
    elif loglevel == 'WARNING':
        logging.getLogger('PhiSpy').warning(msg.strip())
    elif loglevel == 'DEBUG':
        logging.getLogger('PhiSpy').debug(msg.strip())
    else:
        logging.getLogger('PhiSpy').info(msg.strip())

    if not quiet:
        if stderr:
            message(msg,c, "stderr")
        if stdout:
            message(msg, c, "stdout")

