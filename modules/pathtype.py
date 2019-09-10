# ==============
# Type classes
# PathType
# This argparse Type allows for verification the passed argument is a path, but
# not opening a file immediately.
#
# This can be useful for multiple things.
# - 1) In case of an unrelated error causes the script to crash or hang
#      the file will not be left open, as it was never opened.
#      exist without having ever created an output file.
# - 2) Can be used to specify
#
# ==============
# Example Usage
# Check for a valid and existing file.

# parser = argparse.ArgumentParser()
#
# # Add an argument that must be an existing file, but can also be specified as a dash ('-') in the command.
# parser.add_argument('existing_file', type = PathType(type='file', dash_ok=True, exists = True))
#
# # Add an argument for a folder that must NOT exist. Note: folders can not be "dash_ok".
# parser.add_argument('non_existant_folder', type = PathType(type='dir', exists = False))
#
# # Add an argument for EITHER a folder or a file, but can be dash_ok, and don't check if it exists.
# parser.add_argument('maybe_existant_file_or_folder',type=PathType(type=('dir','file').__contains__,dash_ok=True,exists=None))

# ==============

# Code taken from https://stackoverflow.com/questions/11415570/directory-path-types-with-argparse
# This is EXACTLY what I was going to build otherwise.
#
# This is also suggested to be an official change to argparse.
# https://mail.python.org/pipermail/stdlib-sig/2015-July/000990.html

import os
from argparse import ArgumentTypeError

# Existance Options
EXIST_CHECK  = True  # File/folder must exist.
EXIST_INVERT = False # Must *not* exist.
EXIST_EITHER = None  # Allowed to reguardless.

# Types allowed for files.
TYPE_FILE = 'file'
TYPE_DIR  = 'dir'
TYPE_SYM  = 'symlink'
TYPE_ALL  = None # Any thing is fine (don't check)

# Whether folder is not empty
EMPTY_CHECK  = True
EMPTY_INVERT = False
EMPTY_EITHER = None

class PathType(object):
    def __init__(self, exists = EXIST_CHECK, type = TYPE_FILE, empty = EMPTY_CHECK, dash_ok = True):
        """
           exists:
                True: a path that does exist
                False: a path that does not exist, in a valid parent directory
                None: don't care
           type: 'file', 'dir', 'symlink', None, a list of these,
                 or a function returning True for valid paths
                 If None, the type of file/directory/symlink is not checked.
           dash_ok: whether to allow "-" as stdin/stdout
        """
        assert exists in (EXIST_CHECK, EXIST_INVERT, EXIST_EITHER)
        assert empty  in (EMPTY_CHECK, EMPTY_INVERT, EMPTY_EITHER)
        # Make sure type is file, dir, sym, None, list, or callable.
        if isinstance(type, (list, tuple)):
            # Type is a list, make sure that it includes only "file", "dir", or "sym"
            for t in type:
                assert t in (TYPE_FILE, TYPE_DIR, TYPE_SYM)
            # Type is the contains check for this lists.
            type = type.__contains__
        elif not callable(type):
            # Otherwise, make sure this is valid.
            assert type in (TYPE_FILE, TYPE_DIR, TYPE_SYM, EXIST_EITHER)
        # else; type is a callable object, and this is ok.
        self._exists = exists
        self._type = type
        self._empty = empty
        self._dash_ok = dash_ok

    def __call__(self, string):
        if string == '-':
            # the special argument "-" means sys.std[in/out]
            if self._type == TYPE_DIR:
                raise ArgumentTypeError('standard input/output (-) not allowed as directory path')
            elif self._type == TYPE_SYM:
                raise ArgumentTypeError('standard input/output (-) not allowed as symlink path')
            elif not self._dash_ok:
                raise ArgumentTypeError('standard input/output (-) not allowed')
            return string # No reason to check anything else if this works.
        # If the file must exist.
        if self._exists == EXIST_CHECK:
            if not os.path.exists(string):
                raise ArgumentTypeError("path does not exist: %r" % string)
            # If must be file and not: Fail
            if self._type == TYPE_FILE:
                if not os.path.isfile(string):
                    raise ArgumentTypeError("path is not a file: %r" % string)
            # If must be dir and not: Fail
            elif self._type == TYPE_DIR:
                if not os.path.isdir(string):
                    raise ArgumentTypeError("path is not a directory: %r" % string)
            # If must be sym and not: Fail
            elif self._type == TYPE_SYM:
                if not os.path.islink(string):
                    raise ArgumentTypeError("path is not a symlink: %r" % string)
            # Otherwise, call type.
            elif not self._type(string):
                raise ArgumentTypeError("path not valid: %r" % string)
        else:
            if self._exists == False and os.path.exists(string):
                raise ArgumentTypeError("path exists: %r" % string)
            p = os.path.dirname(string) or os.path.curdir
            if not os.path.exists(p):
                raise ArgumentTypeError("parent directory does not exist: %r" % p)
            elif not os.path.isdir(p):
                raise ArgumentTypeError("parent path is not a directory: %r" % p)
        # If the file must be empty
        if self._empty == EMPTY_CHECK:
            if os.path.exists(string) and os.listdir(string):
                raise ArgumentTypeError("folder is not empty: %r" % string)
        return string
