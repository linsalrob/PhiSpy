"""
A custom error class for PhiSpy to make the code more accessible
"""
class Error(Exception):
    """
    Base class for exceptions in this module.
    """
    pass

class ColorNotFoundError(Error):
    """
    Exception raised for sequences not being paired properly.

    :param message: explanation of the error
    """

    def __init__(self, message):
        self.message = message
