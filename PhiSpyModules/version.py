import pkg_resources
# this comes from setuptools

try:
    __version__ = pkg_resources.get_distribution('PhiSpy').version
except Exception:
    __version__ = 'unknown'
