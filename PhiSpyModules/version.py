try:
    from importlib.metadata import version, PackageNotFoundError
except ImportError:
    # Python < 3.8
    from importlib_metadata import version, PackageNotFoundError

try:
    __version__ = version('PhiSpy')
except PackageNotFoundError:
    __version__ = "unknown"
except Exception:
    __version__ = 'unknown'
