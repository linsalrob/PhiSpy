from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version('PhiSpy')
except PackageNotFoundError:
    __version__ = "unknown"
except Exception:
    __version__ = 'unknown'
