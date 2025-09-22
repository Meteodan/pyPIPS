
try:
    from importlib.metadata import version
except ImportError:
    # Python < 3.8
    from importlib_metadata import version

try:
    __version__ = version("pypips")
except Exception:
    __version__ = "unknown"
