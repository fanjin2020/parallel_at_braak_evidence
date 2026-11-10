"""Runtime checks for the Linux analysis entry points."""
import os
import platform
import shutil
import sys

def require_linux():
    if platform.system() != 'Linux':
        raise RuntimeError('Run this workflow on Linux, not Windows or macOS.')

def require_analysis_python():
    require_linux()
    if sys.version_info[:2] != (3,12):
        raise RuntimeError('Use Python 3.12 in the intended analysis environment.')

def resolve_rscript(value):
    resolved=shutil.which(str(value))
    if not resolved or not os.access(resolved,os.X_OK):
        raise FileNotFoundError('Rscript is not executable; activate the R 4.4.2 environment or pass --rscript.')
    return resolved
