"""Run a package R entry point on Linux without a shell or temporary wrapper."""
import argparse
from pathlib import Path
import subprocess
from linux_runtime import require_linux,resolve_rscript

def run(rscript,script,arguments):
    source=Path(__file__).resolve().parent
    if Path(script).name!=script or not (source/script).is_file():
        raise ValueError('R entry point must be a script in this package')
    if any('\n' in str(x) or '\r' in str(x) for x in arguments):
        raise ValueError('Newlines in command arguments are not supported')
    require_linux()
    return subprocess.call([resolve_rscript(rscript),'--vanilla',str(source/script),*[str(x) for x in arguments]])

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--rscript',default='Rscript');p.add_argument('--script',required=True)
    p.add_argument('arguments',nargs=argparse.REMAINDER);a=p.parse_args()
    raise SystemExit(run(a.rscript,a.script,a.arguments[1:] if a.arguments[:1]==['--'] else a.arguments))
