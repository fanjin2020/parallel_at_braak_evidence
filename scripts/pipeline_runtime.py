"""Execution, resume state and atomic publication for local analyses."""
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import time


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(data, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    temp.replace(path)


def inventory(root, suffixes=None):
    root = Path(root)
    if root.is_file():
        return {root.name: sha(root)}
    return {p.relative_to(root).as_posix(): sha(p) for p in sorted(root.rglob('*'))
            if p.is_file() and '__pycache__' not in p.parts
            and (suffixes is None or p.suffix.lower() in suffixes)}


def signature(path, value, resume=False):
    path = Path(path)
    if path.exists():
        if not resume:
            raise ValueError('Existing run: specify --resume with the same output directory.')
        if json.loads(path.read_text(encoding='utf-8')) != value:
            raise ValueError('Inputs, code, environment or settings changed. Use a new run directory.')
    else:
        if resume:
            raise ValueError('Cannot resume: run signature is absent.')
        write_json(path, value)


def recover_publication(out):
    """Finish an interrupted final-directory commit without rerunning models."""
    out = Path(out)
    path = out/'COMPLETED.json'
    if path.is_file():
        record = json.loads(path.read_text(encoding='utf-8'))
    elif (out/'audit/PUBLISH.json').is_file():
        record = json.loads((out/'audit/PUBLISH.json').read_text(encoding='utf-8'))
        if not (out/record['output_directory']).is_dir():
            return None
    else:
        return None
    if record['output_directory'] not in {'submission','results_to_review'}:
        raise ValueError('Unrecognized completed output directory')
    if inventory(out/record['output_directory']) != record['outputs']:
        raise ValueError('Completed submission files changed.')
    if not path.exists():write_json(path,record)
    return record


def publish(stage, out, name, report):
    stage, out = Path(stage), Path(out)
    if name not in {'submission','results_to_review'}:
        raise ValueError('Unrecognized publication directory')
    if (out/name).exists():
        raise ValueError('Final output already exists without completion record.')
    record = dict(completed=True,output_directory=name,**report,outputs=inventory(stage))
    write_json(out/'audit/PUBLISH.json',record)
    stage.rename(out/name)
    write_json(out/'COMPLETED.json',record)
    return record


class Runner:
    def __init__(self, audit, cwd, env=None):
        self.audit = Path(audit)
        self.audit.mkdir(parents=True, exist_ok=True)
        self.cwd = Path(cwd)
        self.env = env or os.environ.copy()

    def completed(self, name, command, root):
        path = self.audit / (name + '.json')
        if not path.exists():
            return False
        record = json.loads(path.read_text(encoding='utf-8'))
        if record.get('command') != list(map(str, command)):
            raise ValueError('Step command changed: ' + name)
        if record.get('status') != 'complete':
            return False
        for relative, digest in record['outputs'].items():
            target = Path(root) / relative
            if not target.is_file() or sha(target) != digest:
                raise ValueError('Completed output changed: ' + name + '/' + relative)
        print(name + ': verified checkpoint', flush=True)
        return True

    def run(self, name, command, root, changed_only=False):
        command = list(map(str, command))
        root = Path(root)
        if self.completed(name, command, root):
            return
        # Primary jobs share a working root; record only files created/changed by this job.
        def snapshot():
            return {p.relative_to(root).as_posix(): (p.stat().st_size, p.stat().st_mtime_ns)
                    for p in root.rglob('*') if p.is_file() and
                    not any(x in p.parts for x in ['logs', 'audit', '__pycache__'])}
        before = snapshot() if changed_only else {}
        record = dict(status='running', command=command, outputs={})
        write_json(self.audit / (name + '.json'), record)
        log = self.audit / (name + '.log')
        print(name + ': running; log=' + str(log), flush=True)
        start = time.monotonic()
        with log.open('a', encoding='utf-8') as stream:
            proc = subprocess.Popen(command, cwd=self.cwd, env=self.env, stdout=stream,
                                    stderr=subprocess.STDOUT, start_new_session=os.name == 'posix')
            try:
                while True:
                    try:
                        code = proc.wait(timeout=30)
                        break
                    except subprocess.TimeoutExpired:
                        print(f'  {name}: {int(time.monotonic()-start)} s; log={log}', flush=True)
            except BaseException:
                if os.name == 'posix':
                    os.killpg(proc.pid, signal.SIGTERM)
                else:
                    proc.terminate()
                try:
                    proc.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    if os.name == 'posix': os.killpg(proc.pid, signal.SIGKILL)
                    else: proc.kill()
                    proc.wait()
                record['status'] = 'interrupted'
                write_json(self.audit / (name + '.json'), record)
                raise
        record.update(status='complete' if code == 0 else 'failed', exit_code=code,
                      seconds=round(time.monotonic()-start, 2))
        if code == 0:
            after = snapshot()
            record['outputs'] = {name: sha(root / name) for name, stat in after.items()
                                 if not changed_only or before.get(name) != stat}
        write_json(self.audit / (name + '.json'), record)
        if code:
            print(log.read_text(encoding='utf-8', errors='replace')[-5000:], flush=True)
            raise RuntimeError('Step failed: ' + name + '; resume this run after checking the log.')
        print(name + ': complete', flush=True)


def code_identity(root):
    """Record the current source tree; no packaged checksum is required."""
    root=Path(root)
    paths=[root/"run.py",root/"run.sh",root/"requirements.txt"]
    for folder in ["scripts","workflows","config"]:
        paths.extend(p for p in (root/folder).rglob("*") if p.is_file() and "__pycache__" not in p.parts)
    return {p.relative_to(root).as_posix():sha(p) for p in sorted(paths) if p.is_file()}


def require_analysis_python():
    import platform,sys
    if platform.system()!="Linux" or sys.version_info[:2]!=(3,12):
        raise RuntimeError("Run the analyses on Linux with Python 3.12 (paper2-r442).")


def resolve_rscript(value="Rscript"):
    import shutil
    command=shutil.which(str(value))
    if command is None:raise FileNotFoundError("Rscript not found; activate paper2-r442.")
    return command


def environment_versions(rscript="Rscript"):
    from importlib.metadata import version
    import sys
    names=["numpy","pandas","scipy","statsmodels","patsy","matplotlib","h5py","openpyxl","PyYAML"]
    packages={name:version(name) for name in names}
    if packages['statsmodels']!='0.14.6':
        raise RuntimeError('The manuscript models require statsmodels 0.14.6; activate paper2-r442.')
    r=subprocess.check_output([resolve_rscript(rscript),"--vanilla","-e",
        'if(getRversion()!="4.4.2")stop("Use R 4.4.2"); '
        'packages<-c("yaml","readxl","jsonlite","statmod","dtangle","BiocVersion","limma","edgeR",'
        '"Biobase","AnnotationDbi","org.Hs.eg.db","huex10sttranscriptcluster.db","oligo","pd.huex.1.0.st.v2","preprocessCore"); '
        'missing<-packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]; '
        'if(length(missing))stop("Missing R packages: ",paste(missing,collapse=", ")); '
        'if(!startsWith(as.character(packageVersion("BiocVersion")),"3.20."))stop("Use Bioconductor 3.20"); '
        'if(as.character(packageVersion("oligo"))!="1.70.1")stop("Use oligo 1.70.1"); '
        'cat(R.version.string); for(p in packages)cat("\\n",p,as.character(packageVersion(p)))'],text=True)
    return dict(python=sys.version,packages=packages,R=r,threads=1)
