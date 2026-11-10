"""Validate a user-supplied ADNI archive without network access."""
import hashlib
import json
from pathlib import PurePosixPath
import tarfile

def check_bundle(path):
    required={'STUDYSUM.csv','REGISTRY.csv','ADVERSE.csv','RECADV.csv','TREATDIS.csv','DATADIC.csv','DATA_DOWNLOADED_DATE.csv'}
    allowed=required|{'OTDISC.csv','inventory.csv','R_sessionInfo.txt','SHA256.json'}
    with tarfile.open(path,'r:gz') as z:
        entries={}
        for m in z:
            if m.isdir():continue
            parts=PurePosixPath(m.name).parts
            if not m.isfile() or len(parts)!=2 or parts[0]!='ADNI_vital_local_20261101_PRIVATE' or parts[1] not in allowed or parts[1] in entries:
                raise ValueError('Unexpected private archive entry')
            entries[parts[1]]=m
        if not (required|{'SHA256.json'}).issubset(entries):raise ValueError('Incomplete private archive')
        checks=json.load(z.extractfile(entries['SHA256.json']))
        if set(checks)!=set(entries)-{'SHA256.json'}:raise ValueError('Incomplete checksum manifest')
        for name,digest in checks.items():
            with z.extractfile(entries[name]) as stream:
                if hashlib.file_digest(stream,'sha256').hexdigest()!=digest:raise ValueError('Private input checksum mismatch: '+name)
