"""Checkpointed, local-only helpers for the reviewer extension."""
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import pandas as pd

class NotEstimable(ValueError):
    pass

def dump(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp')
    def convert(x):
        if isinstance(x,np.generic):return convert(x.item())
        if isinstance(x,Path):return str(x)
        if isinstance(x,float) and not math.isfinite(x):return None
        if isinstance(x,dict):return {k:convert(v) for k,v in x.items()}
        if isinstance(x,(tuple,list)):return [convert(v) for v in x]
        return x
    tmp.write_text(json.dumps(convert(value),ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    tmp.replace(path)

def csv(path,rows,columns=None):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    df=rows if isinstance(rows,pd.DataFrame) else pd.DataFrame(rows,columns=columns)
    if not len(df.columns):df=pd.DataFrame(columns=['status','reason'])
    df.to_csv(path,index=False,lineterminator='\n')

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def checkpoint(out,key,fn):
    file=Path(out)/'checkpoints'/(key+'.json')
    if file.exists():return json.loads(file.read_text(encoding='utf-8'))
    result=fn();dump(file,result);return result

def safe_fit(fn):
    try:return dict(status='completed',value=fn(),reason='')
    except (NotEstimable,np.linalg.LinAlgError) as e:
        return dict(status='not_estimable',reason=str(e),value=None)

def bootstrap_people(df,rng):
    groups=dict(tuple(df.groupby('rid',sort=False)));ids=list(groups)
    blocks=[]
    for i,key in enumerate(rng.choice(ids,len(ids),replace=True)):
        d=groups[key].copy();d['original_rid']=key;d['rid']=str(i);blocks.append(d)
    return pd.concat(blocks,ignore_index=True)
