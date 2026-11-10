from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd


class NotEstimable(Exception):
    pass


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False), encoding='utf-8')


def write_csv(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False, encoding='utf-8-sig')


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def finite_matrix(x):
    x = np.asarray(x, float)
    if not np.isfinite(x).all() or np.linalg.matrix_rank(x) != x.shape[1]:
        raise NotEstimable('Nonfinite or rank-deficient design; no terms were removed')
    return x
