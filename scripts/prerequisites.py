import importlib,importlib.metadata as metadata,shutil,subprocess
from pathlib import Path

def check(rscript='Rscript'):
    if not shutil.which(rscript):raise RuntimeError('Activate paper2-r442: Rscript not found')
    packages=['numpy','pandas','scipy','statsmodels','patsy','matplotlib','h5py','openpyxl','PyYAML']
    for n in packages:importlib.import_module('yaml' if n=='PyYAML' else n)
    versions={n:metadata.version(n) for n in packages}
    if versions['statsmodels']!='0.14.6':raise RuntimeError('The manuscript models use statsmodels 0.14.6')
    subprocess.run([rscript,'--vanilla','-e','stopifnot(getRversion()=="4.4.2"); for(p in c("limma","edgeR","statmod")) if(!requireNamespace(p,quietly=TRUE)) stop("Missing ",p)'],check=True)
    return versions
