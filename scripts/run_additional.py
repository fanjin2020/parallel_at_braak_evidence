"""Run the additional analyses alongside the original repository."""

import argparse
from pathlib import Path
import subprocess
import sys

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--output-dir',type=Path,default=Path('results/extensions'))
    p.add_argument('--bootstrap',type=int,default=2000,help='Participant bootstrap replicates (minimum 100).')
    p.add_argument('--selection-bootstrap',type=int,default=1000,help='Full-process weighting bootstrap replicates (minimum 200).')
    p.add_argument('--seed',type=int,default=20260915)
    a=p.parse_args()
    if a.bootstrap<100 or a.selection_bootstrap<200:
        p.error('Use at least 100 participant and 200 selection bootstrap replicates.')
    root=a.project_root.resolve()
    out=a.output_dir if a.output_dir.is_absolute() else root/a.output_dir
    out=out.resolve()
    if out==root or out==root/'data' or root/'data' in out.parents:
        p.error('Output must be outside the input data directory.')
    if out.exists() and any(out.iterdir()):p.error('Choose an empty output directory.')
    required=['adni_pet_aligned_master.csv','raw/CDR.csv','raw/DXSUM.csv','raw/ADSL.csv','raw/DATA_DOWNLOADED_DATE.csv']
    for name in required:
        if not (root/'data/private/adni'/name).is_file():p.error('Missing input: '+name)
    for name in ['numpy','pandas','scipy','statsmodels','patsy']:__import__(name)
    here=Path(__file__).resolve().parent
    jobs=[('cross_sectional',True,[]),('tau_sensitivity',True,[]),('spline_stability',True,[]),
          ('longitudinal',False,[]),('clinical_tables',False,[]),
          ('longitudinal_sensitivity',False,['--bootstrap',str(a.bootstrap),'--seed',str(a.seed)]),
          ('selection_sensitivity',False,['--bootstrap',str(a.selection_bootstrap),'--seed',str(a.seed)]),
          ('stage_sensitivity',False,[])]
    for name,master,extra in jobs:
        source=['--master',str(root/'data/private/adni/adni_pet_aligned_master.csv')] if master else ['--project-root',str(root)]
        print('\nRunning '+name,flush=True)
        subprocess.run([sys.executable,str(here/(name+'.py')),*source,'--output-dir',str(out/name),*extra],check=True)
    print('\nResults: '+str(out))

if __name__=='__main__':main()
