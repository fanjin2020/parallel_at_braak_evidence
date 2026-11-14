"""Model tasks; no figures, publication archives or schema-directory scans."""
import argparse
import json
from pathlib import Path
from analysis_utils import csv, dump


def observation(root, out):
    import pandas as pd
    from observation_models import run
    from observation_statistics import annotate_effects
    from observation_dictionary import observation_definition
    plan = json.loads((Path(__file__).resolve().parents[1] / 'config/observation_plan.json').read_text())
    run(root, out, plan)
    # Constant-variable SMD handling and bootstrap status annotations are part of the fit module.
    effects = pd.read_csv(out/'effects.csv', float_precision='round_trip')
    draws = pd.read_csv(out/'bootstrap.csv', float_precision='round_trip')
    checked, counts = annotate_effects(effects, draws)
    csv(out/'Observation_Effects.csv', checked)
    parts = []
    for name in ['coverage','observation_models','probability_summary','weight_summary','balance','bootstrap_summary']:
        d = pd.read_csv(out/(name+'.csv'))
        d.insert(0, 'section', name)
        parts.append(d)
    csv(out/'Observation_Diagnostics.csv', pd.concat(parts, ignore_index=True))
    dump(out/'analysis_plan.json', plan)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('task', choices=['external','clinical','mnar','prediction','observation','presentation'])
    p.add_argument('--primary', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--models', type=Path)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    if a.task == 'external':
        from external_tissue_comparison import run
        run(a.primary, a.output)
    elif a.task == 'clinical':
        from longitudinal_observation import run
        run(a.primary, a.output, bootstrap=500)
    elif a.task == 'mnar':
        from covariate_mnar import run
        run(a.primary, a.output, imputations=100, diagnostic=False)
    elif a.task == 'prediction':
        from prediction_calibration import run
        run(a.primary, a.output, bootstrap=200, repeats=5)
    elif a.task == 'observation':
        observation(a.primary, a.output)
    else:
        from manuscript_figures import prepare
        prepare(a.primary, a.models, a.output/'publication_tables')


if __name__ == '__main__': main()
