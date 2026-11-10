"""Derive compact publication evidence from complete fitted results."""
from pathlib import Path
import csv
import json
import numpy as np
import pandas as pd
from pipeline_runtime import inventory, sha, write_json
from source_schema import read_csv_exact

PROCESSING_MODELS = [
    'original_provided_provided', 'technical_qc_pass_provided',
    'technical_qc_pass_reprocessed', 'exclude_rle_flagged_provided',
    'exclude_rle_flagged_reprocessed', 'exclude_technical_failure_donors_provided',
    'exclude_technical_failure_donors_reprocessed',
]


def bh(p):
    p = np.asarray(p, dtype=float)
    if not len(p) or not np.isfinite(p).all() or np.any((p < 0) | (p > 1)):
        raise ValueError('Invalid full-background P values')
    order = np.argsort(p, kind='stable')
    ranked = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p)+1))[::-1])[::-1]
    out = np.empty(len(p)); out[order] = np.minimum(1, ranked)
    return out


def processing_counts(base):
    """Count all features, never just the selected targets or GO sets."""
    base = Path(base); qc = base/'cel_qc'
    cohorts = pd.read_csv(qc/'model_cohorts.csv')
    states = pd.read_csv(qc/'analysis_status.csv')
    if cohorts.model.duplicated().any() or set(cohorts.model) != set(PROCESSING_MODELS):
        raise ValueError('Table S20A requires all seven processing cohorts')
    status = states.set_index('component').status
    if status.index.duplicated().any() or not status.reindex(PROCESSING_MODELS).eq('completed').all():
        raise ValueError('Table S20A processing models are incomplete')
    target = pd.read_csv(qc/'target_processing_comparison.csv', dtype={'feature_id':str})
    baseline = pd.read_csv(base/'targeted/qc_joint/features_baseline.csv', dtype={'id':str})
    family = set(baseline.id)
    if len(baseline) != 22011 or len(family) != len(baseline):
        raise ValueError('Invalid QC baseline full-feature family')
    rows = []
    for model in PROCESSING_MODELS:
        path = qc/('features_'+model+'.csv')
        if not path.is_file():
            raise FileNotFoundError('Missing full-feature evidence for S20A: '+str(path)+
                                    '. Retain this model output; target-only tables cannot replace it.')
        d = pd.read_csv(path, dtype={'feature_id':str}, float_precision='round_trip')
        if d.feature_id.duplicated().any() or set(d.feature_id) != family:
            raise ValueError('Processing feature background differs: '+model)
        values = d[['log2FC','ci_low','ci_high','p','q']].to_numpy(float)
        if not np.isfinite(values).all() or not ((d.ci_low <= d.log2FC) & (d.log2FC <= d.ci_high)).all():
            raise ValueError('Invalid full-feature estimates: '+model)
        if not np.allclose(bh(d.p), d.q, rtol=1e-8, atol=1e-12):
            raise ValueError('Full-feature BH mismatch: '+model)
        t = target.loc[target.model.eq(model)]
        if len(t) != 17 or t.feature_id.duplicated().any():
            raise ValueError('Incomplete processing target set: '+model)
        matched = d.set_index('feature_id').loc[t.feature_id]
        if not np.allclose(matched[['log2FC','ci_low','ci_high','p','q']],
                           t[['log2FC','ci_low','ci_high','p','q']], rtol=1e-8, atol=1e-12):
            raise ValueError('Target/full-background mismatch: '+model)
        cohort = cohorts.set_index('model').loc[model]
        if not t.n_samples.eq(cohort.n_samples).all() or not t.n_donors.eq(cohort.n_donors).all():
            raise ValueError('Processing cohort count mismatch: '+model)
        if model == 'technical_qc_pass_reprocessed':
            z = baseline.set_index('id').loc[d.feature_id]
            if not np.allclose(z[['log2FC','ci_low','ci_high','p','q']], values,
                               rtol=1e-7, atol=1e-12):
                raise ValueError('Data 2 and processing QC baseline differ')
        rows.append(dict(model=model,n_samples=int(cohort.n_samples),n_donors=int(cohort.n_donors),
                         n_features=len(d),n_features_q_lt_0_05=int(d.q.lt(.05).sum()),
                         q_threshold=.05,multiple_testing='BH over all tested features',
                         source_result='cel_qc/'+path.name))
    return pd.DataFrame(rows)


def display_q(value):
    value = float(value)
    if not np.isfinite(value) or not 0 <= value <= 1:
        raise ValueError('Invalid display q')
    return format(value, '.2e') if 0 < value < .0001 else format(value, '.4f')


def table4_display(target):
    models = ['baseline','scan_date','joint_CA_joint_quantile','joint_NG_joint_quantile']
    if len(target) != 64 or target.duplicated(['gene_symbol','model']).any():
        raise ValueError('Table 4 requires 16 targets and four unique models')
    rows = []
    for gene, group in target.groupby('gene_symbol', sort=True):
        if gene == 'UTY' or set(group.model) != set(models):
            raise ValueError('Unexpected Table 4 target/model set')
        g = group.set_index('model'); b = g.loc['baseline']
        rows.append(dict(gene_symbol=gene,entrez_id=str(b.entrez_id),
                         baseline_log2FC=format(float(b.log2FC),'.3f'),
                         baseline_95CI=f'{float(b.ci_low):.3f} to {float(b.ci_high):.3f}',
                         baseline_q=display_q(b.q),
                         joint_CA_q=display_q(g.loc['joint_CA_joint_quantile','q']),
                         joint_NG_q=display_q(g.loc['joint_NG_joint_quantile','q'])))
    return pd.DataFrame(rows)


def reporting_artifacts(base, tables):
    """Refresh reporting only; do not refit models or change the input workspace."""
    from prepare_reporting_tables import source_definitions
    tables = Path(tables); tables.mkdir(parents=True, exist_ok=True)
    processing_counts(base).to_csv(tables/'processing_feature_counts.csv', index=False)
    target = pd.read_csv(tables/'QC_nonY_targets.csv', dtype={'entrez_id':str})
    table4_display(target).to_csv(tables/'Table4_display.csv', index=False)
    source_definitions().to_csv(tables/'data_sources.csv', index=False)


def validate_source_plan(root, ids):
    """Keep every source citation required by the current manuscript contract."""
    path = Path(root)/'config/manuscript_source_requirements.csv'
    with path.open(encoding='utf-8-sig',newline='') as f:
        requirements = list(csv.DictReader(f))
    missing = sorted({r['record_id'] for r in requirements} - set(ids))
    if missing:
        raise ValueError('Manuscript-cited sources missing from export: '+', '.join(missing))
    return len(requirements)


def validate_stage_figure(base, presentation):
    left = pd.read_csv(Path(base)/'enhancement/clinical/stage_propagated_intervals.csv')
    right = pd.read_csv(Path(presentation)/'figure_sources/FigureS19_source.csv')
    keys = ['outcome','contrast']; cols = ['requested','valid','estimate','ci_low','ci_high']
    if left.duplicated(keys).any() or right.duplicated(keys).any():
        raise ValueError('Duplicate Figure S8 source keys')
    a = left.set_index(keys).sort_index(); b = right.set_index(keys).sort_index()
    counts = a[['requested','valid']].to_numpy(float)
    if (not np.isfinite(counts).all() or not np.equal(counts,np.floor(counts)).all()
            or (a.requested<=0).any() or (a.valid<0).any() or (a.valid>a.requested).any()):
        raise ValueError('Invalid threshold-propagation bootstrap counts')
    if not a.index.equals(b.index) or not np.allclose(a[cols], b[cols], rtol=1e-12, atol=1e-14):
        raise ValueError('Figure S8 source differs from Table S5D')
    if not a.status.eq('ok').all() or (a.valid/a.requested < .9).any():
        raise ValueError('Incomplete threshold-propagation bootstrap')


def final_metadata(base, survival, out, code_root, source_records, checks):
    """Public run facts and hashes; no private paths or participant records."""
    base, survival, out = map(Path,[base,survival,out])
    dictionary = read_csv_exact(out/'Variable_Dictionary.csv')
    targets = read_csv_exact(survival/'GSE174367_target_results.csv')
    if targets.status.eq('unmeasured_filtered_or_ambiguous').any():
        raise ValueError('Legacy combined single-nucleus status needs verified reclassification; do not guess its cause')
    mg = targets.loc[targets.cell_type.eq('MG') & targets.model.eq('age_sex_batch_RIN')]
    stage = read_csv_exact(base/'enhancement/clinical/stage_propagated_intervals.csv')
    cv = read_csv_exact(base/'enhancement/clinical/landmark_0_refitted_CV.csv')
    payload = dict(
        code_manifest_sha256=sha(Path(code_root)/'CODE_SHA256SUMS.csv'),
        source_records=source_records,
        source_validation=checks,
        document_checked=False,
        scope='Generated data and figures only. Word files and the response letter are not edited or certified by this run.',
        dictionary=dict(field_instances=len(dictionary),scopes=int(dictionary.scope.nunique()),
                        roles=dictionary.role.value_counts().to_dict(),
                        definition_status=dictionary.definition_status.value_counts().to_dict()),
        threshold_bootstrap=stage[['outcome','contrast','requested','valid']].to_dict('records'),
        pet_day_refitted_bootstrap=cv[['family','requested','valid']].to_dict('records'),
        microglial_untested=mg.loc[mg.status.ne('tested'),['gene_symbol','status']].to_dict('records'),
    )
    write_json(out/'PUBLICATION_SUMMARY.json',payload)
    entries = [dict(file=p,sha256=h) for p,h in inventory(out).items() if p!='SUBMISSION_SHA256SUMS.csv']
    pd.DataFrame(entries).to_csv(out/'SUBMISSION_SHA256SUMS.csv',index=False)
    return payload
