"""Publication summaries calculated directly from fitted result tables."""
from pathlib import Path
import pandas as pd
from export_io import read_csv_exact
from pipeline_runtime import inventory, code_identity, write_json


PROCESSING_MODELS = [
    'original_provided_provided', 'technical_qc_pass_provided',
    'technical_qc_pass_reprocessed', 'exclude_rle_flagged_provided',
    'exclude_rle_flagged_reprocessed', 'exclude_technical_failure_donors_provided',
    'exclude_technical_failure_donors_reprocessed',
]


def processing_counts(base):
    qc=Path(base)/'cel_qc'
    cohorts=pd.read_csv(qc/'model_cohorts.csv').set_index('model')
    rows=[]
    for model in PROCESSING_MODELS:
        path=qc/('features_'+model+'.csv')
        values=pd.read_csv(path,float_precision='round_trip')
        cohort=cohorts.loc[model]
        rows.append(dict(model=model,n_samples=int(cohort.n_samples),n_donors=int(cohort.n_donors),
            n_features=len(values),n_features_q_lt_0_05=int(values.q.lt(.05).sum()),
            q_threshold=.05,multiple_testing='BH over all tested features',source_result='cel_qc/'+path.name))
    return pd.DataFrame(rows)


def display_q(value):
    value=float(value)
    return format(value,'.2e') if 0<value<.0001 else format(value,'.4f')


def table4_display(target):
    rows=[]
    for gene,group in target.groupby('gene_symbol',sort=True):
        g=group.set_index('model');b=g.loc['baseline']
        rows.append(dict(gene_symbol=gene,entrez_id=str(b.entrez_id),
            baseline_log2FC=format(float(b.log2FC),'.3f'),
            baseline_95CI=f'{float(b.ci_low):.3f} to {float(b.ci_high):.3f}',
            baseline_q=display_q(b.q),joint_CA_q=display_q(g.loc['joint_CA_joint_quantile','q']),
            joint_NG_q=display_q(g.loc['joint_NG_joint_quantile','q'])))
    return pd.DataFrame(rows)


def reporting_artifacts(base,tables):
    from prepare_reporting_tables import source_definitions
    tables=Path(tables);tables.mkdir(parents=True,exist_ok=True)
    processing_counts(base).to_csv(tables/'processing_feature_counts.csv',index=False)
    target=pd.read_csv(tables/'QC_nonY_targets.csv',dtype={'entrez_id':str})
    table4_display(target).to_csv(tables/'Table4_display.csv',index=False)
    source_definitions().to_csv(tables/'data_sources.csv',index=False)


def final_metadata(base,survival,out,code_root,source_records):
    base,survival,out=map(Path,[base,survival,out])
    dictionary=read_csv_exact(out/'Variable_Dictionary.csv')
    targets=read_csv_exact(survival/'GSE174367_target_results.csv')
    mg=targets.loc[targets.cell_type.eq('MG')&targets.model.eq('age_sex_batch_RIN')]
    stage=read_csv_exact(base/'enhancement/clinical/stage_propagated_intervals.csv')
    cv=read_csv_exact(base/'enhancement/clinical/landmark_0_refitted_CV.csv')
    import hashlib,json
    identity=hashlib.sha256(json.dumps(code_identity(code_root),sort_keys=True).encode()).hexdigest()
    payload=dict(code_version='V5',code_identity_sha256=identity,source_records=source_records,
        scope='Aggregate results, figure sources and variable definitions.',
        dictionary=dict(field_instances=len(dictionary),scopes=int(dictionary.scope.nunique()),
            roles=dictionary.role.value_counts().to_dict(),definition_status=dictionary.definition_status.value_counts().to_dict()),
        threshold_bootstrap=stage[['outcome','contrast','requested','valid']].to_dict('records'),
        pet_day_refitted_bootstrap=cv[['family','requested','valid']].to_dict('records'),
        microglial_untested=mg.loc[mg.status.ne('tested'),['gene_symbol','status']].to_dict('records'))
    write_json(out/'PUBLICATION_SUMMARY.json',payload)
    pd.DataFrame([dict(file=p,sha256=h) for p,h in inventory(out).items()
                  if p!='SUBMISSION_SHA256SUMS.csv']).to_csv(out/'SUBMISSION_SHA256SUMS.csv',index=False)
    return payload
