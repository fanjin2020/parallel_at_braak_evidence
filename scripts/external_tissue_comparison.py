"""Compare complete external GO families with every planned QC model."""
from pathlib import Path
import numpy as np
import pandas as pd
from analysis_utils import csv

MODELS=['baseline','scan_date','joint_CA_joint_quantile','joint_NG_joint_quantile']

def compare(discovery,external,cohort,external_model):
    if discovery.go_id.duplicated().any() or external.go_id.duplicated().any():raise ValueError('Duplicated GO identifiers')
    for d in [discovery,external]:
        if not np.isfinite(d[['p','q']]).all().all() or not d.direction.isin(['Up','Down']).all():
            raise ValueError('Invalid full-family GO statistics')
    z=discovery[['go_id','direction','p','q']].merge(external[['go_id','direction','p','q']],on='go_id',how='left',suffixes=('_QC','_external'),validate='one_to_one')
    z['external_testable']=z.p_external.notna();z['QC_significant']=z.q_QC.lt(.05)
    z['same_direction']=z.external_testable & z.direction_QC.eq(z.direction_external)
    z['external_significant']=z.external_testable & z.q_external.lt(.05)
    z['cohort']=cohort;z['external_model']=external_model
    selected=z.QC_significant;covered=selected & z.external_testable
    summary=dict(cohort=cohort,external_model=external_model,QC_tested=len(discovery),external_tested=len(external),
        QC_significant=int(selected.sum()),QC_significant_testable=int(covered.sum()),QC_significant_unmeasured=int((selected&~z.external_testable).sum()),
        same_direction_significant=int((covered&z.same_direction&z.external_significant).sum()),
        opposite_direction_significant=int((covered&~z.same_direction&z.external_significant).sum()),
        testable_not_external_significant=int((covered&~z.external_significant).sum()),
        interpretation='Related AD/control comparison, not matched Braak replication' if cohort=='GSE44770' else 'Braak V-VI versus 0 in related inferior temporal cortex; not all discovery regions')
    return z,summary

def run(root,out):
    root=Path(root);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(root/'targeted/qc_joint/GO_all_models.csv')
    e1=pd.read_csv(root/'revision/external/external_pathway_results.csv')
    e2=pd.read_csv(root/'enhancement/external/external_GO_all_models.csv')
    tables=[];summaries=[]
    for model in MODELS:
        qc=d.loc[d.model.eq(model)]
        if qc.empty:raise ValueError('Missing QC model '+model)
        externals=[('GSE44770','age_sex_PMI_pH_RIN_batch',e1)]+[('GSE106241',name,z) for name,z in e2.groupby('model')]
        for cohort,label,e in externals:
            table,summary=compare(qc,e,cohort,label);table['QC_model']=model;summary['QC_model']=model
            tables.append(table);summaries.append(summary)
    csv(out/'GO_external_QC_full.csv',pd.concat(tables,ignore_index=True));csv(out/'GO_external_QC_summary.csv',summaries)
    csv(out/'analysis_status.csv',[dict(component='external_QC_GO',status='completed',reason='All full-family external q values retained; no subset BH recalculation')])
    (out/'methods.txt').write_text('Discovery selection uses each final QC model, not the historical 181-term set. External BH remains over each complete eligible external family. Unmeasured sets are retained with missing statistics, not interpreted as null. Related cohorts and RNA-quality models remain separate.\n',encoding='utf-8')
