"""Collect unrounded manuscript result sources without refitting or changing estimates."""
from pathlib import Path
import argparse
import csv
import shutil
import zipfile
import math
import statistics

SOURCES={
 'Table1A':['extensions/clinical_tables/table1a_clinical_by_at.csv'],
 'Table1B':['extensions/clinical_tables/table1b_longitudinal_cohorts.csv'],
 'Primary_cross_sectional':['manuscript_tables/table2_adni_confirmatory_associations.csv'],
 'FDG_APOE':['manuscript_tables/table3_adni_exploratory_results.csv'],
 'Robust_genes':['manuscript_tables/table4_braak_robust_genes.csv'],
 'GO_terms':['manuscript_tables/table5_stable_go_pathways.csv'],
 'Additional_cross_sectional':['extensions/cross_sectional/extension_results.csv'],
 'Longitudinal':['extensions/longitudinal/longitudinal_results.csv'],
 'S1':['descriptive_tables/S1_donor_characteristics.csv','descriptive_tables/donor_age_summary.csv'],
 'S2':['gse131617/Supplementary_Table_S2_feature_level_result.csv'],
 'S3':['descriptive_tables/S3_ATminus_Tplus_description.csv'],
 'S4':['adni/adni_results.csv'], 'S5':['adni/adni_results.csv'],
 'S6':['adni/adni_results.csv'],
 'S7':['gse_sensitivity/uty_sensitivity.csv','gse131617/braak_robust_genes.csv'],
 'S8':['pathways/stable_pathways.csv','manuscript_tables/table5_stable_go_pathways.csv'],
 'S9':['extensions/clinical_tables/table_s9_followup_comparison.csv'],
 'S10':['extensions/selection_sensitivity/selection_sensitivity_results.csv','extensions/selection_sensitivity/bootstrap_summary.csv'],
 'S11':['extensions/selection_sensitivity/weight_diagnostics.csv','extensions/selection_sensitivity/balance_against_eligible.csv','extensions/selection_sensitivity/selection_probability_summary.csv'],
 'S12':['extensions/stage_sensitivity/standardized_stage_change.csv','extensions/stage_sensitivity/stage_bounded_models.csv'],
 'S13':['extensions/stage_sensitivity/stage_time_calibration.csv'],
 'S14':['extensions/stage_sensitivity/tau_support_by_diagnosis.csv'],
 'Spline':['extensions/tau_sensitivity/nonlinearity_results.csv','extensions/tau_sensitivity/influence_sensitivity.csv','extensions/spline_stability/spline_stability_summary.csv'],
 'Longitudinal_sensitivity':['extensions/longitudinal_sensitivity/cluster_resampling_summary.csv','extensions/longitudinal_sensitivity/calendar_opportunity_sensitivity.csv'],
 'Expression_sources':['gse131617/braak_feature_results.csv','gse131617/braak_region_sensitivity.csv','gse131617/feature_annotation.csv'],
 'Pathway_sources':['pathways/ranked_pathways.csv'],
 'Sample_sizes':['adni/adni_model_counts.csv','extensions/longitudinal/sample_flow.csv','extensions/longitudinal/followup_coverage.csv'],
 'Tau_curves':['extensions/tau_sensitivity/adjusted_tau_curves.csv','extensions/spline_stability/spline_common_range_comparison.csv'],
 'Clinical_description':['extensions/clinical_tables/clinical_baseline_numeric.csv','extensions/clinical_tables/diagnosis_by_analysis_cohort.csv','extensions/clinical_tables/missingness_by_cohort.csv','extensions/clinical_tables/attrition_comparison_extended_numeric.csv','extensions/cross_sectional/age_coverage.csv'],
 'Longitudinal_model_results':['extensions/longitudinal/model_coefficients.csv','extensions/longitudinal/diagnosis_at_tau.csv','extensions/longitudinal/attrition_comparison.csv','extensions/longitudinal_sensitivity/model_diagnostics.csv','extensions/longitudinal_sensitivity/quadratic_time_tests.csv','extensions/longitudinal_sensitivity/standardized_change_contrasts.csv','extensions/longitudinal_sensitivity/sensitivity_model_coefficients.csv','extensions/longitudinal_sensitivity/time_bin_calibration.csv','extensions/longitudinal_sensitivity/calendar_opportunity_summary.csv'],
 'Stage_comparison':['extensions/stage_sensitivity/stage_model_comparison.csv','extensions/stage_sensitivity/stage_model_coefficients.csv','extensions/stage_sensitivity/sample_flow.csv'],
 'Selection_model':['extensions/selection_sensitivity/selection_model_coefficients.csv','extensions/selection_sensitivity/diagnosis_counts.csv'],
}

def read(path):
    with path.open(encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f))

def write(path,rows,fields):
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)

S2_PRIMARY=['feature_id','log2_fold_change','average_expression','moderated_t','p_value','fdr_bh','B']
S2_REGIONAL=['direction_robust','log2_fold_change_leave_out_EC',
             'log2_fold_change_leave_out_FC','log2_fold_change_leave_out_TC']
S2_FIELDS=S2_PRIMARY+S2_REGIONAL+['entrez_id','mapping_status']

def build_s2(primary_path, regional_path, annotation_path, output_path):
    """Join this run's public sources by feature ID; preserve numerical strings."""
    def indexed(path, required):
        if not path.is_file():
            raise FileNotFoundError('Required S2 source missing from this run: '+str(path))
        with path.open(encoding='utf-8-sig', newline='') as handle:
            reader=csv.DictReader(handle)
            missing=set(required)-set(reader.fieldnames or [])
            if missing: raise ValueError(f'{path}: missing columns {sorted(missing)}')
            rows=list(reader)
        result={}
        for row in rows:
            key=row['feature_id']
            if not key or key in result:
                raise ValueError(f'{path}: blank or duplicate feature_id {key!r}')
            result[key]=row
        if not result: raise ValueError('Empty S2 source: '+str(path))
        return result
    primary=indexed(primary_path,S2_PRIMARY)
    regional=indexed(regional_path,S2_PRIMARY+S2_REGIONAL)
    annotation=indexed(annotation_path,['feature_id','entrez_id','mapping_status'])
    if set(primary)!=set(regional):
        raise ValueError('S2 primary and regional feature IDs differ.')
    if set(primary)-set(annotation):
        raise ValueError('S2 annotation does not cover every primary feature ID.')
    rows=[]
    for key, source in primary.items():
        region=regional[key]; annot=annotation[key]
        for field in S2_PRIMARY[1:]:
            a,b=source[field],region[field]
            if a==b: continue
            try:
                x,y=float(a),float(b)
                equal=(math.isnan(x) and math.isnan(y)) or math.isclose(x,y,rel_tol=1e-12,abs_tol=0.0)
            except (ValueError,TypeError): equal=False
            if not equal: raise ValueError(f'S2 primary/regional mismatch: {key}, {field}')
        direction=region['direction_robust'].strip().lower()
        if direction not in {'true','false'}:
            raise ValueError(f'Invalid direction_robust for {key}: {direction!r}')
        status=annot['mapping_status']; entrez=annot['entrez_id']
        ids=entrez.split(';') if entrez else []
        valid=(status=='unmapped_entrez' and not ids or
               status=='unique_entrez_id' and len(ids)==1 or
               status=='ambiguous_multiple_entrez' and len(ids)>1)
        if not valid or any(not value.isdigit() for value in ids) or len(set(ids))!=len(ids):
            raise ValueError(f'Inconsistent Entrez annotation for {key}')
        row={field:source[field] for field in S2_PRIMARY}
        row.update({field:region[field] for field in S2_REGIONAL})
        row['direction_robust']='True' if direction=='true' else 'False'
        row.update(entrez_id=entrez,mapping_status=status)
        rows.append(row)
    write(output_path,rows,S2_FIELDS)
    return rows

def collect(run,out):
    base=run/'results'
    if out==base or base in out.parents:
        raise ValueError('Collection directory must be outside the analysis results directory.')
    if out.exists() and any(out.iterdir()): raise ValueError('Collection directory must be empty.')
    if out.with_suffix('.zip').exists(): raise FileExistsError(out.with_suffix('.zip'))
    annotation=run/'data/public/gse131617/feature_annotation.csv'
    build_s2(base/'gse131617/braak_feature_results.csv',
             base/'gse131617/braak_region_sensitivity.csv',annotation,
             base/SOURCES['S2'][0])

    shutil.copy2(annotation,base/'gse131617/feature_annotation.csv')

    master=read(run/'data/private/adni/adni_pet_aligned_master.csv')
    manifest=read(run/'data/public/gse131617/sample_manifest.csv')
    def values(rows,key):
        ans=[]
        for r in rows:
            try: v=float(r[key])
            except (ValueError,TypeError): continue
            if math.isfinite(v): ans.append(v)
        return ans
    desc=[]
    subset=[r for r in master if r['at_stage']=='A-T+']
    for key in ['age_at_tau','sex_male','education_years','apoe4_carrier','hippocampus_icv','cdrsb','hippocampus_icv_days_from_tau','cdrsb_days_from_tau']:
        v=values(subset,key)
        desc.append(dict(variable=key,n=len(v),mean=statistics.mean(v) if v else '',sd=statistics.stdev(v) if len(v)>1 else '',median=statistics.median(v) if v else '',minimum=min(v) if v else '',maximum=max(v) if v else ''))
    donors={}
    donor_keys=['braak_stage','sex_male','age_at_death','pmi_hours','apoe4_carrier']
    for r in manifest:
        key=r['donor_key']; row={k:r[k] for k in donor_keys}
        if key in donors and donors[key]!=row: raise ValueError('Inconsistent donor covariates in manifest')
        donors[key]=row
    donor_rows=[]
    donor_age_rows=[]
    for stage in ['0','I-II','III-IV','V-VI']:
        group=[r for r in donors.values() if r['braak_stage']==stage]
        age=values(group,'age_at_death');pmi=values(group,'pmi_hours');sex=values(group,'sex_male');apoe=values(group,'apoe4_carrier')
        donor_rows.append(dict(braak_stage=stage,donors=len(group),male=sex.count(1),female=sex.count(0),apoe4_carriers=apoe.count(1),age_median=statistics.median(age) if age else '',pmi_median=statistics.median(pmi) if pmi else ''))
        donor_age_rows.append(dict(braak_stage=stage,donors=len(group),n_age=len(age),
                                   age_mean=statistics.mean(age) if age else '',
                                   age_sd=statistics.stdev(age) if len(age)>1 else '',
                                   age_median=statistics.median(age) if age else ''))
    generated=base/'descriptive_tables';generated.mkdir(exist_ok=True)
    write(generated/'S1_donor_characteristics.csv',donor_rows,list(donor_rows[0]))
    write(generated/'donor_age_summary.csv',donor_age_rows,list(donor_age_rows[0]))
    write(generated/'S3_ATminus_Tplus_description.csv',desc,list(desc[0]))
    missing=sorted({n for names in SOURCES.values() for n in names if not (base/n).is_file()})
    if missing: raise FileNotFoundError('Required result files missing:\n'+'\n'.join(missing))
    if out.exists() and any(out.iterdir()): raise ValueError('Collection directory must be empty.')
    files=sorted({base/name for names in SOURCES.values() for name in names})

    forbidden={'rid','ptid','subjid','subject_id','participant_id','patient_id'}
    for f in files:
        with f.open(encoding='utf-8-sig',newline='') as h: fields=next(csv.reader(h),[])
        if forbidden.intersection(x.lower() for x in fields):
            raise ValueError('Participant-level identifier column found; not bundled: '+str(f))
    out.mkdir(parents=True,exist_ok=True)
    for f in files:
        rel=f.relative_to(base); dest=out/'source_tables'/rel
        dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(f,dest)
    mapping=[{'paper_item':label,'source_file':'source_tables/'+name} for label,names in SOURCES.items() for name in names]
    write(out/'paper_table_sources.csv',mapping,['paper_item','source_file'])
    (out/'README.txt').write_text(
      'Numerical sources for the manuscript tables and analyses.\n'
      'paper_table_sources.csv maps manuscript items to source files.\n'
      'Source filenames retain their analysis-specific table numbers.\n'
      'S4-S6 sources contain several analyses; select the matching outcome/analysis rows.\n'
      'These files are result sources, not typeset Word tables or final figure PDFs.\n'
      'S2 has 13 columns joined from this run: primary results, regional sensitivity and feature annotation.\n'
      'Numeric source strings are preserved. Public feature_annotation.csv is included.\n'
      'The original 7-column feature results are also included.\n'
      'Input participant files and local configuration are excluded. Review before public sharing.\n',encoding='utf-8')
    archive=out.with_suffix('.zip')
    if archive.exists(): raise FileExistsError(archive)
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
        for f in sorted(out.rglob('*')):
            if f.is_file(): z.write(f,Path(run.name+'_'+out.name)/f.relative_to(out))
    print(f'Collected {len(files)} result files: {archive}')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-dir',type=Path,required=True)
    p.add_argument('--output-dir',type=Path)
    a=p.parse_args();run=a.run_dir.resolve()
    collect(run,a.output_dir.resolve() if a.output_dir else run/'paper_results')

if __name__=='__main__': main()
