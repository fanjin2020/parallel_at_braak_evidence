"""One numbered figure set and allowlisted aggregate data for the manuscript."""
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import shutil
import sys
import zipfile
import numpy as np
import pandas as pd
from source_schema import check_public, read_csv_exact, csv_bytes, PATH_PATTERN
from pipeline_runtime import sha, write_json

ROOT=Path(__file__).resolve().parents[1]
SURVIVAL=['MNAR_scenarios','MNAR_sample_flow','GSE174367_model_status','GSE174367_target_results',
          'GSE174367_GO_all_models','GSE174367_interaction_status','GSE174367_interaction_targets',
          'GSE174367_interaction_GO','ADNI_ignored_record_counts','Attrition_reason_counts']


def fingerprint(frame):
    """Numeric digest at 10 significant figures; text formatting is checked separately."""
    numeric=hashlib.sha256();text=hashlib.sha256()
    number=re.compile(r'^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$')
    cols=sorted(frame.columns)
    for i,row in enumerate(frame[cols].itertuples(index=False,name=None)):
        for c,v in zip(cols,row):
            value=str(v);key=f'{i}\t{c}\t'
            if number.fullmatch(value):
                numeric.update((key+format(float(value),'.10g')+'\n').encode())
            else:text.update((key+value+'\n').encode())
    return dict(rows=len(frame),columns=len(cols),numeric_sha256=numeric.hexdigest(),text_sha256=text.hexdigest())


def records(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as stream:return list(csv.DictReader(stream))


def family(source):
    if 'figure_sources/' in source:return 'Figure_sources'
    if '/mi_matched/' in source:return 'Multiple_imputation'
    if 'sex_interaction/' in source:return 'Sex_interaction'
    if 'external/' in source or 'GSE' in source:return 'External_tissue'
    if 'GO_' in source or 'pathway' in source:return 'Pathways'
    if 'qc_' in source.lower() or '/tissue/' in source or source.startswith('publication/'):return 'Tissue_QC'
    if 'mnar/' in source or 'missingness/' in source:return 'Missing_data'
    if any(x in source for x in ['longitudinal/','followup','visit','selection_']):return 'Longitudinal'
    if 'prediction' in source or 'landmark_' in source:return 'Prediction'
    if 'definition' in source:return 'Definitions'
    return 'Clinical'


def group_tables(items):
    """Group within an analysis family, with at most 60 data columns per CSV."""
    groups=[]
    for item in items:
        columns=set(item['data'].columns);kind=family(item['source'])
        choices=[g for g in groups if g['family']==kind and len(columns|g['columns'])<=60]
        if choices:
            group=max(choices,key=lambda g:len(columns&g['columns']))
            group['columns']|=columns;group['items'].append(item)
        else:groups.append(dict(family=kind,columns=columns,items=[item]))
    return groups


def export(base,models,presentation,observation,survival,out,audit,draw_figures=True):
    base,models,presentation,observation,survival,out,audit=map(Path,[base,models,presentation,observation,survival,out,audit])
    out.mkdir(parents=True,exist_ok=True);audit.mkdir(parents=True,exist_ok=True)
    from publication_validation import reporting_artifacts, validate_source_plan, validate_stage_figure, final_metadata
    print('Checking all-feature processing results and deriving compact display tables...',flush=True)
    reporting_artifacts(base,presentation/'publication_tables')
    if draw_figures:
        print('Rendering the single publication figure set...',flush=True)
        from manuscript_figures import PublicationFigures
        PublicationFigures(base,models,presentation,out/'Figures').run()
        sys.path.insert(0,str(ROOT/'workflows/tissue_survival'))
        from reporting import figures
        figures(survival,out/'Figures')
        wanted={f'Figure{i}' for i in range(1,5)}|{f'FigureS{i}' for i in range(1,25)}
        actual={p.stem for p in (out/'Figures').glob('*.pdf')}
        if actual!=wanted or len(list((out/'Figures').iterdir()))!=84:raise ValueError('Final figure set incomplete')
    validate_stage_figure(base,presentation)
    print('Validating and grouping publication source data...',flush=True)
    def resolve(source):
        if source.startswith('review/'):return models/source.removeprefix('review/')
        if source.startswith('publication/'):return presentation/'publication_tables'/source.removeprefix('publication/')
        if source.startswith('closeout/sex_interaction/'):return models/'sex_interaction'/Path(source).name
        if source.startswith('closeout/figure_sources/'):return presentation/'figure_sources'/Path(source).name
        return base/source
    public=read_csv_exact(base/'data/public/gse131617/sample_manifest.csv')
    items=[];index=[];metadata={};regression=[];dictionary_frames=[]
    sourceplan=records(ROOT/'config/publication_sources.csv')
    expected={r['record_id']:r for r in records(ROOT/'config/result_reference.csv')}
    def compare(key,frame):
        observed=fingerprint(frame);reference=expected.get(key)
        ok=reference is not None and all(str(observed[k])==reference[k] for k in observed)
        regression.append(dict(record_id=key,match=ok,**observed))
    for spec in sourceplan:
        path=resolve(spec['source']);key=spec['record_id']
        if not path.is_file():raise FileNotFoundError('Required publication source: '+spec['source'])
        if path.suffix=='.json':
            payload=path.read_text(encoding='utf-8')
            if re.search(PATH_PATTERN,payload):raise ValueError('Local path in public metadata')
            name='metadata/'+key+'.json';metadata[name]=payload
            index.append(dict(record_id=key,file='Supplementary_Data_1.zip/'+name,source=spec['source'],rows=1,
                              columns_json='[]',source_sha256=sha(path),scope=spec['reason']))
            continue
        frame=read_csv_exact(path);check_public(frame,spec['source'],public)
        compare(key,frame)
        dictionary_frames.append((key,frame))
        items.append(dict(**spec,data=frame,sha256=sha(path)))
    groups=group_tables(items);n_by_family={}
    with zipfile.ZipFile(out/'Supplementary_Data_1.zip','x',zipfile.ZIP_DEFLATED) as z:
        for group in groups:
            kind=group['family'];n_by_family[kind]=n_by_family.get(kind,0)+1
            name=f'data/{kind}_{n_by_family[kind]:02d}.csv';parts=[]
            for item in group['items']:
                data=item['data'].copy();columns=list(data)
                data.insert(0,'_source_row',np.arange(1,len(data)+1).astype(str))
                data.insert(0,'_record_id',item['record_id']);parts.append(data)
                index.append(dict(record_id=item['record_id'],file='Supplementary_Data_1.zip/'+name,
                    source=item['source'],rows=len(data),columns_json=json.dumps(columns),
                    source_sha256=item['sha256'],scope=item['reason']))
            merged=pd.concat(parts,ignore_index=True).fillna('')
            for part in parts:
                if len(part):
                    restored=merged.loc[merged._record_id.eq(part._record_id.iloc[0]),part.columns].reset_index(drop=True)
                    if not restored.equals(part.reset_index(drop=True)):raise ValueError('Grouped source changed')
            z.writestr(name,csv_bytes(merged))
        for name,payload in metadata.items():z.writestr(name,payload)
    extra=[]
    for name in SURVIVAL:extra.append((name,read_csv_exact(survival/(name+'.csv'))))
    for name in ['Observation_Effects','Observation_Diagnostics']:
        d=read_csv_exact(observation/(name+'.csv'))
        d=d.drop(columns=['interpretation','inference_policy','source_status','diagnostic_only','note'],errors='ignore')
        extra.append((name,d))
    with zipfile.ZipFile(out/'Supplementary_Data_3.zip','x',zipfile.ZIP_DEFLATED) as z:
        for name,data in extra:
            key='D3_'+name;check_public(data,name);blob=csv_bytes(data)
            compare(key,data);dictionary_frames.append((key,data))
            z.writestr(name+'.csv',blob)
            index.append(dict(record_id=key,file='Supplementary_Data_3.zip/'+name+'.csv',
                source='Supplementary_Data_3/'+name+'.csv',rows=len(data),columns_json=json.dumps(list(data)),
                source_sha256=hashlib.sha256(blob).hexdigest(),scope='Final aggregate result'))
    feature=read_csv_exact(base/'targeted/qc_joint/features_baseline.csv')
    check_public(feature,'features_baseline');compare('Supplementary_Data_2',feature)
    feature.to_csv(out/'Supplementary_Data_2.csv',index=False)
    dictionary_frames.append(('Supplementary_Data_2',feature))
    from manuscript_bindings import table_bindings, retained_counts
    tm=table_bindings(ROOT/'config/submission_table_bindings.csv',{r['source']:r['record_id'] for r in index})
    fm=records(ROOT/'config/figure_sources.csv')
    ids={r['record_id'] for r in index}|{'Supplementary_Data_2'}
    if len(ids)!=len(index)+1:raise ValueError('Duplicate source ID')
    for row in tm+fm:
        if row['record_id'] not in ids:raise ValueError('Unresolved table/figure: '+str(row))
    required_citations=validate_source_plan(ROOT,ids)
    pd.DataFrame(index).to_csv(out/'SOURCE_INDEX.csv',index=False)
    for file,rows in [('TABLE_SOURCE_MAP.csv',tm),('FIGURE_SOURCE_MAP.csv',fm)]:
        pd.DataFrame(rows).to_csv(out/file,index=False)
    # Keep indexes once at the top level, not repeated in every archive.
    from used_dictionary import build
    print('Writing the used-variable dictionary and final integrity checks...',flush=True)
    count=build(base,survival,dictionary_frames,out/'Variable_Dictionary.csv')
    pd.DataFrame(regression).to_csv(audit/'MANUSCRIPT_COMPARISON.csv',index=False)
    counts=retained_counts(base)
    write_json(audit/'SOURCE_COUNTS.json',[dict(check=k,expected=e,observed=v,match=e==v) for k,e,v in counts])
    report=dict(manuscript_match=all(r['match'] for r in regression) and all(e==v for _,e,v in counts),checked_result_tables=len(regression),
                source_records=len(index),data1_csv_files=len(groups),data3_csv_files=len(extra),
                dictionary_rows=count,figures=28 if draw_figures else 0,required_source_citations=required_citations,
                source_validation_passed=True,document_checked=False)
    (out/'READ_FIRST.txt').write_text(
        'Figures 1-4 and S1-S24: one set in PDF, PNG and SVG.\n'
        'Data 1: aggregate table/figure sources, grouped only when schemas are compatible.\n'
        'Data 2: all QC-retained baseline expression-feature tests, including nonsignificant tests.\n'
        'Data 3: survival, single-nucleus and observation-process results.\n'
        'SOURCE_INDEX: use _record_id and _source_row to recover each source; columns_json records its columns.\n'
        'Variable_Dictionary: used input/model variables and exported result fields only.\n'
        'PUBLICATION_SUMMARY: actual bootstrap counts, dictionary size and single-nucleus status definitions.\n'
        'SUBMISSION_SHA256SUMS: file integrity for this single publication set; do not mix runs.\n'
        'Word documents and the response letter are not updated or certified by this program.\n'
        'No participant-level ADNI records. The run audit folder is local-only.\n',encoding='utf-8')
    final_metadata(base,survival,out,ROOT,len(index),dict(required_citations=required_citations,
        processing_models=7,full_feature_BH_verified=True,table4_rounding_verified=True,
        FigureS8_TableS5D_same_source=True,reference_match=report['manuscript_match']))
    write_json(audit/'EXPORT_CHECK.json',report)
    return report
