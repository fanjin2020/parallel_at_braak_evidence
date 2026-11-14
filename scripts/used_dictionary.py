"""Dictionary restricted to selected source/model variables and exported results."""
import json
from pathlib import Path
import sys
import pandas as pd
from result_dictionary import result_definition

ROOT=Path(__file__).resolve().parents[1]
COMMON=['RID','HAS_QC_ERROR']
PET=['RID','SCANDATE','LONIUID','TRACER','qc_flag','META_TEMPORAL_SUVR']
RAW={
 'UCBERKELEY_AMY_6MM':(['RID','SCANDATE','LONIUID','TRACER','qc_flag','AMYLOID_STATUS_COMPOSITE_REF','CENTILOIDS'],'prepare_adni.py; analyze_tau_incremental.py'),
 'UCBERKELEY_TAU_6MM':(PET+['SITEID','ORIGPROT'],'prepare_adni.py; pet_staging_covariates.py'),
 'UCBERKELEY_TAUPVC_6MM':(PET,'pet_staging_covariates.py'),
 'UCSFFSX7':(['RID','EXAMDATE','ST29SV','ST88SV','ST10CV','FIELD_STRENGTH','OVERALLQC','FSVER'],'prepare_adni.py; pet_staging_covariates.py; pet_day_clinical_models.py'),
 'CDR':(COMMON+['VISDATE','CDRSB'],'prepare_adni.py; longitudinal.py'),
 'DXSUM':(COMMON+['EXAMDATE','DIAGNOSIS'],'longitudinal.py; pet_day_clinical_models.py'),
 'PTDEMOG':(['RID','PTGENDER','PTEDUCAT','PTDOB','VISDATE'],'prepare_adni.py'),
 'APOERES':(['RID','GENOTYPE'],'prepare_adni.py'),
 'ADSL':(['SUBJID','RACE','ETHNIC'],'pet_staging_covariates.py'),
 'UCBERKELEYFDG_8mm':(['RID','EXAMDATE','ROINAME','MEAN'],'prepare_adni.py'),
 'GDSCALE':(COMMON+['VISDATE','GDTOTAL'],'prepare_confounders.py'),
 'BACKMEDS':(COMMON+['VISDATE','KEYMED'],'prepare_confounders.py'),
 'MEDHIST':(COMMON+['VISDATE','MHPSYCH','MH2NEURL','MH4CARD','MH5RESP','MH6HEPAT','MH7DERM','MH8MUSCL','MH9ENDO','MH10GAST','MH11HEMA','MH12RENA','MH17MALI'],'prepare_confounders.py'),
 'MODHACH':(COMMON+['VISDATE','HMSCORE','HMSTROKE'],'prepare_confounders.py'),
 'RECCMEDS':(COMMON+['VISDATE','CMMED','CMMEDO','CMBGN','CMEND','CMCONT'],'prepare_confounders.py'),
 'INITHEALTH':(COMMON+['VISDATE','IHPRESENT','IHSYMPTOM'],'prepare_confounders.py'),
}
EVENT={
 'STUDYSUM':COMMON+['COLPROT','SDPRIMARY','SDDATE'],
 'REGISTRY':COMMON+['COLPROT','VISTYPE','RGSOURCE','RGCONDCT','EXAMDATE','RGREASON'],
 'ADVERSE':COMMON+['COLPROT','SAEDEATH','AEOUTCOME','AEHDTHDT'],
 'RECADV':COMMON+['COLPROT','AEHDEATH','AEHDTHDT'],
 'TREATDIS':COMMON+['COLPROT','WDREASON','EXAMDATE'],
 'OTDISC':COMMON+['COLPROT','OTDATE','OTREASON'],
}


def build(base,survival,frames,dest):
    rules=json.loads((ROOT/'config/variable_definitions.json').read_text(encoding='utf-8'))['fields']
    snapshot=pd.read_csv(ROOT/'config/adni_source_dictionary_snapshot.csv',dtype=str).fillna('')
    official=dict(tuple(snapshot.groupby(['source_table','field'],sort=False)))
    rows=[]
    def add(scope,field,role,evidence,definition=None):
        known=definition or rules.get(field)
        if known and not str(known[0]).strip():known=None
        rows.append(dict(scope=scope,field=field,role=role,definition=known[0] if known else 'Definition not yet verified for this used field',
            unit=known[1] if known else 'not_inferred',coding_or_missing=known[2] if known else 'Do not infer coding from the field name',
            definition_status='documented' if known else 'source_definition_unverified',code=evidence))
    for table,(fields,code) in RAW.items():
        path=next((base/'data/private/adni'/folder/(table+'.csv') for folder in ['raw','confounder_raw']
                   if (base/'data/private/adni'/folder/(table+'.csv')).is_file()),None)
        if path is None:continue
        columns=list(pd.read_csv(path,nrows=0).columns)
        if table=='UCBERKELEY_TAU_6MM':
            # The actual code checks all regional measurements for same-scan ambiguity.
            fields=fields+[c for c in columns if c.endswith(('_SUVR','_VOLUME'))]
        for field in dict.fromkeys(fields):
            if field not in columns:continue
            d=official.get((table.upper(),field));meaning=None
            if d is not None:
                def values(c):return '; '.join(dict.fromkeys(d[c].loc[d[c].ne('')]))
                if values('definition'):meaning=[values('definition'),values('units'),values('coding')]
            if meaning is None and field.startswith('CTX_') and field.endswith(('_SUVR','_VOLUME')):
                meaning=['Regional tau measurement used in regional aggregation or exact-scan duplicate checks',
                         'SUVR' if field.endswith('_SUVR') else 'source volume units',
                         'Volume-weighted cortical composites; nonfinite or nonpositive weights are not replaced']
            add('ADNI/'+table,field,'source_input',code,meaning)
    event_path=survival/'ADNI_event_variable_dictionary.csv'
    if event_path.is_file():
        d=pd.read_csv(event_path,dtype=str).fillna('')
        for table,fields in EVENT.items():
            for field in fields:
                part=d.loc[d.TBLNAME.eq(table)&d.FLDNAME.eq(field)]
                if part.empty:continue
                def values(c):return '; '.join(dict.fromkeys(part[c])) if c in part else ''
                add('ADNI/'+table,field,'source_input','vital_followup.py',
                    [values('TEXT'),values('UNITS'),values('CODE')])
    # Only variables actually present in the model input, never all raw schemas.
    for path,scope,code in [
        (base/'data/private/adni/adni_pet_aligned_master.csv','ADNI/model','prepare_adni.py; adni_models.py'),
        (base/'cel_qc/retained_sample_manifest.csv','GSE131617/model','qc_joint_models.R; sex_interaction.R'),
        (base/'data/public/cel_metadata/sample_technical_metadata.csv','GSE131617/technical','source_io.py'),
        (survival/'GSE174367_donors.csv','GSE174367/model','morabito_models.R; celltype_interaction.R')]:
        if path.is_file():
            for field in pd.read_csv(path,nrows=0).columns:
                if field in rules:add(scope,field,'model_variable',code)
    # Public covariate definitions are tied to the actual preparation code and cohort.
    public_models=[
        (base/'data/public/GSE44770/manifest.csv','GSE44770/model','prepare_external_gse44770.py',{
            'sample_id':['GEO sample accession matched to expression columns','identifier','One sample per donor'],
            'donor_key':['Donor title from the GEO sample metadata','identifier','Duplicate donor keys stop preparation'],
            'group':['Disease group used in the external comparison','category','A becomes case; N becomes reference; other diagnoses excluded'],
            'age_at_death':['Age at death supplied by the cohort','years','Numeric conversion; missing values retained'],
            'sex_male':['Sex indicator derived from source gender','indicator','M=1; F=0'],
            'pmi_hours':['Postmortem interval supplied by the cohort','hours','Missing values retained'],
            'tissue_ph':['Brain tissue pH supplied by the cohort','pH','Missing values retained'],
            'rin':['RNA integrity number supplied by the cohort','RIN scale','Not inferred from scan date'],
            'batch':['Experimental batch supplied by the cohort','category','Source batch prefixed with batch_; NA remains missing']}),
        (base/'data/public/GSE106241/manifest.csv','GSE106241/model','prepare_braak_external.py',{
            'sample_id':['GEO accession matched to expression columns','identifier','Unique sample'],
            'donor_key':['Donor title matched to the published matrix columns','identifier','Unique donor'],
            'brain_region':['Brain region supplied by the cohort','category','One region required by preparation'],
            'braak_stage':['Published Braak stage','ordinal stage','0 through 6'],
            'braak_group':['Grouped Braak stage used in the model','category','0; I-II; III-IV; V-VI'],
            'age_at_death':['Age supplied in the postmortem metadata','years','Numeric conversion'],
            'sex_male':['Sex indicator from source sex coding','indicator','Source 1 female becomes 0; source 2 male becomes 1'],
            'apoe4_count':['Published additive APOE epsilon4 allele count','alleles','0, 1, or 2; missing retained'],
            'rna_quality_category':['Published categorical RNA quality','category','Kept categorical; not relabeled as numeric RIN']}),
        (survival/'GSE174367_donors.csv','GSE174367/model','morabito_pseudobulk.py',{
            'donor_id':['SampleID shared by nuclei from the same donor','identifier','Donor covariates must be constant across nuclei'],
            'age':['Source Age for the donor','years','Required for the primary model'],
            'Sex':['Source Sex for the donor','category','F or M'],
            'sex_male':['Indicator derived from source Sex','indicator','M=1; F=0'],
            'batch':['Source Batch for the donor','category','Required for the primary model'],
            'rin':['Source RIN for the donor','RIN scale','Required for the primary model'],
            'pmi':['Source PMI for the donor','source PMI units','Complete cases used in the PMI-adjusted sensitivity model'],
            'Tangle.Stage':['Published tangle stage string','category','Stage 0 through Stage 6'],
            'braak':['Integer parsed from Tangle.Stage','ordinal stage','I-II reference; V-VI case; no invented stage-0 donors']})]
    for path,scope,code,fields in public_models:
        if path.is_file():
            for field in pd.read_csv(path,nrows=0).columns:
                if field in fields:
                    rows[:]=[r for r in rows if not (r['scope']==scope and r['field']==field)]
                    add(scope,field,'model_variable',code,fields[field])
    for field,meaning in {
        'Barcode':['Cell barcode used to align metadata to the count matrix','identifier','Must be unique and present in the H5 matrix'],
        'SampleID':['Donor identifier used for pseudobulk aggregation','identifier','Not an independent-cell replicate'],
        'Cell.Type':['Author-assigned cell type used for aggregation','category','MG, ASC, ODC, OPC, EX, INH, PER.END'],
        'matrix/data':['Raw nonnegative integer counts in sparse CSC format','counts','Normalized expression is rejected'],
        'matrix/indices':['Zero-based gene-row positions for nonzero counts','index','Bounds checked against the number of features'],
        'matrix/indptr':['Sparse CSC column boundaries','index','Length, monotonicity and terminal count checked'],
        'matrix/barcodes':['Matrix cell barcodes matched to metadata Barcode','identifier','Metadata-matched cells enter aggregation'],
        'matrix/features/id':['Gene feature identifiers','identifier','Must be unique'],
        'matrix/features/name':['Gene symbols accompanying feature identifiers','gene symbol','Mapping follows the gene annotation rules'],
        'matrix/shape':['Declared matrix dimensions','features by cells','Checked against gene and barcode lengths']}.items():
        add('GSE174367/source',field,'source_input','morabito_pseudobulk.py',meaning)
    sys.path.insert(0,str(ROOT/'workflows/tissue_survival'))
    from result_schema import definition as sn_definition
    from observation_dictionary import observation_definition
    # The same field can have different units across model families: keep source identity.
    for key,frame in frames:
        for field in frame.columns:
            meaning=None
            if key=='QC_PROCESSING_COUNTS':
                meaning={
                    'model':['Array processing and QC analysis specification','label','Seven supplied/reprocessed and exclusion variants'],
                    'n_samples':['Number of arrays included in this processing fit','count','Taken from model_cohorts for the corresponding fitted model'],
                    'n_donors':['Number of independent tissue donors in this fit','count','Donor blocking retained'],
                    'n_features':['Number of tested microarray features','count','Entire 22,011-feature background, not selected targets'],
                    'n_features_q_lt_0_05':['Count of features with full-background BH q below 0.05','count','Counted from the complete fitted feature table using q < 0.05'],
                    'q_threshold':['Significance cutoff used only to count fitted features','probability','Strictly q < 0.05'],
                    'multiple_testing':['Multiplicity family and adjustment method','text','BH across all tested features within each processing model'],
                    'source_result':['Relative analysis-output path used for the count','relative path','Not a local machine path or a participant identifier'],
                }.get(field)
            elif key=='MAIN_TABLE4_DISPLAY':
                meaning={
                    'gene_symbol':['Selected non-Y target symbol','symbol','16 targets; UTY reported separately'],
                    'entrez_id':['Target NCBI Gene identifier','identifier','Fixed feature-to-gene mapping'],
                    'baseline_log2FC':['QC baseline expression contrast formatted to three decimals','log2 expression','Display only; full precision remains in the QC target source'],
                    'baseline_95CI':['QC baseline nominal 95% confidence limits formatted to three decimals','log2 expression','Display only; not a new confidence-interval calculation'],
                    'baseline_q':['QC baseline full-feature BH q formatted for Table 4','probability','Four decimals; values below 0.0001 use scientific notation'],
                    'joint_CA_q':['Joint scan-date and CA-composition full-feature BH q','probability','Four decimals; values below 0.0001 use scientific notation'],
                    'joint_NG_q':['Joint scan-date and NG-composition full-feature BH q','probability','Four decimals; values below 0.0001 use scientific notation'],
                }.get(field)
            elif key.startswith('D3_GSE') or key.startswith('D3_MNAR'):
                meaning=sn_definition(key.removeprefix('D3_')+'.csv',field)
            elif key.startswith('D3_Observation'):
                meaning=[observation_definition(field),'see model specification','See status and interval columns']
            elif field in rules:meaning=rules[field]
            if key=='CL_119' and field.endswith(('_effect','_ci_low','_ci_high')):
                meaning=['Braak/sex contrast or confidence bound in standardized GO-score expression','GO-score units','Not log2 fold change']
            meaning=result_definition(key,field) or meaning
            add(key,field,'exported_result','SOURCE_INDEX.csv',meaning)
    data=pd.DataFrame(rows).drop_duplicates()
    data.to_csv(dest,index=False)
    return len(data)
