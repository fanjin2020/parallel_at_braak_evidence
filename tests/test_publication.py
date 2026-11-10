"""Publication regression and export-only checks using synthetic inputs."""
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import os
import json
import sys
import tempfile
import unittest
import zipfile
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from publication_validation import (bh,processing_counts,PROCESSING_MODELS,display_q,
    table4_display,validate_source_plan,validate_stage_figure)
from prepare_reporting_tables import source_definitions
from export_existing import locate_run,run as export_only,existing_sources
from pipeline_runtime import inventory,sha,write_json


def processing_fixture(root):
    qc=root/'cel_qc';qc.mkdir(parents=True)
    target=root/'targeted/qc_joint';target.mkdir(parents=True)
    p=np.linspace(.0001,.9,22011)
    features=pd.DataFrame(dict(feature_id=[str(i) for i in range(22011)],log2FC=.1,
        ci_low=-.2,ci_high=.4,p=p,q=bh(p)))
    features.rename(columns={'feature_id':'id'}).to_csv(target/'features_baseline.csv',index=False)
    cohorts=[];states=[];targets=[]
    for i,model in enumerate(PROCESSING_MODELS):
        features.to_csv(qc/('features_'+model+'.csv'),index=False)
        n=[213,212,212,197,197,210,210][i]
        cohorts.append(dict(model=model,n_samples=n,n_donors=71 if n>=212 else 70))
        states.append(dict(component=model,status='completed',reason=''))
        t=features.head(17).copy();t['model']=model;t['n_samples']=n;t['n_donors']=cohorts[-1]['n_donors']
        targets.append(t)
    pd.DataFrame(cohorts).to_csv(qc/'model_cohorts.csv',index=False)
    pd.DataFrame(states).to_csv(qc/'analysis_status.csv',index=False)
    pd.concat(targets).to_csv(qc/'target_processing_comparison.csv',index=False)
    return qc


def completed_fixture(root):
    work=root/'.work'
    for name in ['models','presentation','observation','survival/results','primary/.work/targeted/qc_joint']:
        (work/name).mkdir(parents=True)
    base=work/'primary/.work'
    pd.DataFrame({'id':['test'],'p':[.5],'q':[.5]}).to_csv(base/'targeted/qc_joint/features_baseline.csv',index=False)
    (base/'source.csv').write_text('value\n1\n')
    public=root/'results_to_review';public.mkdir()
    pd.DataFrame([dict(record_id='test',file='Supplementary_Data_1.zip/data/test.csv',source='source.csv',
        source_sha256=sha(base/'source.csv'))]).to_csv(public/'SOURCE_INDEX.csv',index=False)
    (public/'Supplementary_Data_2.csv').write_bytes((base/'targeted/qc_joint/features_baseline.csv').read_bytes())
    with zipfile.ZipFile(public/'Supplementary_Data_3.zip','w'):pass
    write_json(root/'COMPLETED.json',dict(completed=True,output_directory=public.name,outputs=inventory(public)))
    return base


class PublicationTests(unittest.TestCase):
    def test_public_dictionary_metadata_and_provenance(self):
        path=ROOT/'docs/Variable_Dictionary.csv'
        record=json.loads((ROOT/'docs/Variable_Dictionary.provenance.json').read_text(encoding='utf-8'))
        fields=pd.read_csv(path,dtype=str).fillna('')
        self.assertEqual(sha(path),record['sha256'])
        self.assertFalse(record['contains_participant_observations'])
        self.assertEqual(list(fields.columns),['scope','field','role','definition','unit','coding_or_missing','definition_status','code'])
        self.assertEqual(len(fields),record['rows'])
        self.assertEqual(fields.scope.nunique(),record['scopes'])
        self.assertFalse(fields.duplicated(['scope','field']).any())
        self.assertEqual(fields.role.value_counts().to_dict(),record['roles'])
        self.assertEqual(fields.definition_status.value_counts().to_dict(),record['definition_status'])
        manifest=pd.read_csv(ROOT/'CODE_SHA256SUMS.csv',dtype=str)
        for relative in ['docs/Variable_Dictionary.csv','docs/Variable_Dictionary.provenance.json']:
            entry=manifest.loc[manifest.file.eq(relative)]
            self.assertEqual(len(entry),1)
            self.assertEqual(entry.iloc[0].sha256,sha(ROOT/relative))

    def test_correct_rounding_and_nonzero_small_q(self):
        self.assertEqual(display_q(.079159461995531),'0.0792')
        self.assertEqual(display_q(.0258516662887989),'0.0259')
        self.assertNotEqual(display_q(1e-10),'0.0000')
        with self.assertRaises(ValueError):display_q(float('nan'))

    def test_table4_display_preserves_inputs(self):
        rows=[]
        for gene in ['CD14','MTHFD2']+[f'GENE{i}' for i in range(14)]:
            for model in ['baseline','scan_date','joint_CA_joint_quantile','joint_NG_joint_quantile']:
                q=.079159461995531 if gene=='CD14' else .0258516662887989
                rows.append(dict(gene_symbol=gene,model=model,entrez_id='123',log2FC=-.2,ci_low=-.3,ci_high=-.1,q=q))
        frame=pd.DataFrame(rows);before=frame.copy(deep=True)
        d=table4_display(frame).set_index('gene_symbol')
        self.assertEqual(d.loc['CD14','joint_CA_q'],'0.0792')
        self.assertEqual(d.loc['MTHFD2','joint_CA_q'],'0.0259')
        pd.testing.assert_frame_equal(frame,before)

    def test_full_feature_counts_not_selected_target_counts(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);processing_fixture(root)
            d=processing_counts(root)
            self.assertEqual(len(d),7)
            self.assertTrue(d.n_features.eq(22011).all())
            self.assertTrue(d.n_features_q_lt_0_05.eq(0).all())

    def test_missing_full_feature_file_stops(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);qc=processing_fixture(root)
            (qc/('features_'+PROCESSING_MODELS[0]+'.csv')).unlink()
            with self.assertRaisesRegex(FileNotFoundError,'full-feature evidence'):processing_counts(root)

    def test_bad_feature_q_stops(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);qc=processing_fixture(root);path=qc/('features_'+PROCESSING_MODELS[0]+'.csv')
            d=pd.read_csv(path);d.loc[100,'q']=.01;d.to_csv(path,index=False)
            with self.assertRaisesRegex(ValueError,'BH mismatch'):processing_counts(root)

    def test_source_list_includes_single_nucleus(self):
        sources=source_definitions().source.tolist()
        self.assertEqual(sources.count('GSE174367'),1)
        self.assertEqual(len(sources),7)

    def test_cited_sources_cannot_be_silently_pruned(self):
        d=pd.read_csv(ROOT/'config/manuscript_source_requirements.csv')
        ids=set(d.record_id)
        self.assertEqual(validate_source_plan(ROOT,ids),len(d))
        for key in ['SD238','SD240','CL_124','QC_PROCESSING_COUNTS']:
            with self.assertRaisesRegex(ValueError,key):validate_source_plan(ROOT,ids-{key})

    def test_stage_figure_is_from_same_bootstrap(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);base=root/'base';pres=root/'presentation'
            (base/'enhancement/clinical').mkdir(parents=True);(pres/'figure_sources').mkdir(parents=True)
            d=pd.DataFrame([dict(outcome='cdrsb',contrast='B_vs_A',requested=500,valid=499,estimate=.7,ci_low=.2,ci_high=1.3,status='ok')])
            d.to_csv(base/'enhancement/clinical/stage_propagated_intervals.csv',index=False)
            path=pres/'figure_sources/FigureS19_source.csv';d.to_csv(path,index=False)
            validate_stage_figure(base,pres)
            d['requested']=1000;d['valid']=998;d.to_csv(path,index=False)
            with self.assertRaisesRegex(ValueError,'Figure S8 source differs'):validate_stage_figure(base,pres)
            d['requested']=0;d['valid']=0
            d.to_csv(base/'enhancement/clinical/stage_propagated_intervals.csv',index=False)
            d.to_csv(path,index=False)
            with self.assertRaisesRegex(ValueError,'Invalid threshold'):validate_stage_figure(base,pres)

    def test_latest_complete_fit_run_not_export_only_run(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);runs=root/'runs';runs.mkdir()
            a=runs/'first';b=runs/'second';completed_fixture(a);completed_fixture(b)
            os.utime(a/'COMPLETED.json',(100,100));os.utime(b/'COMPLETED.json',(200,200))
            c=runs/'export';c.mkdir();write_json(c/'COMPLETED.json',{'completed':True})
            self.assertEqual(locate_run(root),b)
            self.assertEqual(locate_run(a),a)

    def test_export_only_never_starts_models_or_changes_old_run(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);previous=root/'runs/old';completed_fixture(previous)
            (root/'CODE_SHA256SUMS.csv').write_text('synthetic test manifest')
            before=inventory(previous)
            args=SimpleNamespace(resume=False,check=False,prepared_inputs=None,from_run=root,
                output_dir=root/'runs/new',keep_work=False)
            def fake_export(base,models,presentation,observation,survival,stage,audit):
                stage.mkdir();(stage/'result.csv').write_text('verified\n1\n')
                return dict(manuscript_match=True)
            with patch('publication_export.export',side_effect=fake_export), \
                 patch('subprocess.Popen',side_effect=AssertionError('Unexpected subprocess/model run')):
                export_only(args,root)
            self.assertEqual(inventory(previous),before)
            record=json.loads((root/'runs/new/COMPLETED.json').read_text())
            self.assertEqual(record['model_fits_started'],0)
            self.assertEqual(record['downloads_started'],0)
            self.assertTrue((root/'runs/new/submission/result.csv').is_file())

    def test_changed_old_source_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);base=completed_fixture(root)
            (base/'source.csv').write_text('value\n2\n')
            with self.assertRaisesRegex(ValueError,'Stored source changed'):existing_sources(root)


if __name__=='__main__':unittest.main()
