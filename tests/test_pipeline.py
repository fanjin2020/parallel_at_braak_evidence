"""Small synthetic regression checks; no clinical records or network access."""
from pathlib import Path
import json
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from pipeline_runtime import Runner, signature, sha, publish, recover_publication, write_json, inventory
from publication_export import group_tables, fingerprint
from observation_statistics import dispersion, rule_value
import pandas as pd


class PipelineTests(unittest.TestCase):
    def test_relocated_cache_files_are_protected(self):
        from offline_sources import resolved_inputs
        from public_sources import SOURCES
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);cache=root/'cache';cache.mkdir();external=root/'external';external.mkdir()
            paths={}
            for key,relative,_ in SOURCES:
                p=external/(key+'.input');p.write_text('synthetic');paths[key]=str(p)
            ontology=external/'ontology.obo';ontology.write_text('synthetic');paths['go_ontology']=str(ontology)
            cel=external/'CEL';cel.mkdir();paths['cel_dir']=str(cel)
            write_json(cache/'public_sources_manifest.json',{'paths':paths})
            actual=resolved_inputs(cache)
            self.assertIn(cel.resolve(),actual)
            self.assertIn(ontology.resolve(),actual)
            self.assertEqual(len(actual),len(SOURCES)+2)

    def test_publication_commit_and_tamper(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);stage=root/'.stage';stage.mkdir()
            (stage/'data.csv').write_text('x\n1\n')
            publish(stage,root,'submission',{'manuscript_match':True})
            self.assertTrue(recover_publication(root)['completed'])
            (root/'submission/data.csv').write_text('x\n2\n')
            with self.assertRaises(ValueError):recover_publication(root)

    def test_interrupted_commit_is_recovered(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);final=root/'results_to_review';final.mkdir()
            (final/'data.csv').write_text('x\n1\n')
            write_json(root/'audit/PUBLISH.json',dict(completed=True,output_directory=final.name,outputs=inventory(final)))
            self.assertIsNotNone(recover_publication(root))
            self.assertTrue((root/'COMPLETED.json').is_file())

    def test_signature_rejects_changes(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'sig.json'
            signature(p,{'input':'a'})
            signature(p,{'input':'a'},True)
            with self.assertRaises(ValueError):signature(p,{'input':'b'},True)
            with self.assertRaises(ValueError):signature(p,{'input':'a'},False)

    def test_completed_job_is_not_rerun(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);out=root/'fit';out.mkdir()
            r=Runner(root/'audit',root)
            cmd=[sys.executable,'-c','from pathlib import Path;p=Path("fit/value.txt");assert not p.exists();p.write_text("123")']
            r.run('fit',cmd,out)
            r.run('fit',cmd,out)
            self.assertEqual((out/'value.txt').read_text(),'123')

    def test_tampered_checkpoint_stops(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);out=root/'fit';out.mkdir()
            r=Runner(root/'audit',root)
            cmd=[sys.executable,'-c','from pathlib import Path;Path("fit/value.txt").write_text("123")']
            r.run('fit',cmd,out);(out/'value.txt').write_text('changed')
            with self.assertRaises(ValueError):r.run('fit',cmd,out)

    def test_failed_job_can_be_retried(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);out=root/'fit';out.mkdir()
            r=Runner(root/'audit',root)
            cmd=[sys.executable,'-c','from pathlib import Path;import sys;p=Path("fit/attempt");old=p.exists();p.write_text("retry");sys.exit(0 if old else 2)']
            with self.assertRaises(RuntimeError):r.run('fit',cmd,out)
            r.run('fit',cmd,out)
            self.assertEqual(json.loads((root/'audit/fit.json').read_text())['status'],'complete')

    def test_grouping_column_cap(self):
        frames=[dict(source='revision/clinical/x.csv',data=pd.DataFrame(columns=[f'x{i}' for i in range(40)])),
                dict(source='revision/clinical/y.csv',data=pd.DataFrame(columns=[f'y{i}' for i in range(40)]))]
        self.assertEqual(len(group_tables(frames)),2)

    def test_reference_detects_numeric_changes(self):
        a=pd.DataFrame({'x':['0.125','2'],'model':['a','b']});b=a.copy();b.loc[0,'x']='0.126'
        self.assertNotEqual(fingerprint(a)['numeric_sha256'],fingerprint(b)['numeric_sha256'])

    def test_reference_detects_label_changes(self):
        a=pd.DataFrame({'x':['.1'],'contrast':['male - female']});b=a.copy();b.loc[0,'contrast']='female - male'
        self.assertNotEqual(fingerprint(a)['text_sha256'],fingerprint(b)['text_sha256'])

    def test_constant_balance_not_zero_evidence(self):
        self.assertFalse(dispersion([1.,1.,1.])['smd_applicable'])
        self.assertTrue(dispersion([1.,2.,3.])['smd_applicable'])
        self.assertIsNone(rule_value(float('nan')))
        with self.assertRaises(ValueError):rule_value('unknown')


if __name__=='__main__':unittest.main()
