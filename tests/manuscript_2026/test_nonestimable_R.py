"""Synthetic integration check: never pool silently when a cohort has zero IQR."""
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
import pandas as pd

MODULE=Path(__file__).resolve().parents[2]/'src/validation/manuscript_2026'
COMPONENTS=['CCL8+CCL13+ Myeloid DC','Myeloid DC','B cells','T cells','Tregs']

@unittest.skipUnless(shutil.which('Rscript'), 'Rscript required for integration guard')
class NonestimableModel(unittest.TestCase):
    def test_zero_IQR_writes_failure_and_no_pool(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);inp=root/'input';inp.mkdir();out=root/'out'
            sample=[f'synthetic_{i:02}' for i in range(24)]
            scores=pd.DataFrame({'Mixture':sample,**{c:[1.]*24 for c in COMPONENTS},'Absolute score (sig.score)':[5.]*24})
            scores.to_csv(inp/'cologne_absolute_scores_LOCAL_ONLY.csv',index=False)
            shares=pd.DataFrame({'source_id':sample,**{c:[.2]*24 for c in COMPONENTS}})
            shares.to_csv(inp/'cologne_normalized_shares.csv',index=False)
            response=[]
            for cohort in ['moffitt','zuma7','stahl','li_cart']:
                for i in range(22):
                    response.append(dict(cohort=cohort,endpoint='synthetic',marker='CCL8_CCL13',source_id=sample[i],patient_id=sample[i],fraction=.2,bad=i%2))
            pd.DataFrame(response).to_csv(inp/'response_analysis_data.csv',index=False)
            chemo=[]
            for cohort in ['monti','rchop','chop']:
                for i in range(8):
                    chemo.append(dict(cohort=cohort,source_id=sample[i],patient_id=sample[i],OS_Time=i+1,OS_Status=i%2,fraction=.2))
            pd.DataFrame(chemo).to_csv(inp/'chemo_analysis_data.csv',index=False)
            survival=[]
            for marker in ['CCL8','CCL13','CCL8_CCL13','CCL18']:
                survival.append(dict(cohort='li_cart',endpoint='EFS',marker=marker,source_id='synthetic',patient_id='synthetic',time_months=1,event=1,fraction=.2,high=0,group='Low'))
            pd.DataFrame(survival).to_csv(inp/'survival_analysis_data.csv',index=False)
            process=subprocess.run(['Rscript',str(MODULE/'analyze.R'),str(inp),str(out)],capture_output=True,text=True)
            self.assertNotEqual(process.returncode,0)
            self.assertTrue((out/'model_statistics.csv').is_file(),process.stderr)
            models=pd.read_csv(out/'model_statistics.csv')
            self.assertEqual(len(models),7)
            self.assertFalse(models.estimable.any())
            self.assertTrue(models.reason.str.contains('IQR = 0').all())
            self.assertEqual(json.loads((out/'R_validation.json').read_text())['status'],'FAIL')
            self.assertFalse((out/'pooled_statistics.csv').exists())

if __name__=='__main__':
    unittest.main()
