"""Synthetic endpoint, denominator, isolation and survival regression guards."""
import sys
from pathlib import Path
import tempfile
import unittest
import numpy as np
import pandas as pd

MODULE=Path(__file__).resolve().parents[2]/'src/validation/manuscript_2026'
sys.path.insert(0,str(MODULE))
from prepare import landmark, median_groups, normalize_scores, COMPONENTS, unique
from verify import survival_rows, reml, logrank
from common import new_directory, REPOSITORY

class ValidationGuards(unittest.TestCase):
    def test_six_month_boundary_and_early_censoring(self):
        actual=landmark([5,5,6,6,7,7,-1,np.nan,3],[1,0,1,0,1,0,1,1,2])
        np.testing.assert_allclose(actual,[1,np.nan,1,0,0,0,np.nan,np.nan,np.nan],equal_nan=True)

    def test_zero_median_unequal_groups(self):
        cutoff, high, groups=median_groups(pd.Series([0,0,0,0,.01,.03]))
        self.assertEqual(cutoff,0)
        self.assertEqual(high.sum(),2)
        self.assertEqual(list(groups),['Low']*4+['High']*2)

    def scores(self):
        return pd.DataFrame({'Mixture':['synthetic_A','synthetic_B'],**{c:[2.,1.] for c in COMPONENTS},'Absolute score (sig.score)':[10.,5.]})

    def test_all_components_in_denominator(self):
        result=normalize_scores(self.scores())
        np.testing.assert_allclose(result,.2)
        np.testing.assert_allclose(result.sum(axis=1),1)

    def test_invalid_component_scores_and_duplicate_samples(self):
        for mode in ['negative','zero','nan','mismatch','duplicate']:
            d=self.scores()
            if mode=='negative': d.loc[0,COMPONENTS[0]]=-1
            if mode=='zero': d.loc[0,COMPONENTS]=0
            if mode=='nan': d.loc[0,COMPONENTS[0]]=np.nan
            if mode=='mismatch': d.loc[0,'Absolute score (sig.score)']=20
            if mode=='duplicate': d.loc[1,'Mixture']='synthetic_A'
            with self.subTest(mode=mode),self.assertRaises(ValueError): normalize_scores(d)

    def test_patient_uniqueness(self):
        with self.assertRaises(ValueError): unique(pd.DataFrame({'id':['synthetic','synthetic']}),'id','test')

    def test_outputs_cannot_be_created_in_public_checkout(self):
        with self.assertRaises(ValueError): new_directory(REPOSITORY/'should_not_exist_patient_data')

    def test_outputs_never_overwrite(self):
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaises(ValueError): new_directory(temp)

    def test_km_tied_event_and_censor(self):
        d=pd.DataFrame({'time_months':[1,1,2,3],'event':[1,0,1,0]})
        actual=survival_rows(d)
        np.testing.assert_allclose(actual[:,[0,1,4,5,6]],[[1,.75,4,1,1],[2,.375,2,1,0],[3,.375,1,0,1]])

    def test_reml_boundary_and_direction(self):
        result=reml([.2,.2,.2],[.3,.3,.3])
        self.assertEqual(result['tau2'],0)
        self.assertAlmostEqual(result['beta'],.2)
        reverse=reml([-.2,-.2,-.2],[.3,.3,.3])
        self.assertAlmostEqual(result['beta'],-reverse['beta'])
        self.assertAlmostEqual(result['p'],reverse['p'])

    def test_logrank_direction_invariant(self):
        d=pd.DataFrame({'time_months':[1,2,2,3,4,5],'event':[1,1,0,1,1,0],'high':[1,1,0,1,0,0]})
        result=logrank(d);d['high']=1-d.high
        self.assertAlmostEqual(result,logrank(d))

if __name__=='__main__':
    unittest.main()
