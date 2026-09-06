"""Actual-support design and independent Gaussian score checks (no outcomes fit)."""
import importlib.util
import os
from pathlib import Path
import sys
import unittest
import numpy as np
from scipy.stats import multivariate_normal
from sklearn.metrics import adjusted_rand_score
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts/analysis'))
import run_hypervolume_actual_support_blocked as m


class Algebra(unittest.TestCase):
    def setUp(self):
        self.rng=np.random.default_rng(123)
        self.x=self.rng.normal(size=(2,40,20,3))
        self.y=self.rng.random((2,40,20))>.5

    def test_fast_predictions_match_direct_mahalanobis(self):
        mu,cov,valid=m.base.fit_classes(self.x,self.y)
        pred=m.predict_qda(self.x[:,20:],mu[:,:20],cov[:,:20])
        scores=[]
        for k in (0,1):
            diff=self.x[:,None,20:]-mu[:,:20,k,None,None]
            q=np.einsum('btnpi,btij,btnpj->btnp',diff,np.linalg.inv(cov[:,:20,k]),diff)
            scores.append(-.5*(q+np.linalg.slogdet(cov[:,:20,k])[1][:,:,None,None]))
        np.testing.assert_array_equal(pred,scores[1]>scores[0])

    def test_fast_score_matches_independent_scipy_sklearn(self):
        x=self.x[:1];y=self.y[:1]
        mu,cov,valid=m.base.fit_classes(x,y)
        total=0.;nonconstant=0
        for i in range(20):
            for j in range(20,40):
                a=multivariate_normal.logpdf(x[0,j],mean=mu[0,i,0],cov=cov[0,i,0])
                b=multivariate_normal.logpdf(x[0,j],mean=mu[0,i,1],cov=cov[0,i,1])
                pred=b>a
                if valid[0,i] and valid[0,j] and pred.any() and (~pred).any():
                    total+=adjusted_rand_score(y[0,j],pred);nonconstant+=1
        got=m.score_domain(x,y)
        self.assertAlmostEqual(float(got['transfer'][0]),total/400,places=12)
        self.assertAlmostEqual(float(got['informative_pair_fraction'][0]),nonconstant/400,places=12)

    def test_parent_score_equivalence(self):
        old=m.base.N_TRAIN
        try:
            m.base.N_TRAIN=20
            original=m.base.score_domain(self.x,self.y)
        finally:m.base.N_TRAIN=old
        new=m.score_domain(self.x,self.y)
        np.testing.assert_allclose(new['transfer'],original['transfer'],atol=1e-12,rtol=0)
        np.testing.assert_allclose(new['separation'],original['separation'],atol=1e-12,rtol=0)

    def test_label_flip_invariance(self):
        y=self.y.copy();y[:,::2]=~y[:,::2]
        np.testing.assert_allclose(m.score_domain(self.x,self.y)['transfer'],m.score_domain(self.x,y)['transfer'],atol=1e-12,rtol=0)

    def test_invalid_species_not_dropped(self):
        y=np.zeros_like(self.y)
        got=m.score_domain(self.x,y)
        for k in got:np.testing.assert_array_equal(got[k],0)

    def test_heldout_labels_cannot_change_predictions(self):
        mu,cov,_=m.base.fit_classes(self.x[:,:20],self.y[:,:20])
        a=m.predict_qda(self.x[:,20:],mu,cov)
        yy=self.y.copy();yy[:,20:]=~yy[:,20:]
        mu2,cov2,_=m.base.fit_classes(self.x[:,:20],yy[:,:20])
        np.testing.assert_array_equal(a,m.predict_qda(self.x[:,20:],mu2,cov2))

    def test_constant_prediction_is_zero_information(self):
        y=np.array([False,True,False,True]);p=np.ones(4,dtype=bool)
        self.assertEqual(float(m.base.adjusted_rand_binary(y,p)),0)

    def test_sector_wrap_and_boundaries(self):
        np.testing.assert_array_equal(m.sector_id(np.array([-180,-90,0,90,180,179.99])),[0,1,2,3,0,3])

    def test_spherical_distance(self):
        xyz=m.xyz_of(np.array([0,0]),np.array([0,90]))
        self.assertAlmostEqual(float(m.chord_to_km(np.linalg.norm(xyz[0]-xyz[1]))),np.pi*m.RADIUS_KM/2,places=8)

    def test_seed_and_stages_distinct(self):
        self.assertEqual(m.seed_for('x',1),m.seed_for('x',1))
        self.assertNotEqual(m.seed_for('x',1),m.seed_for('y',1))

    def test_scenario_counts(self):
        self.assertEqual(len(m.NUISANCE),10);self.assertEqual(len(m.POSITIVE),6)
        self.assertEqual((10+16)*m.N_REPS*16,104000)


class ActualGeometry(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path=os.environ.get('HYPERVOLUME_GEOMETRY')
        if not path:raise unittest.SkipTest('set HYPERVOLUME_GEOMETRY for actual geometry checks')
        cls.geo=m.Geometry(Path(path));cls.s=cls.geo.schedule('unit_test',0)

    def test_source_and_capacity(self):
        self.assertEqual(len(self.geo.df),36900)
        self.assertEqual(int(self.geo.mask.sum()),21424)
        self.assertEqual([r['evaluation_species'] for r in self.geo.capacity],[60,62,59,26])
        self.assertEqual(len(self.geo.train_sid&self.geo.test_sid),0)

    def test_species_roles_and_pairing(self):
        sid=self.geo.sid[self.s[:,:,:,0]]
        for d in range(4):
            for f in range(4):
                self.assertEqual(len(set(sid[d,f,:20])),20)
                self.assertEqual(len(set(sid[d,f,20:])),20)
                self.assertTrue(set(sid[d,f,:20])<=self.geo.train_sid)
                self.assertTrue(set(sid[d,f,20:])<=self.geo.test_sid)
                np.testing.assert_array_equal(sid[d,f],sid[0,f])

    def test_twenty_distinct_photos(self):
        for row in self.s.reshape(-1,20):self.assertEqual(len(set(row)),20)

    def test_mask_and_region_holdout(self):
        self.assertTrue(self.geo.mask[self.s[:2]].all())
        for d in range(4):
            for f in range(4):
                self.assertTrue((self.geo.sector[self.s[d,f,20:]]==f).all())
                if d%2==1:self.assertFalse((self.geo.sector[self.s[d,f,:20]]==f).any())

    def test_buffer_distance_direct(self):
        for d in (1,3):
            for f in range(4):
                tr=self.geo.xyz[self.s[d,f,:20].ravel()];te=self.geo.xyz[self.s[d,f,20:].ravel()]
                distance=m.chord_to_km(np.linalg.norm(tr[:,None]-te[None,:],axis=-1))
                self.assertGreaterEqual(distance.min(),500.-1e-8)

    def test_targets_paired_between_training_modes(self):
        np.testing.assert_array_equal(self.s[0,:,20:],self.s[1,:,20:])
        np.testing.assert_array_equal(self.s[2,:,20:],self.s[3,:,20:])

    def test_reproducible_schedule_and_independent_stages(self):
        np.testing.assert_array_equal(self.s,self.geo.schedule('unit_test',0))
        self.assertFalse(np.array_equal(self.s,self.geo.schedule('other_unit_test',0)))

    def test_coherent_photo_labels(self):
        y=m.labels_for(self.geo,'unit_test','environmental_shared',1.,.5,0)
        yy=m.labels_for(self.geo,'unit_test','environmental_shared',1.,.5,0)
        np.testing.assert_array_equal(y,yy)
        indexed=y[self.s]
        for d in (0,2):np.testing.assert_array_equal(indexed[d,:,20:],indexed[d+1,:,20:])
        self.assertEqual(y.shape,(36900,))

if __name__=='__main__':unittest.main()
