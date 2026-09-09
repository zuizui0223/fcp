"""Independent probability/algebra checks for the pooled-response extension."""
from pathlib import Path
import sys
import itertools
import unittest
import numpy as np
from scipy.special import ndtr, expit, logsumexp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/analysis'))
import run_hypervolume_pooled_response as m

class ProbabilityTests(unittest.TestCase):
    def setUp(self):
        self.rng=np.random.default_rng(77109)
        self.x=self.rng.normal(size=(2,4,6,2))
        self.y=self.rng.integers(0,2,size=(3,2,4,6)).astype(bool)

    def slow_likelihood(self,x,y,link,axis):
        terms=[]
        projection=x@m.AXES[axis]
        for a in m.SLOPES[link]:
            p=ndtr(a*np.tanh(projection/.5)) if link=='saturating_probit' else expit(a*projection)
            for sign in (1,-1):
                prob=p if sign==1 else 1-p
                terms.append(np.prod(np.where(y,prob,1-prob)))
        return np.log(sum(terms)/len(terms))

    def test_saturating_likelihood_matches_explicit_products(self):
        ll=m.axis_log_likelihood(self.y,m.likelihood_features(self.x,'saturating_probit'))
        for k in (0,11,35,47):
            self.assertAlmostEqual(ll[1,0,2,k],self.slow_likelihood(self.x[0,2],self.y[1,0,2],'saturating_probit',k),places=12)

    def test_logit_likelihood_matches_explicit_products(self):
        ll=m.axis_log_likelihood(self.y,m.likelihood_features(self.x,'linear_logit'))
        for k in (0,11,35,47):
            self.assertAlmostEqual(ll[1,0,2,k],self.slow_likelihood(self.x[0,2],self.y[1,0,2],'linear_logit',k),places=12)

    def test_whole_species_label_swap_invariant(self):
        flip=self.y.copy();flip[:,:,1]=~flip[:,:,1]
        for model in m.SLOPES:
            feat=m.likelihood_features(self.x,model)
            np.testing.assert_allclose(m.axis_log_likelihood(self.y,feat),m.axis_log_likelihood(flip,feat),atol=1e-12,rtol=0)

    def test_all_binary_data_probabilities_sum_to_one(self):
        x=self.x[:1,:1,:3]
        labels=np.array(list(itertools.product([False,True],repeat=3)))[:,None,None,:]
        for model in m.SLOPES:
            ll=m.axis_log_likelihood(labels,m.likelihood_features(x,model))
            np.testing.assert_allclose(np.exp(ll).sum(0),1,atol=1e-12,rtol=0)

    def test_zero_coordinates_yield_no_orientation_information(self):
        for model in m.SLOPES:
            ll=m.axis_log_likelihood(self.y,m.likelihood_features(np.zeros_like(self.x),model))
            np.testing.assert_allclose(m.normalized_log_evidence(ll),0,atol=1e-12,rtol=0)

    def test_evidence_axis_average_is_one(self):
        lr=m.normalized_log_evidence(self.rng.normal(size=(3,2,4,m.K)))
        np.testing.assert_allclose(np.exp(lr).mean(-1),1,atol=1e-12,rtol=0)

    def test_posterior_matches_explicit_enumeration(self):
        lr=m.normalized_log_evidence(self.rng.normal(size=(5,m.K)))
        post=m.train_posterior(lr)
        weights=np.array([[np.prod([1-f+f*np.exp(lr[s,k]) for s in range(5)]) for k in range(m.K)] for f in m.FRACTIONS])
        weights/=weights.sum()
        np.testing.assert_allclose(post,weights,atol=1e-13,rtol=1e-12)

    def test_predictive_gain_matches_explicit_enumeration(self):
        lr=m.normalized_log_evidence(self.rng.normal(size=(7,m.K)))
        post=m.train_posterior(lr[:5]);got=m.predictive_gain(post,lr[5:])
        expected=np.log([sum(post[fi,k]*(1-f+f*np.exp(lr[s,k])) for fi,f in enumerate(m.FRACTIONS) for k in range(m.K)) for s in (5,6)])
        np.testing.assert_allclose(got,expected,atol=1e-13,rtol=0)

    def test_uniform_training_cannot_manufacture_target_evidence(self):
        post=m.train_posterior(np.zeros((20,m.K)))
        lr=m.normalized_log_evidence(self.rng.normal(size=(20,m.K))*4)
        np.testing.assert_allclose(m.predictive_gain(post,lr),0,atol=1e-12,rtol=0)

    def test_target_labels_do_not_update_training_posterior(self):
        lr=m.normalized_log_evidence(self.rng.normal(size=(40,m.K)))
        post=m.train_posterior(lr[:20]);saved=post.copy()
        m.predictive_gain(post,lr[20:])
        m.predictive_gain(post,m.normalized_log_evidence(self.rng.normal(size=(20,m.K))*20))
        np.testing.assert_array_equal(post,saved)

    def test_zero_sharing_predictive_ratio_is_one(self):
        post=np.zeros((5,m.K));post[0]=1/m.K
        target=m.normalized_log_evidence(self.rng.normal(size=(7,m.K))*3)
        np.testing.assert_allclose(m.predictive_gain(post,target),0,atol=1e-13)

    def test_joint_predictive_density_normalized(self):
        # The predictive score is a per-target marginal density ratio, not a
        # joint likelihood treating repeated folds as independent.
        x=self.x[:1,:1,:3];labels=np.array(list(itertools.product([False,True],repeat=3)))[:,None,None,:]
        ll=m.axis_log_likelihood(labels,m.likelihood_features(x,'linear_logit'))[:,0,0]
        lr=m.normalized_log_evidence(ll)
        post=m.train_posterior(m.normalized_log_evidence(self.rng.normal(size=(20,m.K))))
        density0=np.exp(logsumexp(ll,axis=-1)-np.log(m.K))
        gain=m.predictive_gain(post,lr)
        self.assertAlmostEqual(float(np.sum(density0*np.exp(gain))),1.,places=12)

    def test_single_source_average_is_normalized(self):
        lr=m.normalized_log_evidence(self.rng.normal(size=(2,3,40,m.K)))
        score,pf=m.marginal_scores(lr)
        self.assertEqual(score['joint_train'].shape,(2,3))
        self.assertTrue(np.all((pf>=0)&(pf<=1)))
        sep=np.exp(lr[0,0,:20]).mean(0)/m.K
        expect=np.log(np.exp(lr[0,0,20:])@sep).mean()/20
        self.assertAlmostEqual(expect,score['single_source_average'][0,0],places=13)

    def test_invalid_arguments_rejected(self):
        with self.assertRaises(ValueError):m.likelihood_features(self.x,'unknown')
        with self.assertRaises(ValueError):m.likelihood_features(self.x*float('nan'),'linear_logit')
        with self.assertRaises(ValueError):m.axis_log_likelihood(self.y.astype(float)+1,m.likelihood_features(self.x,'linear_logit'))
        with self.assertRaises(ValueError):m.normalized_log_evidence(np.zeros((3,7)))

    def test_seed_and_predeclared_designs(self):
        self.assertEqual(m.SEED,2026090704);self.assertEqual(m.REPS,250)
        self.assertEqual(len(m.MODELS),5);self.assertEqual(len(m.NUISANCE),10)
        self.assertEqual(len(m.POSITIVE),4)
        self.assertNotEqual(m.seed_for('calibration',0),m.seed_for('evaluation',0))

    def test_source_hashes(self):
        for name,h in m.PARENT_HASHES.items():
            self.assertEqual(m.digest(Path(m.__file__).with_name(name)),h)

class ActualGeometryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path=Path(__file__).resolve().parents[1]/'inputs/geometry_only.csv'
        if not path.exists():raise unittest.SkipTest('coordinate export not present')
        cls.geo=m.support.Geometry(path)

    def test_fixed_capacity(self):
        self.assertEqual(int(self.geo.mask.sum()),21424)
        self.assertEqual(len(self.geo.species),369)
        self.assertEqual(len(self.geo.train_sid&self.geo.test_sid),0)

    def test_schedule_roles_masks_unique_photos(self):
        geo=self.geo;idx=m.schedule(geo,'unit_test',0)
        for f in range(4):
            tr=set(geo.sid[idx[f,:20,0]]);te=set(geo.sid[idx[f,20:,0]])
            self.assertTrue(tr<=geo.train_sid);self.assertTrue(te<=geo.test_sid)
            self.assertFalse(tr&te);self.assertTrue(geo.mask[idx[f]].all())
            self.assertTrue((geo.sector[idx[f,20:]]==f).all())
            for j in range(40):
                self.assertEqual(len(set(idx[f,j])),20)
                self.assertEqual(len(set(geo.sid[idx[f,j]])),1)
                self.assertTrue(set(idx[f,j])<=set(geo.pools[(1,f,int(geo.sid[idx[f,j,0]]))]))

    def test_schedules_reproducible_and_stage_disjoint_rng(self):
        a=m.schedule(self.geo,'unit_test',0)
        np.testing.assert_array_equal(a,m.schedule(self.geo,'unit_test',0))
        self.assertFalse(np.array_equal(a,m.schedule(self.geo,'unit_test',1)))

    def test_no_truth_parameters_in_label_return(self):
        y=m.labels_for(self.geo,'unit_test','environmental_shared',1.,1.,0)
        self.assertIsInstance(y,np.ndarray);self.assertEqual(y.dtype,np.bool_)
        self.assertEqual(y.shape,(36900,))

if __name__=='__main__':unittest.main()
