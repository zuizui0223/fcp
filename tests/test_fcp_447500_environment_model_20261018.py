"""Mathematical and selection checks for the executed two-stage model."""
import unittest
import numpy as np
import pandas as pd
from scipy.optimize._numdiff import approx_derivative
from scripts.analysis.run_fcp_447500_environment_model_20261018 import (
    firth_objective,fit_firth,reml_meta,thin_photos,eligible)

class EnvironmentModelTests(unittest.TestCase):
    def test_firth_gradient_independent_finite_difference(self):
        rng=np.random.default_rng(1);x=np.column_stack([np.ones(60),rng.normal(size=(60,3))]);y=rng.integers(0,2,60)
        beta=np.array([-.2,.5,-.8,.2])
        numeric=approx_derivative(lambda b:firth_objective(b,x,y)[0],beta,method='3-point').ravel()
        np.testing.assert_allclose(firth_objective(beta,x,y)[1],numeric,atol=2e-7)

    def test_intercept_only_firth_has_analytic_half_count_solution(self):
        y=np.r_[np.ones(7),np.zeros(43)];fit=fit_firth(np.ones((50,1)),y)
        self.assertEqual(fit['status'],'ok')
        self.assertAlmostEqual(fit['beta'][0],np.log(7.5/43.5),places=6)
        self.assertAlmostEqual(fit['cov'][0,0],1/(50*(7.5/51)*(43.5/51)),places=6)

    def test_complete_separation_finite_and_response_complement_symmetry(self):
        z=np.linspace(-2,2,60);x=np.column_stack([np.ones(60),z]);y=(z>0).astype(float)
        a=fit_firth(x,y);b=fit_firth(x,1-y)
        self.assertEqual(a['status'],'ok');self.assertEqual(b['status'],'ok')
        self.assertTrue(np.isfinite(a['beta']).all());self.assertGreater(a['beta'][1],0)
        np.testing.assert_allclose(a['beta'],-b['beta'],atol=1e-5)
        np.testing.assert_allclose(a['cov'],b['cov'],atol=1e-5)

    def test_rank_deficiency_not_silently_dropped(self):
        z=np.arange(50);x=np.column_stack([np.ones(50),z,2*z])
        self.assertEqual(fit_firth(x,z%2)['status'],'rank_deficient')

    def test_reml_equal_variance_matches_analytic_sample_variance(self):
        y=np.array([-2.,-1.,0.,1.,2.]);v=np.repeat(.25,5)
        result,mean,var,prob=reml_meta(y,v)
        self.assertAlmostEqual(result['tau2'],y.var(ddof=1)-.25,places=6)
        self.assertAlmostEqual(result['reml_mean'],0.,places=7)
        self.assertAlmostEqual(result['equal_species_unshrunk_se_independence'],np.sqrt(v.sum())/5,places=9)
        self.assertTrue(result['tau2_profile_ci_low']<result['tau2']<result['tau2_profile_ci_high'])
        np.testing.assert_allclose(mean,.9*y,atol=1e-6)
        self.assertAlmostEqual(prob.mean(),.5,places=7)

    def test_reml_zero_boundary_mhksj_does_not_shrink_uncertainty(self):
        result,_,_,_=reml_meta(np.zeros(6),np.ones(6))
        self.assertEqual(result['tau2'],0.)
        self.assertEqual(result['hksj_q_used'],1.)
        self.assertAlmostEqual(result['reml_mean_se_mhksj'],np.sqrt(1/6))

    def test_thinning_response_blind_order_invariant_and_selection(self):
        g=pd.DataFrame({'photo_id':np.arange(50),'observer_id':np.repeat(np.arange(25),2),
                        'latitude':np.repeat(np.arange(25),2)+.1,'longitude':np.zeros(50),
                        'white':[True]*5+[False]*45})
        self.assertTrue(eligible(g));self.assertFalse(eligible(g,10))
        for mode in ('observer_thinned','spatial_thinned'):
            a=thin_photos(g,mode);b=thin_photos(g.sample(frac=1,random_state=3).assign(white=lambda d:~d.white),mode)
            self.assertEqual(set(a.photo_id),set(b.photo_id));self.assertEqual(len(a),25)
            self.assertFalse(eligible(a))

if __name__=='__main__':unittest.main()
