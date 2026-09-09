"""Independent algebra, optimization, and denominator checks; no empirical colour."""
import sys
import unittest
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import multivariate_normal, chi2
from sklearn.metrics import adjusted_rand_score
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/analysis'))
import run_hypervolume_support_response_diagnostic as m

class DiagnosticTests(unittest.TestCase):
    def setUp(self):
        self.rng=np.random.default_rng(12345)
        self.x=self.rng.normal(size=(2,40,20,2))
        self.y=self.rng.random((2,40,20))<.5
        self.mu,self.cov,self.valid=m.base.fit_classes(self.x,self.y)

    def test_exact_gaussian_log_density(self):
        odds=m.gaussian_log_odds(self.x[:,20:],self.mu[:,:20],self.cov[:,:20])
        for b,i,j in [(0,0,0),(0,12,8),(1,3,19),(1,19,17)]:
            direct=multivariate_normal.logpdf(self.x[b,j+20],self.mu[b,i,1],self.cov[b,i,1])-multivariate_normal.logpdf(self.x[b,j+20],self.mu[b,i,0],self.cov[b,i,0])
            np.testing.assert_allclose(odds[b,i,j],direct,rtol=1e-12,atol=1e-12)

    def test_parent_predictions_unchanged(self):
        odds=m.gaussian_log_odds(self.x[:,20:],self.mu[:,:20],self.cov[:,:20])
        np.testing.assert_array_equal(odds>0,m.parent.predict_qda(self.x[:,20:],self.mu[:,:20],self.cov[:,:20]))

    def test_prior_bayes_identity(self):
        n=self.y[:,:20].sum(-1); prior=np.log(n/(20-n))
        density_odds=m.gaussian_log_odds(self.x[:,20:],self.mu[:,:20],self.cov[:,:20])
        posterior=expit(density_odds+prior[...,None,None])
        f1=multivariate_normal.pdf(self.x[0,20],self.mu[0,0,1],self.cov[0,0,1])
        f0=multivariate_normal.pdf(self.x[0,20],self.mu[0,0,0],self.cov[0,0,0])
        p=n[0,0]/20
        np.testing.assert_allclose(posterior[0,0,0],p*f1/(p*f1+(1-p)*f0),atol=1e-12)

    def test_overlap_independent_direct(self):
        w=m.overlap_weights(self.x)
        for b,i,j in [(0,0,3),(1,12,8)]:
            x=self.x[b,i]; mu=x.mean(0); cov=m.base.shrink_covariance(np.cov(x,rowvar=False)); delta=self.x[b,j+20]-mu
            q=np.sum(delta*(delta@np.linalg.inv(cov)),axis=1)
            self.assertEqual(w[b,i,j],np.mean(q<=chi2.ppf(.95,2)))
        self.assertTrue(((w>=0)&(w<=1)).all())

    def test_overlap_translation_invariant(self):
        np.testing.assert_array_equal(m.overlap_weights(self.x),m.overlap_weights(self.x+np.array([3.,-7.])))

    def test_zero_support_retains_targets(self):
        ari=np.ones((1,2,2));valid=np.ones_like(ari,dtype=bool);w=np.array([[[1.,0.],[1.,0.]]])
        self.assertEqual(float(m.aggregate(ari,valid,w)[0]),.5)

    def test_invalid_not_removed_from_sources(self):
        ari=np.ones((1,2,2));valid=np.array([[[True,True],[False,False]]]);w=np.ones_like(ari)
        self.assertEqual(float(m.aggregate(ari,valid,w)[0]),.5)

    def test_equal_weights_exact_mean(self):
        ari=self.rng.normal(size=(3,20,20));valid=self.rng.random((3,20,20))<.8
        np.testing.assert_allclose(m.aggregate(ari,valid,np.ones_like(ari)),np.where(valid,ari,0).mean((1,2)),atol=1e-15)

    def test_logit_independent_optimizer(self):
        beta,grad=m.logistic_fit(self.x[:,:20],self.y[:,:20]);self.assertLess(grad,1e-6)
        for b,i in [(0,0),(0,9),(1,3)]:
            z=np.c_[np.ones(20),self.x[b,i]]; y=self.y[b,i];pen=np.array([1e-6,1.,1.])
            def f(v):
                a=z@v
                return np.logaddexp(0,a).sum()-y@a+.5*np.sum(pen*v*v)
            def jac(v):return z.T@(expit(z@v)-y)+pen*v
            opt=minimize(f,np.zeros(3),jac=jac,method='BFGS',options={'gtol':1e-10})
            self.assertAlmostEqual(f(beta[b,i]),f(opt.x),places=10)
            np.testing.assert_allclose(beta[b,i],opt.x,atol=2e-7)

    def test_logit_colour_reversal(self):
        a,_=m.logistic_fit(self.x[:,:20],self.y[:,:20]);b,_=m.logistic_fit(self.x[:,:20],~self.y[:,:20])
        np.testing.assert_allclose(a,-b,atol=1e-10)

    def test_logit_constant_labels_numeric_guard(self):
        for y in [np.zeros((2,20),dtype=bool),np.ones((2,20),dtype=bool)]:
            beta,g=m.logistic_fit(self.x[0,:2],y);self.assertLess(g,1e-6);self.assertTrue(np.isfinite(beta).all())

    def test_seed_new_and_stage_separated(self):
        self.assertEqual(m.seed_for('test',1),m.seed_for('test',1))
        self.assertNotEqual(m.seed_for('labels','calibration',1),m.seed_for('labels','evaluation',1))
        self.assertNotEqual(m.seed_for('labels',1),m.parent.seed_for('labels',1))

    def test_ari_sklearn_nonconstant(self):
        pred=self.rng.random((2,20,20,20))<.5
        ari=m.base.adjusted_rand_binary(self.y[:,None,20:],pred)
        for b,i,j in [(0,0,0),(1,19,19)]:
            self.assertAlmostEqual(ari[b,i,j],adjusted_rand_score(self.y[b,j+20],pred[b,i,j]),places=12)

if __name__=='__main__':unittest.main()
