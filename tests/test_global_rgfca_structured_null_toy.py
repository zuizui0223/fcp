"""Deterministic checks of toy identities, marginals, validation and generators."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np
P = Path(__file__).resolve().parents[1]/'scripts/analysis/run_global_rgfca_structured_null_toy.py'
spec = importlib.util.spec_from_file_location('toy', P)
toy = importlib.util.module_from_spec(spec); spec.loader.exec_module(toy)

class ToyTests(unittest.TestCase):
    def test_shape(self):
        x = toy.shape()
        self.assertAlmostEqual(x.mean(), 0, places=13)
        self.assertAlmostEqual(x.std(), 1, places=13)
        self.assertAlmostEqual(np.mean([np.roll(x,i) for i in range(toy.CELLS)],axis=0).max(),0,places=13)
    def test_covariance_identity(self):
        x = np.random.default_rng(4).normal(size=(7,5,13)); x-=x.mean(axis=2,keepdims=True)
        values = toy.decomposition(x)
        direct = [2*np.triu(v@v.T/13,k=1).sum()/25 for v in x]
        np.testing.assert_allclose(values[:,2], direct, rtol=1e-12,atol=1e-14)
        np.testing.assert_allclose(values[:,0], values[:,1]+values[:,2], atol=1e-14)
    def test_identical_species(self):
        x = np.tile(toy.shape(),(1,toy.SPECIES,1))
        np.testing.assert_allclose(toy.decomposition(x)[0], [1,1/40,39/40], atol=1e-13)
    def test_centering(self):
        x = np.random.default_rng(5).normal(size=(2,4,10))
        np.testing.assert_allclose(toy.decomposition(x),toy.decomposition(x+9),atol=1e-14)
    def test_scramble_preserves_marginals(self):
        x=toy.generate(np.random.default_rng(2),4,1.0,10)
        y=np.random.default_rng(3).permuted(x,axis=2)
        np.testing.assert_array_equal(np.sort(x,axis=2),np.sort(y,axis=2))
        np.testing.assert_allclose(toy.decomposition(x)[:,1],toy.decomposition(y)[:,1],atol=1e-14)
    def test_reproducible(self):
        np.testing.assert_array_equal(toy.generate(toy.stream(1,2),3,0.5,4),toy.generate(toy.stream(1,2),3,0.5,4))
    def test_expectation_no_sharing(self):
        self.assertEqual(toy.expected(2,0)[2],0)
        self.assertEqual(toy.expected(2,1)[2],0)
        self.assertAlmostEqual(toy.expected(1,10)[2],90/1600)
    def test_input_validation(self):
        for x in [np.ones((2,2)),np.full((2,2,3),np.nan),np.ones((2,1,3))]:
            with self.assertRaises(ValueError): toy.decomposition(x)
        with self.assertRaises(ValueError): toy.wilson(2,1)
    def test_wilson_boundaries(self):
        self.assertLessEqual(toy.wilson(0,100)[0],1e-15)
        self.assertGreaterEqual(toy.wilson(100,100)[1],1-1e-15)
        self.assertLess(toy.wilson(50,100)[0],.5)
        self.assertGreater(toy.wilson(50,100)[1],.5)

if __name__=='__main__': unittest.main()
