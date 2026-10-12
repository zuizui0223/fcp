import importlib.util
import math
from pathlib import Path
import sys
import unittest
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts' / 'analysis'))

class EntropyContract(unittest.TestCase):
    def setUp(self):
        spec = importlib.util.find_spec('fcp_hesperis_entropy_20261012')
        self.assertIsNotNone(spec, 'summary entropy decomposition module missing')
        self.m = __import__('fcp_hesperis_entropy_20261012')

    def test_uniform_and_zero_entries(self):
        self.assertAlmostEqual(self.m.entropy([1/3]*3), math.log(3))
        self.assertEqual(self.m.entropy([1, 0, 0]), 0)

    def test_reject_bad_probability(self):
        for p in ([.3, .3], [-.1, 1.1], [float('nan'), 1]):
            with self.assertRaises(ValueError): self.m.entropy(p)

    def test_identical_stands_have_no_between_component(self):
        p = [.2, .3, .5]
        d = self.m.decompose(p, self.m.entropy(p))
        self.assertAlmostEqual(d['H_between'], 0)

    def test_fixed_different_stands_have_all_between_component(self):
        d = self.m.decompose([.5, .5], 0)
        self.assertAlmostEqual(d['between_fraction'], 1)

    def test_reject_inconsistent_or_mixed_cohort_summaries(self):
        with self.assertRaises(ValueError): self.m.decompose([.5, .5], 1)
        with self.assertRaises(ValueError):
            self.m.summary_result({'year':2021, 'n_stands':132, 'alpha_n_stands':105,
                'white':.129, 'pink':.333, 'purple':.538, 'mean_shannon':.71})

    def test_exact_identity_from_observed_stand_vectors(self):
        vectors = [[.9,.05,.05], [.1,.4,.5], [.3,.3,.4]]
        mean = [sum(v[k] for v in vectors)/3 for k in range(3)]
        alpha = sum(self.m.entropy(v) for v in vectors)/3
        d = self.m.decompose(mean, alpha)
        kl = sum(sum(p*math.log(p/q) for p,q in zip(v,mean) if p) for v in vectors)/3
        self.assertAlmostEqual(d['H_between'], kl)

    def test_rounding_extrema_cover_simplex_grid(self):
        lo, hi = self.m.entropy_bounds([.15,.30,.55], .02)
        vals=[]
        for i in range(21):
            for j in range(21):
                p1=.13+.002*i; p2=.28+.002*j; p3=1-p1-p2
                if .53 <= p3 <= .57: vals.append(self.m.entropy([p1,p2,p3]))
        self.assertLessEqual(lo, min(vals)+1e-12)
        self.assertGreaterEqual(hi, max(vals)-1e-12)

    def test_uniform_entropy_maximum_is_interior(self):
        lo,hi = self.m.entropy_bounds([1/3]*3,.1)
        self.assertAlmostEqual(hi,math.log(3))
        self.assertLess(lo,hi)

    def test_infeasible_rounding_box_is_rejected(self):
        with self.assertRaises(ValueError): self.m.entropy_bounds([.1,.1,.1],.01)

    def test_single_colour_has_undefined_fraction(self):
        d=self.m.decompose([1,0],0)
        self.assertIsNone(d['between_fraction'])

    def test_source_n_not_independent_study_replications(self):
        d=self.m.summary_result({'year':2021, 'n_stands':132, 'alpha_n_stands':132,
            'white':.129, 'pink':.333, 'purple':.538, 'mean_shannon':.71})
        self.assertEqual(d['n_stands'],132)
        self.assertEqual(d['inference'],'descriptive_summary_identity_no_p_value')
        self.assertNotIn('p_value',d)
        self.assertLess(d['rounding_bounds']['between_fraction'][0],d['between_fraction'])
        self.assertGreater(d['rounding_bounds']['between_fraction'][1],d['between_fraction'])

if __name__=='__main__': unittest.main()

class RecoveryContract(unittest.TestCase):
    def setUp(self):
        spec=importlib.util.find_spec('recover_hesperis_public_data_20261012')
        self.assertIsNotNone(spec, 'bounded public-data retrieval module missing')
        self.m=__import__('recover_hesperis_public_data_20261012')

    def test_html_challenge_is_not_csv(self):
        with self.assertRaises(ValueError):
            self.m.profile_bytes('data.csv',b'<html><title>Validating...</title></html>')

    def test_schema_and_rows_are_observed_not_inferred(self):
        out=self.m.profile_bytes('data.csv',b'stand,year,white,pink,purple\nA,2021,3,2,1\n')
        self.assertEqual(out['rows'],1)
        self.assertEqual(out['headers'],['stand','year','white','pink','purple'])
        self.assertFalse(out['census_semantics_verified'])

    def test_ragged_csv_is_not_recovered_data(self):
        with self.assertRaises(ValueError):
            self.m.profile_bytes('data.csv',b'x,y\n1,2,3\n')
