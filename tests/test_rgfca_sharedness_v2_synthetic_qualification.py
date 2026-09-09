from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "scripts/analysis"
sys.path.insert(0, str(ANALYSIS))
spec = importlib.util.spec_from_file_location("v2", ANALYSIS / "run_rgfca_sharedness_v2_synthetic_qualification.py")
v2 = importlib.util.module_from_spec(spec); assert spec.loader; spec.loader.exec_module(v2)
bench_spec = importlib.util.spec_from_file_location("bench", ANALYSIS / "benchmark_rgfca_sharedness_v2_exact_kernel.py")
bench = importlib.util.module_from_spec(bench_spec); assert bench_spec.loader; bench_spec.loader.exec_module(bench)


def test_frozen_geometry_census_and_split():
    _, mapping, frame, xyz, train, test = v2.load_inputs()
    assert len(frame) == 45000 and xyz.shape == (150, 300, 3)
    assert len(train) == 75 and len(test) == 75 and not (set(train) & set(test))
    assert mapping["execution"]["observed_flower_colour_opened"] is False


def test_arm_census_is_frozen():
    _, mapping, *_ = v2.load_inputs()
    assert len(v2.nuisance_arms(mapping)) == 21
    assert len(v2.positive_scenarios(mapping)) == 19
    assert len(v2.evaluation_arms(mapping)) == 78


def test_positive_grid_contains_both_threshold_sd_values_for_hard_cells():
    _, mapping, *_ = v2.load_inputs()
    p = v2.positive_scenarios(mapping)
    for a, f in [(1.0,1.0),(1.0,0.5),(2.0,0.5)]:
        vals = sorted(x["threshold_sd"] for x in p if x["amplitude"]==a and x["shared_fraction"]==f)
        assert vals == [0.0, 0.25]


def test_exact_kernel_matches_reference():
    rng=np.random.default_rng(8); n=14
    xyz=bench.unit_vectors(rng,(n,3)); flower=rng.normal(size=(n,3)); background=rng.normal(size=(n,3))
    a=bench.exact_axis_scores_vectorized(xyz,flower,background); b=bench.exact_axis_scores_reference(xyz,flower,background)
    assert np.array_equal(np.isfinite(a),np.isfinite(b))
    k=np.isfinite(a); assert np.allclose(a[k],b[k],atol=1e-11,rtol=1e-11)


def test_zero_log_evidence_predictive_score_is_zero():
    _,_,_,_,train,test=v2.load_inputs(); lr=np.zeros((150,96))
    score, mean_f, entropy, max_mass=v2.predictive_score(lr,train,test)
    assert abs(score) < 1e-12
    assert np.isfinite([mean_f,entropy,max_mass]).all()


def test_evaluation_does_not_change_training_posterior_summary():
    _,_,_,_,train,test=v2.load_inputs(); rng=np.random.default_rng(9)
    lr=rng.normal(size=(150,96)); a=v2.predictive_score(lr,train,test)
    lr2=lr.copy(); lr2[test]+=rng.normal(scale=4,size=(len(test),96)); b=v2.predictive_score(lr2,train,test)
    # posterior summaries depend on training only; score itself may change.
    assert np.allclose(a[1:],b[1:],atol=1e-12,rtol=1e-12)


def test_retention_exact_counts_and_scenario_independent_default():
    _,mapping,_,xyz,_,_=v2.load_inputs(); root=int(mapping["execution"]["deterministic_seed_root"])
    for r,k in [(0.4,120),(0.6,180),(0.8,240)]:
        a=v2.retention_indices(root,"calibration",3,r,xyz,False)
        b=v2.retention_indices(root,"calibration",3,r,xyz,False)
        assert a.shape==(150,k) and np.array_equal(a,b)
        assert all(len(np.unique(row))==k for row in a)


def test_calibration_arms_are_all_structured_nulls():
    _,mapping,*_=v2.load_inputs()
    for a in v2.nuisance_arms(mapping):
        assert a["kind"]=="nuisance"
        assert float(a.get("shared_fraction",0.0))==0.0


def test_seed_namespaces_differ_calibration_evaluation():
    _,mapping,*_=v2.load_inputs(); root=int(mapping["execution"]["deterministic_seed_root"])
    assert v2.seed_for(root,"world","calibration","x",0) != v2.seed_for(root,"world","evaluation","x",0)


def test_runner_contains_no_old_6_or_34_inputs():
    text=(ANALYSIS/"run_rgfca_sharedness_v2_synthetic_qualification.py").read_text()
    forbidden=["jbi_ch1_stage", "34species", "six_species", "data/derived/global_monte_carlo_measured_photos_v1.csv"]
    assert not any(x in text for x in forbidden)
