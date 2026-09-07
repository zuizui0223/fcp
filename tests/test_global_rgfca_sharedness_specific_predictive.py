import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "analysis"))

import run_global_rgfca_sharedness_specific_predictive as run


GEOMETRY = Path("inputs/geometry_only.csv")


def test_frozen_contract_constants_and_arm_counts():
    assert run.SPECIFICATION_COMMIT == "34fe3093cd17c6d20fb8044259984314761cd210"
    assert run.REPS == 250
    assert run.N_TRAIN == 100
    assert run.N_TEST == 100
    assert run.N_PHOTOS == 20
    assert len(run.NUISANCE) == 9
    assert len(run.POSITIVE) == 10
    assert len(run.NUISANCE) * run.REPS == 2250
    assert (len(run.NUISANCE) + len(run.POSITIVE)) * run.REPS == 4750


def test_structured_independent_boundaries_remain_in_null():
    assert all(float(s.get("shared_fraction", 0.0)) == 0.0 for s in run.NUISANCE)
    arms = {(s["world"], float(s["amplitude"])) for s in run.NUISANCE}
    assert ("independent_boundaries", 0.5) in arms
    assert ("independent_boundaries", 1.0) in arms
    assert ("independent_boundaries", 2.0) in arms
    assert ("independent_shifted_boundaries", 2.0) in arms
    assert any(s["world"] == "directionally_clustered_independent" for s in run.NUISANCE)
    assert any(s["world"] == "bimodal_clustered_independent" for s in run.NUISANCE)


def test_fixed_axis_grid_is_deterministic_unit_and_antipodally_canonical():
    a = run.fixed_axes()
    b = run.fixed_axes()
    assert a.shape == (96, 3)
    assert np.array_equal(a, b)
    assert np.allclose(np.linalg.norm(a, axis=1), 1.0, atol=1e-12)
    assert len(np.unique(np.round(a, 12), axis=0)) == 96
    for row in a:
        nz = np.flatnonzero(~np.isclose(row, 0.0))
        assert len(nz) > 0
        assert row[nz[0]] > 0


def test_geometry_recovers_rgfca_only_369_species_frame_and_fixed_split():
    geo = run.Geometry(GEOMETRY)
    assert len(geo.species) == 369
    assert int(geo.mask.sum()) == 21424
    assert len(geo.train_ids) == 184
    assert len(geo.test_ids) == 185
    assert not (set(geo.train_ids) & set(geo.test_ids))
    assert min(len(pool) for pool in geo.pools.values()) >= 40


def test_schedule_is_deterministic_species_disjoint_and_without_photo_replacement():
    geo = run.Geometry(GEOMETRY)
    a = geo.schedule("calibration", 7)
    b = geo.schedule("calibration", 7)
    assert np.array_equal(a, b)
    assert a.shape == (200, 20)
    sid = geo.sid[a[:, 0]]
    assert len(np.unique(sid)) == 200
    assert set(sid[:100]).issubset(geo.train_sid)
    assert set(sid[100:]).issubset(geo.test_sid)
    assert not (set(sid[:100]) & set(sid[100:]))
    for row in a:
        assert len(np.unique(row)) == 20


def test_shared_generator_places_exact_common_boundary_in_prespecified_fraction():
    geo = run.Geometry(GEOMETRY)
    scenario = dict(
        world="partially_shared_boundaries",
        amplitude=1.0,
        shared_fraction=0.5,
        threshold_sd=0.25,
    )
    _, normals, threshold, sign = run.generator_params(geo, "evaluation", scenario, 3)
    assert normals.shape == (369, 3)
    assert threshold.shape == (369,)
    assert sign.shape == (369,)
    rounded = np.round(normals, 12)
    _, counts = np.unique(rounded, axis=0, return_counts=True)
    assert counts.max() == round(369 * 0.5)


def test_labels_are_deterministic_and_independent_null_is_spatially_structured():
    geo = run.Geometry(GEOMETRY)
    scenario = dict(world="independent_boundaries", amplitude=1.0, shared_fraction=0.0, threshold_sd=0.0)
    y1 = run.labels_for(geo, "calibration", scenario, 11)
    y2 = run.labels_for(geo, "calibration", scenario, 11)
    assert np.array_equal(y1, y2)
    assert y1.dtype == bool
    assert y1.shape == (len(geo.df),)
    assert 0 < int(y1.sum()) < len(y1)


def test_axis_log_likelihood_matches_direct_one_species_one_axis_recomputation():
    rng = np.random.default_rng(1234)
    x = rng.normal(size=(200, 20, 3))
    x /= np.linalg.norm(x, axis=2, keepdims=True)
    y = rng.integers(0, 2, size=(1, 200, 20), dtype=np.int8)
    features = run.precompute_features(x)
    ll = run.axis_log_likelihood(y, features)

    species = 0
    axis = 0
    projection = x[species] @ run.AXES[axis]
    terms = []
    for slope in run.SLOPES:
        for offset in run.OFFSETS:
            z = slope * (projection - offset)
            lp = -np.logaddexp(0.0, -z)
            lq = -np.logaddexp(0.0, z)
            yy = y[0, species].astype(float)
            pos = float(np.sum(yy * lp + (1.0 - yy) * lq))
            neg = float(np.sum(yy * lq + (1.0 - yy) * lp))
            terms.append(np.logaddexp(pos, neg) - np.log(2.0))
    direct = float(logsumexp(terms) - np.log(len(run.SLOPES) * len(run.OFFSETS)))
    assert math.isclose(ll[0, species, axis], direct, rel_tol=0.0, abs_tol=1e-10)


def test_evaluation_likelihood_never_updates_training_posterior():
    rng = np.random.default_rng(55)
    ll = rng.normal(size=(2, 200, run.K))
    ll[1, :100] = ll[0, :100]
    ll[1, 100:] = rng.normal(loc=3.0, scale=2.0, size=(100, run.K))
    score, mean_f, entropy, max_mass = run.score_from_ll(ll)
    assert score.shape == (2,)
    assert np.isclose(mean_f[0], mean_f[1], atol=1e-12)
    assert np.isclose(entropy[0], entropy[1], atol=1e-12)
    assert np.isclose(max_mass[0], max_mass[1], atol=1e-12)


def test_uninformative_axis_likelihood_has_zero_predictive_gain():
    ll = np.zeros((1, 200, run.K), dtype=float)
    score, mean_f, entropy, max_mass = run.score_from_ll(ll)
    assert abs(float(score[0])) < 1e-12
    assert 0.0 <= float(mean_f[0]) <= 1.0
    assert np.isclose(float(entropy[0]), 1.0, atol=1e-12)
    assert np.isclose(float(max_mass[0]), 1.0 / run.K, atol=1e-12)


def test_positive_grid_contains_partial_and_full_sharing_without_posthoc_gap():
    keys = {(float(s["amplitude"]), float(s["shared_fraction"])) for s in run.POSITIVE if s["world"] == "partially_shared_boundaries"}
    for pair in [(1.0, 0.1), (1.0, 0.25), (1.0, 0.5), (1.0, 0.75), (1.0, 1.0), (2.0, 0.25), (2.0, 0.5), (2.0, 1.0)]:
        assert pair in keys


def test_no_six_species_or_34species_inputs_are_imported_by_runner_contract():
    source = Path(run.__file__).read_text()
    lowered = source.lower()
    assert "jbi_chapter1" not in lowered
    assert "34species" not in lowered
    assert "six_species" not in lowered
