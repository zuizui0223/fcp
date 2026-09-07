from pathlib import Path

import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "analysis"))

import run_global_rgfca_within_species_spatial_omnibus as run
from fcp_pipeline.global_g3 import species_distance_colour_rho


MEASURED = Path("data/derived/global_monte_carlo_measured_photos_v1.csv")


def random_simplex(rng, n):
    x = rng.gamma(shape=1.5, scale=1.0, size=(n, 4))
    return x / x.sum(axis=1, keepdims=True)


def test_frozen_contract_and_firewall_constants():
    assert run.SPECIFICATION_COMMIT == "92db63bae67a91ddb2526eb5da5922ddf2c1e7df"
    assert run.PROTOCOL == "global-rgfca-within-species-spatial-omnibus-v1"
    assert run.N_PERM == 999
    assert run.EXPECTED_ROWS == 21424
    assert run.EXPECTED_SPECIES == 369


def test_rgfca_pool_is_exact_21424_rows_369_species_and_no_small_species():
    pool = run.load_rgfca_pool(MEASURED)
    assert len(pool) == 21424
    assert pool.species.nunique() == 369
    assert int(pool.groupby("species").size().min()) >= 40
    assert not pool.photo_id.duplicated().any()


def test_fast_observed_rho_matches_frozen_g3_spearman_exactly():
    rng = np.random.default_rng(10)
    n = 17
    lat = rng.uniform(-60, 60, n)
    lon = rng.uniform(-170, 170, n)
    colours = random_simplex(rng, n)
    obs, _, _ = run.species_observed_and_null(lat, lon, colours, species="synthetic", batch_size=20)
    direct = species_distance_colour_rho(lat, lon, colours)
    assert abs(obs - direct) < 2e-12


def test_fast_permutation_rho_matches_direct_scipy_for_fixed_indices():
    rng = np.random.default_rng(11)
    n = 19
    lat = rng.uniform(-55, 55, n)
    lon = rng.uniform(-175, 175, n)
    colours = random_simplex(rng, n)
    _, null, _ = run.species_observed_and_null(lat, lon, colours, species="synthetic-perm", batch_size=31)
    for p in (0, 1, 123, 998):
        perm = np.random.default_rng(run.seed_for("synthetic-perm", p)).permutation(n)
        direct = species_distance_colour_rho(lat, lon, colours[perm])
        assert abs(float(null[p]) - direct) < 2e-12


def test_null_is_deterministic_across_batch_sizes():
    rng = np.random.default_rng(12)
    n = 14
    lat = rng.uniform(-50, 50, n)
    lon = rng.uniform(-160, 160, n)
    colours = random_simplex(rng, n)
    a, na, _ = run.species_observed_and_null(lat, lon, colours, species="deterministic", batch_size=17)
    b, nb, _ = run.species_observed_and_null(lat, lon, colours, species="deterministic", batch_size=64)
    assert a == b
    assert np.array_equal(na, nb)


def test_constant_colour_distance_species_is_defined_as_zero_like_g3():
    lat = np.array([-20.0, -5.0, 10.0, 30.0, 45.0])
    lon = np.array([-150.0, -80.0, 0.0, 70.0, 140.0])
    colours = np.tile(np.array([0.1, 0.2, 0.3, 0.4]), (5, 1))
    obs, null, _ = run.species_observed_and_null(lat, lon, colours, species="constant", batch_size=50)
    assert obs == 0.0
    assert np.array_equal(null, np.zeros(999))
    assert species_distance_colour_rho(lat, lon, colours) == 0.0


def test_species_seed_stream_differs_by_species_and_permutation():
    assert run.seed_for("A", 0) != run.seed_for("A", 1)
    assert run.seed_for("A", 0) != run.seed_for("B", 0)
    assert run.seed_for("A", 0) == run.seed_for("A", 0)


def test_runner_source_contains_no_six_or_34_species_input_route():
    source = Path(run.__file__).read_text().lower()
    assert "jbi_chapter1" not in source
    assert "34species" not in source
    assert "six_species" not in source
