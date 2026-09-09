from pathlib import Path
import json
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "analysis"))

import finalize_global_rgfca_background_control as fin
import build_global_rgfca_background_control_firewall as fw


CONTRACT = Path("docs/supporting/global_rgfca_within_species_spatial_background_control_contract_v2.json")
AMENDMENT = Path("docs/supporting/global_rgfca_background_control_postoutcome_status_amendment_v2a.json")
MEASURED = Path("data/derived/global_monte_carlo_measured_photos_v1.csv")
PREFLIGHT = Path("docs/supporting/global_rgfca_background_control_preflight_v2.json")


def random_simplex(rng, n, k=12):
    x = rng.gamma(1.3, 1.0, size=(n, k))
    return x / x.sum(axis=1, keepdims=True)


def test_contract_is_postoutcome_falsification_and_sharedness_remains_mainline():
    contract = json.loads(CONTRACT.read_text())
    amendment = json.loads(AMENDMENT.read_text())
    preflight = json.loads(PREFLIGHT.read_text())
    assert amendment["status"] == "frozen_after_flower_only_omnibus_outcome_but_before_any_background_colour_recovery"
    assert amendment["known_flower_only_result_before_background_control"]["permutation_p_upper"] == 0.001
    assert amendment["still_unopened_at_this_amendment"]["recovered_background_palette_vectors"] is False
    assert contract["main_line_firewall"]["scientific_main_question"] == "cross-species sharedness of flower-colour transition structure"
    assert contract["main_line_firewall"]["this_diagnostic_is_main_claim"] is False
    assert contract["joint_exact_randomization_null"]["permutations"] == 999
    assert contract["joint_exact_randomization_null"]["master_seed"] == 202609071503
    assert preflight["status"] == "pass_background_control_preflight_without_pixels"
    assert preflight["background_colour_opened"] is False
    assert preflight["matched_rows"] == 21424
    assert preflight["matched_species"] == 369


def test_pairwise_twelve_palette_jsd_is_symmetric_zero_diagonal_and_bounded():
    rng = np.random.default_rng(44)
    p = random_simplex(rng, 17)
    d = fin.pairwise_jsd_matrix(p)
    assert d.shape == (17, 17)
    assert np.allclose(d, d.T, atol=0, rtol=0)
    assert np.max(np.abs(np.diag(d))) < 1e-14
    assert np.min(d) >= 0.0
    assert np.max(d) <= 1.0 + 1e-14


def test_fast_joint_vertex_permutation_matches_direct_scipy():
    rng = np.random.default_rng(45)
    n = 21
    lat = rng.uniform(-55, 55, n)
    lon = rng.uniform(-170, 170, n)
    flower = random_simplex(rng, n)
    background = random_simplex(rng, n)
    geo = fin.pairwise_geo_km(lat, lon)
    fj = fin.pairwise_jsd_matrix(flower)
    bj = fin.pairwise_jsd_matrix(background)
    differential = fj - bj
    obs, null = fin.observed_and_null_differential(
        geo,
        differential,
        species="synthetic-paired",
        master_seed=202609071503,
        permutations=37,
        batch_size=7,
    )
    ui, uj = np.triu_indices(n, k=1)
    direct_obs = float(spearmanr(geo[ui, uj], differential[ui, uj]).statistic)
    assert abs(obs - direct_obs) < 2e-12
    for p in (0, 1, 19, 36):
        perm = np.random.default_rng(fin.seed_for(202609071503, "synthetic-paired", p)).permutation(n)
        # The same vertex permutation is applied to the whole flower-minus-background
        # matrix; equivalently flower and matched background from each photograph move together.
        direct = float(spearmanr(geo[ui, uj], differential[np.ix_(perm, perm)][ui, uj]).statistic)
        assert abs(float(null[p]) - direct) < 2e-12


def test_joint_pairing_differs_from_illegal_independent_flower_background_permutation():
    rng = np.random.default_rng(46)
    n = 18
    flower = random_simplex(rng, n)
    # Make background correlated with its matched flower while still not identical.
    background = 0.85 * flower + 0.15 * random_simplex(rng, n)
    background /= background.sum(axis=1, keepdims=True)
    fj = fin.pairwise_jsd_matrix(flower)
    bj = fin.pairwise_jsd_matrix(background)
    paired = np.random.default_rng(100).permutation(n)
    independently_permuted_background = np.random.default_rng(101).permutation(n)
    legal = (fj - bj)[np.ix_(paired, paired)]
    illegal = fj[np.ix_(paired, paired)] - bj[np.ix_(independently_permuted_background, independently_permuted_background)]
    assert not np.allclose(legal, illegal)


def test_seed_stream_is_species_and_permutation_specific_and_deterministic():
    master = 202609071503
    assert fin.seed_for(master, "A", 0) == fin.seed_for(master, "A", 0)
    assert fin.seed_for(master, "A", 0) != fin.seed_for(master, "A", 1)
    assert fin.seed_for(master, "A", 0) != fin.seed_for(master, "B", 0)


def test_firewall_builder_source_does_not_open_flower_or_background_outcome():
    source = Path(fw.__file__).read_text().lower()
    assert "photo_url_large" in source  # only to construct sealed acquisition key
    assert "background_palette_count_" not in source
    assert "nearest_palette_counts" not in source
    assert "image.open" not in source
    assert "six_species" not in source
    assert "34species" not in source


def test_failure_path_is_fail_closed_before_primary(monkeypatch, tmp_path):
    # The finalizer checks the complete 21,424-row recovery census before joining biological data.
    # Construct that exact census with one failed row; no primary statistic is permitted.
    root = tmp_path / "recovery"
    root.mkdir()
    n = 21424
    ids = np.array([f"X{i:05d}" for i in range(n)], dtype=object)
    status = np.full(n, "exact_matched_background_recovered", dtype=object)
    status[-1] = "image_sha_mismatch_no_replacement"
    splits = np.array_split(np.arange(n), 128)
    for k, idx in enumerate(splits):
        s, p = divmod(k, 4)
        pd.DataFrame({"measurement_id": ids[idx], "recovery_status": status[idx]}).to_csv(
            root / f"recovery_s{s:02d}_p{p:02d}.csv", index=False
        )
    output = tmp_path / "final"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "finalize",
            "--recovery-dir", str(root),
            "--measured", str(MEASURED),
            "--contract", str(CONTRACT),
            "--amendment", str(AMENDMENT),
            "--output-dir", str(output),
        ],
    )
    assert fin.main() == 0
    summary = json.loads((output / "summary.json").read_text())
    assert summary["status"] == "not_evaluable_incomplete_exact_background_recovery"
    assert summary["exact_rows"] == 21423
    assert summary["failed_rows"] == 1
    assert summary["primary_statistic_computed"] is False
    assert summary["replacement_photos_used"] is False
    assert summary["denominator_adapted_after_recovery"] is False
    assert not (output / "species_results_v2.csv").exists()
    assert not (output / "global_null_v2.csv").exists()
