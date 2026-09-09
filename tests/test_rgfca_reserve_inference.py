import numpy as np
import pandas as pd
import pytest
from scipy.stats import spearmanr

from fcp_pipeline.rgfca_reserve_inference import (
    METRICS, controlled_statistic, evaluate_species, pairwise_jsd_matrix,
    permutation_seed, summarize, vertex_permutation,
)
from fcp_pipeline.global_g3 import _upper_triangle_distances_km
from fcp_pipeline.photo_first_measurement import REFERENCE_RGB


def test_masked_null_recomputes_ranks_and_matches_scipy():
    rng = np.random.default_rng(23)
    n = 16
    geo = _upper_triangle_distances_km(rng.uniform(-50, 50, n), rng.uniform(-170, 170, n))
    matrix = pairwise_jsd_matrix(rng.dirichlet(np.ones(4), n))
    u, v = np.triu_indices(n, k=1)
    mask = (u + v) % 3 != 0
    observed, null, error = controlled_statistic(geo, matrix, species="synthetic", metric="observer_pair_exclusion", pair_mask=mask, permutations=31)
    assert np.isclose(observed, spearmanr(geo[mask], matrix[u[mask], v[mask]]).statistic)
    for index in (0, 14, 30):
        p = vertex_permutation(n, species="synthetic", index=index, metric="observer_pair_exclusion")
        assert abs(null[index] - spearmanr(geo[mask], matrix[p[u[mask]], p[v[mask]]]).statistic) < 2e-12
    assert error < 2e-12


def test_quarter_permutations_never_cross_strata():
    quarters = np.array([1, 1, 1, 2, 3, 3, 4, 4])
    for index in range(20):
        p = vertex_permutation(8, species="quarter", index=index, metric="calendar_quarter_stratification", quarters=quarters)
        assert sorted(p) == list(range(8))
        assert np.array_equal(quarters[p], quarters)
        assert p[3] == 3


def test_same_photo_flower_background_pairing():
    rng = np.random.default_rng(15)
    n = 11
    flower = rng.dirichlet(np.ones(12), n)
    background = .8 * flower + .2 * rng.dirichlet(np.ones(12), n)
    diff = pairwise_jsd_matrix(flower) - pairwise_jsd_matrix(background)
    p = vertex_permutation(n, species="paired", index=7, metric="matched_background_differential")
    assert np.allclose(diff[np.ix_(p, p)], pairwise_jsd_matrix(flower[p]) - pairwise_jsd_matrix(background[p]))
    wrong = pairwise_jsd_matrix(flower[p]) - pairwise_jsd_matrix(background)
    assert not np.allclose(wrong, diff[np.ix_(p, p)])


def test_seed_rule_is_literal_contract_not_imported_background_seed():
    import hashlib
    expected = int.from_bytes(hashlib.sha256(b"202609071503|A|9").digest()[:8], "little")
    assert permutation_seed("A", 9, "matched_background_differential") == expected
    assert expected != permutation_seed("B", 9, "matched_background_differential")


def test_constant_geometry_is_not_silently_a_biological_zero():
    with pytest.raises(ValueError, match="not_evaluable_pair_geometry"):
        controlled_statistic(np.zeros(6), np.zeros((4, 4)), species="zero", metric="observer_pair_exclusion", permutations=2)


def test_full_species_execution_includes_all_fixed_controls():
    rng = np.random.default_rng(56)
    n = 40
    x = pd.DataFrame({"inat_taxon_id": [991] * n, "species": ["synthetic full"] * n, "photo_id": np.arange(n),
                      "latitude": rng.uniform(-50, 50, n), "longitude": rng.uniform(-170, 170, n),
                      "observer_id": np.repeat(np.arange(20), 2), "observed_on": [f"2020-{1+(i%4)*3:02d}-15" for i in range(n)]})
    from fcp_pipeline.rgfca_reserve_replication import COLOURS
    for index, col in enumerate(COLOURS):
        x[col] = rng.dirichlet(np.ones(4), n)[:, index]
    x[COLOURS] = x[COLOURS].div(x[COLOURS].sum(axis=1), axis=0)
    for col in REFERENCE_RGB:
        x[f"palette_count_{col}"] = rng.integers(1, 100, n)
        x[f"background_palette_count_{col}"] = rng.integers(1, 100, n)
    row, null = evaluate_species(x)
    assert null.shape == (4, 999) and np.isfinite(null).all()
    assert row["direct_checks"] == 7 and row["maximum_direct_check_error"] < 2e-12
    assert abs(row["rho_primary"] - row["rho_calendar_quarter_stratification"]) < 2e-12
    assert all(f"rho_{k}" in row for k in METRICS)


def test_conjunction_cannot_rescue_a_failed_control():
    n = 250
    rows = pd.DataFrame({f"rho_{k}": np.full(n, .1) for k in METRICS})
    rng = np.random.default_rng(86)
    null = rng.normal(0, .1, (n, 4, 999))
    null[:, -1, :] += .2
    out = summarize(rows, null)
    assert out["directional_photo_association_replicated"] is True
    assert out["flower_specific_robust_replication"] is False
    assert out["metrics"]["matched_background_differential"]["p_upper"] == 1.0


def test_incomplete_null_species_census_rejected():
    with pytest.raises(ValueError, match="incomplete"):
        summarize(pd.DataFrame({"rho_primary": [0] * 249}), np.zeros((249, 4, 999)))


def test_finalizer_checks_census_hashes_and_records_failure(monkeypatch, tmp_path):
    import json
    from scripts.analysis import run_rgfca_reserve_inference as run
    from fcp_pipeline.rgfca_reserve_replication import sha
    pool = pd.DataFrame({"inat_taxon_id": np.repeat(np.arange(250), 40)})
    monkeypatch.setattr(run, "load_pool", lambda: (pool, {"github_run_id": "synthetic-measurement"}))
    monkeypatch.setattr(run, "lineage", lambda: {"synthetic": True})
    root = tmp_path / "shards"
    root.mkdir()
    for index in range(20):
        directory = root / f"shard-{index}"
        directory.mkdir()
        taxa = list(range(250))[index::20]
        frame = pd.DataFrame({"inat_taxon_id": taxa, "photos": 40, "maximum_direct_check_error": 0., "direct_checks": 7,
                              **{f"rho_{k}": .01 for k in METRICS}})
        frame.to_csv(directory / "species.csv", index=False)
        np.savez_compressed(directory / "null.npz", taxon_ids=taxa, null=np.zeros((len(taxa), 4, 999)))
        receipt = {"status": "complete_reserve_species_shard", "shard_index": index,
                   "lineage": {"synthetic": True}, "metrics": list(METRICS), "permutations": 999,
                   "eligible_species": 250, "taxa": taxa,
                   "files_sha256": {f: sha((directory / f).read_bytes()) for f in ("species.csv", "null.npz")}}
        (directory / "receipt.json").write_text(json.dumps(receipt))
    run.finalize(root, tmp_path / "valid")
    result = json.loads((tmp_path / "valid/result.json").read_text())
    assert result["eligible_species"] == 250 and result["eligible_photos"] == 10000
    assert result["direct_checks"] == 1750
    assert result["flower_specific_robust_replication"] is False  # degenerate null never passes
    failure = root / "shard-0/receipt.json"
    data = json.loads(failure.read_text())
    data["status"] = "not_evaluable_invalid_eligible_species"
    failure.write_text(json.dumps(data))
    run.finalize(root, tmp_path / "failed")
    result = json.loads((tmp_path / "failed/result.json").read_text())
    assert result["primary_computed"] is False and result["replacement_used"] is False
    assert not (tmp_path / "failed/species.csv").exists()
    data["status"] = "complete_reserve_species_shard"
    failure.write_text(json.dumps(data))
    (root / "shard-0/species.csv").write_text("tampered")
    with pytest.raises(ValueError, match="hash mismatch"):
        run.finalize(root, tmp_path / "tampered")
