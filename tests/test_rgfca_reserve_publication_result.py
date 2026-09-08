"""Reconstruct already-completed reserve evidence; never repeat image inference.

Exact Git bytes and stored global nulls are independent of workspace newlines.
Only the originally fixed bootstrap is replayed, in its original shard order.
"""
from functools import lru_cache
import hashlib
import io
import json
from pathlib import Path
import subprocess

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
REV = "ef00a78a2eda7bc79f7e6b88f719a8a696dc8b43"
AUTH = "391caecaa1a6d3c3e2407f5f19d0bad7057922e0"
RESULT_PATH = "docs/supporting/rgfca_reserve_replication_result_v1.json"
RESULT_SHA = "1cb9de6feb8bda3bf05f3fbe51abd8eeebb41774a9cd8253034d618b78dcb3f3"
METRICS = ("primary", "observer_pair_exclusion", "calendar_quarter_stratification",
           "matched_background_differential")
EXPECTED = ((.025482606069841617, .001, True), (.02521802297988336, .001, True),
            (.025482606069841617, .001, True), (.004477251893133873, .087, False))


def git_bytes(path, revision=REV):
    return subprocess.check_output(["git", "show", f"{revision}:{path}"], cwd=ROOT)


@lru_cache(None)
def sources():
    raw = git_bytes(RESULT_PATH)
    assert hashlib.sha256(raw).hexdigest() == RESULT_SHA
    assert git_bytes(RESULT_PATH, "HEAD") == raw
    result = json.loads(raw)
    tables = {}
    for name, suffix in (("species.csv", "species"), ("global_null.csv", "global_null")):
        path = f"data/derived/rgfca_reserve_replication_{suffix}_v1.csv"
        raw = git_bytes(path)
        assert hashlib.sha256(raw).hexdigest() == result["files_sha256"][name]
        assert git_bytes(path, "HEAD") == raw
        tables[suffix] = pd.read_csv(io.BytesIO(raw))
    return result, tables["species"], tables["global_null"]


def test_completed_result_and_full_species_photo_shard_census():
    result, species, null = sources()
    assert result["status"] == "complete_reserve_replication_with_fixed_controls"
    assert result["github_run_id"] == "34178957447"
    assert result["measurement_run_id"] == "34091091640"
    assert result["raw_species"] == 500 and result["raw_photos"] == 50000
    assert result["eligible_species"] == len(species) == species.inat_taxon_id.nunique() == 363
    assert result["eligible_photos"] == species.photos.sum() == 20903
    assert species.photos.min() >= 40 and result["primary_computed"] is True
    assert len(result["source_shard_receipts"]) == 20
    taxa = sorted(species.inat_taxon_id.tolist())
    for index, receipt in enumerate(result["source_shard_receipts"]):
        assert receipt["shard_index"] == index and receipt["shards"] == 20
        assert receipt["taxa"] == taxa[index::20]
        assert receipt["metrics"] == list(METRICS) and receipt["permutations"] == 999
        assert receipt["lineage"] == result["lineage"]
        assert receipt["github_run_id"] == result["github_run_id"]
        assert receipt["status"] == "complete_reserve_species_shard"
        assert set(receipt["files_sha256"]) == {"species.csv", "null.npz"}
    assert null.permutation_index.tolist() == list(range(999))
    assert np.isfinite(null[list(METRICS)].to_numpy()).all()


def test_authorized_inputs_and_complete_measurement_audit_unchanged():
    result, _, _ = sources()
    auth_path = "docs/supporting/rgfca_reserve_inference_authorization_v1.json"
    auth = json.loads(git_bytes(auth_path, AUTH))
    assert git_bytes(auth_path) == git_bytes(auth_path, AUTH)
    assert len(auth["files_sha256"]) == 16
    for path, expected in auth["files_sha256"].items():
        raw = git_bytes(path)
        assert hashlib.sha256(raw).hexdigest() == expected, path
        assert raw == git_bytes(path, AUTH) == git_bytes(path, "HEAD"), path
    for path, expected in result["lineage"]["code_sha256"].items():
        assert expected == auth["files_sha256"][path]
    audit = json.loads(git_bytes("docs/supporting/rgfca_reserve_complete_census_audit_v1.json", AUTH))
    assert audit["terminal_partitions_verified"] == 256 and audit["terminal_files_verified"] == 512
    assert audit["eligible_species"] == 363 and audit["eligible_photos"] == 20903
    assert audit["classifiable_photos"] == 24885 and audit["raw_photos"] == 50000
    assert sum(audit["terminal_status_counts"].values()) == 50000
    assert audit["prior_photo_observation_overlap"] == 0
    assert audit["measurement_gate_pass"] is True
    assert audit["inference_run_by_audit"] is False


@pytest.mark.parametrize("metric,expected", zip(METRICS, EXPECTED))
def test_every_observed_mean_and_stored_randomization_tail(metric, expected):
    result, species, null = sources()
    observed = species[f"rho_{metric}"].mean()
    draws = null[metric].to_numpy()
    p = (1 + np.count_nonzero(draws >= observed)) / 1000
    fixed = result["metrics"][metric]
    assert observed == pytest.approx(expected[0], rel=0, abs=2e-12)
    assert observed == pytest.approx(fixed["mean_rho"], rel=0, abs=2e-12)
    assert p == fixed["p_upper"] == expected[1]
    assert draws.mean() == pytest.approx(fixed["null_mean"], rel=0, abs=2e-12)
    assert draws.std(ddof=1) == pytest.approx(fixed["null_sd"], rel=0, abs=2e-12)
    assert np.ptp(draws) > 1e-15 and fixed["null_nondegenerate"] is True
    assert fixed["status"] == "evaluable"
    assert fixed["positive_supported"] is expected[2]
    assert bool(observed > 0 and p < .05) is expected[2]


def test_fixed_bootstrap_reconstructed_in_original_shard_order():
    result, species, _ = sources()
    order = [taxon for receipt in result["source_shard_receipts"] for taxon in receipt["taxa"]]
    observed = species.set_index("inat_taxon_id").loc[order, "rho_primary"].to_numpy()
    rng = np.random.default_rng(202609070903)
    fixed_draws = observed[rng.integers(363, size=(4999, 363))].mean(axis=1)
    interval = np.quantile(fixed_draws, [.025, .975])
    np.testing.assert_allclose(interval, result["conditional_species_bootstrap_primary_mean_95pct"], atol=2e-12, rtol=0)
    np.testing.assert_allclose(interval, [.017005540177786112, .034332127466449106], atol=2e-12, rtol=0)
    assert "not spatially or phylogenetically independent" in result["bootstrap_ceiling"]


def test_primary_replication_does_not_override_failed_conjunction():
    result, species, _ = sources()
    assert result["directional_photo_association_replicated"] is True
    assert result["flower_specific_robust_replication"] is False
    assert not all(row["positive_supported"] for row in result["metrics"].values())
    assert result["cause_or_shared_boundary_inferred"] is False
    # Direct and optimized rank paths can differ at machine precision (observed
    # maximum 1.11e-16); exact artifact hashes above are never tolerance-based.
    np.testing.assert_allclose(species.rho_primary, species.rho_calendar_quarter_stratification,
                               atol=2e-12, rtol=0)
    # The matched JSD differential is not a subtraction of two rank correlations.
    assert not np.allclose(species.rho_matched_background_differential,
                           species.rho_flower12_descriptive - species.rho_background12_descriptive)


def test_complete_tensor_artifact_audit_and_direct_checks_reconciled():
    result, species, _ = sources()
    audit = json.loads((ROOT / "docs/supporting/rgfca_reserve_inference_artifact_audit_v1.json").read_text())
    assert audit["result_sha256"] == RESULT_SHA and audit["result_commit"] == REV
    assert audit["null_tensor_shape"] == [363, 4, 999] and audit["shards_verified"] == 20
    assert audit["metrics"] == result["metrics"]
    assert audit["direct_checks"] == result["direct_checks"] == species.direct_checks.sum() == 2541
    assert audit["maximum_direct_check_error"] == result["maximum_direct_check_error"] == species.maximum_direct_check_error.max()
    assert result["maximum_direct_check_error"] <= 2e-12
    assert audit["new_permutations_or_measurements"] is False
    assert audit["selected_or_replaced_species"] is False
    assert audit["byte_hashes_verified_exactly"] is True


@pytest.mark.parametrize("path", ["README.md", "docs/RGFCA_MANUSCRIPT.md",
                                  "docs/RGFCA_RESEARCH_STATUS.md", "docs/RGFCA_SUPPORTING_EVIDENCE.md",
                                  "docs/RGFCA_RESERVE_REPLICATION_RESULTS.md"])
def test_active_documents_retain_positive_and_unsupported_result_together(path):
    text = (ROOT / path).read_text(encoding="utf-8")
    for token in ("363", "20,903", "0.001", "0.087", "flower-specific", "background causation"):
        assert token in text, (path, token)
    assert "reserve outcomes remain unopened" not in text
    assert "independent replication pending" not in text


def test_background_non_support_is_not_discovery_technical_stop():
    manuscript = (ROOT / "docs/RGFCA_MANUSCRIPT.md").read_text(encoding="utf-8")
    assert "not_evaluable_incomplete_exact_background_recovery" in manuscript
    assert "No background-adjusted statistic or p-value was computed" in manuscript
    assert "all four fixed tests subsequently" in manuscript
    assert "flower_specific_robust_replication = false" in manuscript
