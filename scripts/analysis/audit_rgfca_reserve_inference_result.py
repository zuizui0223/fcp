"""Reconcile every fixed inference shard, then reproduce the recorded decisions.

Uses stored null draws only. No images, new permutations or alternative tests.

Historical run-specific audit, not an importable analysis module. Run from the
repository root with PYTHONPATH=. and without Python -O. Download run 34178957447
artifact rgfca-reserve-replication-full-result-v1 into
.artifacts/reserve-inference-34178957447/final/; put each of its 20
reserve-inference-shard-N artifacts in .../shards/reserve-inference-shard-N/.
Requires full Git history, NumPy, pandas and SciPy. The frozen inference and
source files must be unchanged. Only the audit receipt is written.
"""
import hashlib
import io
import json
from pathlib import Path
import subprocess

if not __debug__:
    raise RuntimeError("Run this assertion-based artifact audit without Python -O")

import numpy as np
import pandas as pd

from fcp_pipeline.rgfca_reserve_inference import METRICS, summarize

ROOT = Path.cwd()
BASE = ROOT / ".artifacts/reserve-inference-34178957447"
REV = "ef00a78a2eda7bc79f7e6b88f719a8a696dc8b43"
AUTH = "391caecaa1a6d3c3e2407f5f19d0bad7057922e0"


def git_bytes(path, rev=REV):
    return subprocess.check_output(["git", "show", f"{rev}:{path}"])


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


result_raw = (BASE / "final/result.json").read_bytes()
assert result_raw == git_bytes("docs/supporting/rgfca_reserve_replication_result_v1.json")
result = json.loads(result_raw)
assert result["status"] == "complete_reserve_replication_with_fixed_controls"
assert result["github_run_id"] == "34178957447"
assert result["measurement_run_id"] == "34091091640"
assert result["primary_computed"] is True
auth = json.loads(git_bytes("docs/supporting/rgfca_reserve_inference_authorization_v1.json", AUTH))
for path, expected in auth["files_sha256"].items():
    assert sha(git_bytes(path, AUTH)) == expected, path
    assert git_bytes(path) == git_bytes(path, AUTH), path
lineage = result["lineage"]
for key, path in (
    ("contract_sha256", "docs/supporting/rgfca_reserve_replication_contract_v1.json"),
    ("metadata_audit_sha256", "docs/supporting/rgfca_reserve_metadata_audit_v1.json"),
    ("measured_sha256", "data/derived/rgfca_reserve_replication_measured_photos_v1.csv"),
    ("measurement_result_sha256", "docs/supporting/rgfca_reserve_replication_measurement_result_v1.json"),
):
    assert lineage[key] == auth["files_sha256"][path]
for path, expected in lineage["code_sha256"].items():
    assert expected == auth["files_sha256"][path]

support = pd.read_csv(io.BytesIO(git_bytes("data/derived/rgfca_reserve_replication_species_support_v1.csv")))
eligible = support.loc[support.evaluable].set_index("inat_taxon_id")
taxa = sorted(eligible.index.tolist())
assert len(taxa) == result["eligible_species"] == 363
assert int(eligible.classifiable_photos.sum()) == result["eligible_photos"] == 20903
assert result["raw_photos"] == 50000 and result["raw_species"] == 500
shards_root = BASE / "shards"
assert {p.name for p in shards_root.iterdir()} == {f"reserve-inference-shard-{i}" for i in range(20)}
tables, tensors = [], []
for index in range(20):
    folder = shards_root / f"reserve-inference-shard-{index}"
    receipt = json.loads((folder / "receipt.json").read_bytes())
    assert receipt == result["source_shard_receipts"][index]
    assert receipt["status"] == "complete_reserve_species_shard"
    assert receipt["lineage"] == lineage
    assert receipt["shard_index"] == index and receipt["shards"] == 20
    assert receipt["permutations"] == 999 and receipt["metrics"] == list(METRICS)
    assert receipt["eligible_species"] == 363 and receipt["github_run_id"] == "34178957447"
    assert set(receipt["files_sha256"]) == {"species.csv", "null.npz"}
    for name, expected in receipt["files_sha256"].items():
        assert sha((folder / name).read_bytes()) == expected
    table = pd.read_csv(folder / "species.csv")
    expected_taxa = taxa[index::20]
    assert table.inat_taxon_id.tolist() == receipt["taxa"] == expected_taxa
    assert table.photos.tolist() == eligible.loc[expected_taxa].classifiable_photos.tolist()
    with np.load(folder / "null.npz", allow_pickle=False) as stored:
        assert stored["taxon_ids"].tolist() == expected_taxa
        tensor = stored["null"].copy()
    assert tensor.shape == (len(expected_taxa), 4, 999) and np.isfinite(tensor).all()
    tensors.append(tensor)
    tables.append(table)

table = pd.concat(tables, ignore_index=True)
null = np.concatenate(tensors)
assert table.inat_taxon_id.nunique() == len(table) == 363
assert null.shape == (363, 4, 999)
for name, repository_path in (
    ("species.csv", "data/derived/rgfca_reserve_replication_species_v1.csv"),
    ("global_null.csv", "data/derived/rgfca_reserve_replication_global_null_v1.csv"),
):
    raw = (BASE / "final" / name).read_bytes()
    assert sha(raw) == result["files_sha256"][name]
    assert raw == git_bytes(repository_path)
final_table = pd.read_csv(BASE / "final/species.csv")
pd.testing.assert_frame_equal(table.sort_values("inat_taxon_id").reset_index(drop=True), final_table,
                              check_exact=False, atol=2e-12, rtol=0)
final_null = pd.read_csv(BASE / "final/global_null.csv")
assert final_null.permutation_index.tolist() == list(range(999))
np.testing.assert_allclose(final_null[list(METRICS)].to_numpy(), null.mean(axis=0).T, atol=2e-12, rtol=0)
rebuilt = summarize(table, null)  # Same fixed 4,999-draw conditional bootstrap; no new permutation.
for metric in METRICS:
    for key, value in rebuilt["metrics"][metric].items():
        if isinstance(value, float):
            assert abs(result["metrics"][metric][key] - value) <= 2e-12, (metric, key)
        else:
            assert result["metrics"][metric][key] == value, (metric, key)
    assert result["metrics"][metric]["p_upper"] == rebuilt["metrics"][metric]["p_upper"]
for key in ("directional_photo_association_replicated", "flower_specific_robust_replication", "cause_or_shared_boundary_inferred"):
    assert result[key] == rebuilt[key]
np.testing.assert_allclose(result["conditional_species_bootstrap_primary_mean_95pct"],
                           rebuilt["conditional_species_bootstrap_primary_mean_95pct"], atol=2e-12, rtol=0)
assert result["direct_checks"] == int(table.direct_checks.sum()) == 2541
assert result["maximum_direct_check_error"] == float(table.maximum_direct_check_error.max())
assert result["maximum_direct_check_error"] <= 2e-12
audit = {
    "status": "verified_complete_fixed_reserve_inference_artifacts",
    "run_id": 34178957447, "authorization_commit": AUTH, "result_commit": REV,
    "shards_verified": 20, "species_verified": 363, "photos_verified": 20903,
    "null_tensor_shape": list(null.shape), "global_null_draws_per_metric": 999,
    "direct_checks": 2541, "maximum_direct_check_error": result["maximum_direct_check_error"],
    "metrics": rebuilt["metrics"],
    "directional_photo_association_replicated": rebuilt["directional_photo_association_replicated"],
    "flower_specific_robust_replication": rebuilt["flower_specific_robust_replication"],
    "conditional_species_bootstrap_primary_mean_95pct": rebuilt["conditional_species_bootstrap_primary_mean_95pct"],
    "byte_hashes_verified_exactly": True,
    "numeric_reconstruction_absolute_tolerance": 2e-12,
    "new_permutations_or_measurements": False,
    "selected_or_replaced_species": False,
    "cause_or_shared_boundary_inferred": False,
    "result_sha256": sha(result_raw),
}
(ROOT / "docs/supporting/rgfca_reserve_inference_artifact_audit_v1.json").write_text(
    json.dumps(audit, indent=2) + "\n", encoding="utf-8", newline="\n")
print(json.dumps(audit, indent=2))
