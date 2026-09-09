"""Closed measurement census only; no statistical inference or image acquisition.

Historical run-specific audit, not an importable analysis module. Run from the
repository root with PYTHONPATH=. and without Python -O. Under
.artifacts/reserve-complete-34091091640/, download that run's measured artifact
into measured/, firewall artifact into firewall/, and 256 terminal artifacts
into terminals/ (one named directory per artifact). Download inference preflight
34091922722 into .artifacts/reserve-inference-preflight-34091922722/.
Export the literal frozen contract from REV with core.autocrlf=false into
exact-inputs/docs/supporting/rgfca_reserve_replication_contract_v1.json under BASE.
Materialize Git blob b6e76efb4db8e7529b5ce5a91c60b6ae68637fb5 as the untracked
fcp_pipeline/photo_first_measurement_execution.py (as the frozen worker does).
Requires full Git history and the frozen reassembly dependencies. Creates
verified copies in flat-terminals/ and writes only the audit receipt. A later
replay checks the recorded pre-inference census, not a new no-peek claim.
"""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

if not __debug__:
    raise RuntimeError("Run this assertion-based artifact audit without Python -O")

import pandas as pd

from fcp_pipeline import rgfca_reserve_replication as reserve

ROOT = Path.cwd()
BASE = ROOT / ".artifacts/reserve-complete-34091091640"
REV = "9f5abe7b45fcdc20ba83adf75d1a8d4a640f622c"
MEASURED = BASE / "measured"
RESULT = "docs/supporting/rgfca_reserve_replication_measurement_result_v1.json"


def git_bytes(path, rev=REV):
    return subprocess.check_output(["git", "show", f"{rev}:{path}"])


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


preflight = json.loads((ROOT / ".artifacts/reserve-inference-preflight-34091922722/tmp/rgfca-reserve-preflight.json").read_bytes())
assert preflight["github_run_id"] == "34091922722"
assert preflight["head_sha"] == "708b305e314087047a39d27956a76c9327585497"
for path, expected in preflight["files_sha256"].items():
    assert sha(git_bytes(path)) == expected, path
    assert git_bytes(path) == git_bytes(path, preflight["head_sha"]), path
execution = ROOT / "fcp_pipeline/photo_first_measurement_execution.py"
assert subprocess.check_output(["git", "hash-object", str(execution)], text=True).strip() == "b6e76efb4db8e7529b5ce5a91c60b6ae68637fb5"

result_raw = (MEASURED / RESULT).read_bytes()
assert result_raw == git_bytes(RESULT)
result = json.loads(result_raw)
assert result["github_run_id"] == "34091091640" and result["inference_run"] is False
for path, expected in result["files_sha256"].items():
    raw = (MEASURED / path).read_bytes()
    assert raw == git_bytes(path), path
    assert sha(raw) == expected, path

# The exact frozen contract is exported using core.autocrlf=false. Only this
# audit's input path changes; neither the contract nor a verifier is rewritten.
reserve.CONTRACT = BASE / "exact-inputs/docs/supporting/rgfca_reserve_replication_contract_v1.json"
assert sha(reserve.CONTRACT.read_bytes()) == preflight["files_sha256"]["docs/supporting/rgfca_reserve_replication_contract_v1.json"]
flat = BASE / "flat-terminals"
flat.mkdir(exist_ok=True)
sources = list((BASE / "terminals").glob("*/b*_partition_s*_p*.*"))
assert len(sources) == 512 and len({p.name for p in sources}) == 512
for source in sources:
    target = flat / source.name
    if target.exists():
        assert target.read_bytes() == source.read_bytes()
    else:
        shutil.copyfile(source, target)

joined, support, rebuilt = reserve.reassemble(BASE / "firewall", flat)
for key, value in rebuilt.items():
    assert result[key] == value, key
audit = json.loads(git_bytes("docs/supporting/rgfca_reserve_metadata_audit_v1.json"))
assert {field: reserve.id_digest(joined[field].unique()) for field in reserve.ID_FIELDS} == audit["ids_sha256"]
for field in ("measurement_id", "photo_id", "observation_id"):
    assert not joined[field].duplicated().any(), field
for path in (reserve.DISCOVERY, *reserve.PRIOR):
    prior = pd.read_csv(__import__("io").BytesIO(git_bytes(path, reserve.SOURCE_COMMIT)), usecols=["photo_id", "observation_id"])
    assert reserve.overlap_counts(joined, prior) == {"photo_id": 0, "observation_id": 0}

stored = pd.read_csv(MEASURED / "data/derived/rgfca_reserve_replication_measured_photos_v1.csv")
for key in (*reserve.ID_FIELDS, "measurement_id", "global_classifiable", "measurement_status"):
    assert joined[key].tolist() == stored[key].tolist(), key
stored_support = pd.read_csv(MEASURED / "data/derived/rgfca_reserve_replication_species_support_v1.csv")
pd.testing.assert_frame_equal(support, stored_support, check_dtype=False, check_exact=True)
eligible = set(support.loc[support.evaluable, "inat_taxon_id"])
eligible_photos = int((joined.global_classifiable & joined.inat_taxon_id.isin(eligible)).sum())
receipt = {
    "status": "complete_census_verified_before_fixed_reserve_inference",
    "measurement_run_id": 34091091640,
    "measurement_head_sha": "1f80af7f5db81d61f28ae2818c130ef058a9f850",
    "result_commit": REV,
    "preflight_run_id": 34091922722,
    "preflight_head_sha": preflight["head_sha"],
    "tested_files_verified": len(preflight["files_sha256"]),
    "terminal_partitions_verified": 256,
    "terminal_files_verified": 512,
    "raw_species": len(support), "raw_photos": len(joined),
    "classifiable_photos": result["classifiable_photos"],
    "eligible_species": len(eligible), "eligible_photos": eligible_photos,
    "measurement_gate_pass": result["measurement_gate_pass"],
    "all_partition_hashes_match_completed_manifest": True,
    "photo_observation_ids_match_frozen_audit": True,
    "prior_photo_observation_overlap": 0,
    "terminal_status_counts": {str(k): int(v) for k, v in joined.measurement_status.value_counts().items()},
    "files_sha256": {**result["files_sha256"], RESULT: sha(result_raw)},
    "inference_run_by_audit": False,
    "images_reacquired": False,
    "threshold_or_cohort_changed": False,
    "interpretation": "Measurement completeness and eligibility only; no ecological effect or p-value computed.",
}
out = ROOT / "docs/supporting/rgfca_reserve_complete_census_audit_v1.json"
out.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8", newline="\n")
print(json.dumps(receipt, indent=2))
