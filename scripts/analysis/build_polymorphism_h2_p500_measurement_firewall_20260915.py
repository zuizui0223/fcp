#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from fcp_pipeline.photo_first_measurement_execution import (
    WORKER_FIELDS,
    ACQUISITION_FIELDS,
    semantic_shard,
    compute_partition,
)

ROOT = Path(__file__).resolve().parents[2]
PROTOCOL = ROOT / "docs" / "POLYMORPHISM_H2_P500_PROSPECTIVE_MEASUREMENT_PROTOCOL_20260915.md"
METADATA = ROOT / "data" / "frozen" / "polymorphism_h2_p500_candidate_metadata_v1.csv.gz"
FULL100 = ROOT / "results" / "polymorphism_h2_p500_candidate_metadata_20260913" / "full100_species.csv"
PARENT_RESULT = ROOT / "results" / "polymorphism_h2_p500_candidate_metadata_20260913" / "result.json"

EXPECTED_P500_SHA = "f54d07fb2338a20a3f92808b58f0ebc948ab0d5af3cb01fb26702f861184b7a4"
EXPECTED_METADATA_SHA = "a2339a3eba7bec8c29e726edc8c71a64a1a98b5cad4367764f7466c436e9f595"
EXPECTED_FULL100_SHA = "568a64e692dc21f4a5f76c5eadb615a63b07d71f01bc9cd91e8496ccd0fa3f9e"
EXPECTED_SPECIES = 499
TARGET_PER_SPECIES = 100
EXPECTED_ROWS = EXPECTED_SPECIES * TARGET_PER_SPECIES
BLIND_SALT = "FCP_H2_P500_PROSPECTIVE_20260915_V1"
BLIND_BATCHES = 2
SEMANTIC_SHARDS = 32
COMPUTE_PARTITIONS = 4


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _blind(label: str, value: object, *, length: int = 24) -> str:
    payload = f"{BLIND_SALT}\x1f{label}\x1f{value}".encode("utf-8")
    return hashlib.sha256(payload).hexdigest().upper()[:length]


def measurement_id(photo_id: object) -> str:
    return "FCPH2-" + _blind("photo", photo_id)


def measurement_batch(mid: str) -> int:
    digest = hashlib.sha256(f"fcp-h2-p500-batch\x1f{mid}".encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big") % BLIND_BATCHES


def bool_series(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.strip().str.lower().isin(["true", "1", "yes"])


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    for path in [PROTOCOL, METADATA, FULL100, PARENT_RESULT]:
        if not path.exists():
            raise RuntimeError(f"missing frozen input: {path}")

    if sha256_file(METADATA) != EXPECTED_METADATA_SHA:
        raise RuntimeError("P500 metadata SHA drifted")
    if sha256_file(FULL100) != EXPECTED_FULL100_SHA:
        raise RuntimeError("P500 full100 species SHA drifted")

    parent = json.loads(PARENT_RESULT.read_text(encoding="utf-8"))
    if parent.get("decision", {}).get("verdict") != "P500_CANDIDATE_METADATA_GATE_PASS":
        raise RuntimeError("parent metadata gate did not pass")
    if parent.get("decision", {}).get("pixels_may_be_authorized_in_separate_next_step") is not True:
        raise RuntimeError("parent metadata result does not authorize a later pixel step")
    if parent.get("P500", {}).get("sha256") != EXPECTED_P500_SHA:
        raise RuntimeError("P500 identity hash drifted")
    if int(parent.get("capacity", {}).get("full100_species", -1)) != EXPECTED_SPECIES:
        raise RuntimeError("parent full100 species denominator drifted")
    if int(parent.get("capacity", {}).get("selected_metadata_rows_full100_species", -1)) != EXPECTED_ROWS:
        raise RuntimeError("parent full100 row denominator drifted")
    firewall = parent.get("outcome_firewall", {})
    forbidden_open = ["image_pixels_opened", "flower_colour_opened", "morph_opened", "palette_opened", "D_opened", "H2_W_opened"]
    if any(firewall.get(k) is not False for k in forbidden_open):
        raise RuntimeError("parent outcome firewall was already opened")

    full = pd.read_csv(FULL100)
    if len(full) != EXPECTED_SPECIES:
        raise RuntimeError("full100 list row count drifted")
    sp_col_full = "species" if "species" in full.columns else "query_species"
    full_species = set(full[sp_col_full].astype(str))
    if len(full_species) != EXPECTED_SPECIES:
        raise RuntimeError("full100 species identities are not unique")

    metadata = pd.read_csv(METADATA, low_memory=False)
    if "full_target_species" not in metadata.columns:
        raise RuntimeError("metadata lacks pre-outcome full_target_species flag")
    metadata = metadata.loc[bool_series(metadata["full_target_species"])].copy()
    sp_col = "query_species" if "query_species" in metadata.columns else "species"
    if sp_col not in metadata.columns:
        raise RuntimeError("metadata lacks species identity")
    metadata["species"] = metadata[sp_col].astype(str)

    if len(metadata) != EXPECTED_ROWS:
        raise RuntimeError(f"pixel-opening denominator drifted: {len(metadata)} != {EXPECTED_ROWS}")
    if metadata["species"].nunique() != EXPECTED_SPECIES:
        raise RuntimeError("pixel-opening species denominator drifted")
    if set(metadata["species"]) != full_species:
        raise RuntimeError("metadata species set differs from frozen full100 set")
    counts = metadata.groupby("species", observed=True).size()
    if not (counts.astype(int) == TARGET_PER_SPECIES).all():
        raise RuntimeError("not every pixel-opening species has exactly 100 frozen rows")

    required = {
        "species", "inat_taxon_id", "observation_id", "photo_id", "photo_url_large",
        "photo_license", "latitude", "longitude",
    }
    missing = sorted(required - set(metadata.columns))
    if missing:
        raise RuntimeError(f"metadata lacks measurement fields: {missing}")
    if metadata["observation_id"].astype(str).duplicated().any():
        raise RuntimeError("duplicate observation IDs in P500 measurement denominator")
    if metadata["photo_id"].astype(str).duplicated().any():
        raise RuntimeError("duplicate photo IDs in P500 measurement denominator")

    metadata = metadata.reset_index(drop=True)
    metadata.insert(0, "measurement_id", [measurement_id(x) for x in metadata["photo_id"]])
    metadata.insert(1, "measurement_batch", [measurement_batch(x) for x in metadata["measurement_id"].astype(str)])
    if metadata["measurement_id"].nunique() != EXPECTED_ROWS:
        raise RuntimeError("blinded measurement IDs are not unique")

    worker_all = pd.DataFrame({
        "measurement_id": metadata["measurement_id"].astype(str),
        "image_filename": metadata["measurement_id"].astype(str) + ".jpg",
        "photo_license": metadata["photo_license"].astype(str),
    })
    acquisition_all = pd.DataFrame({
        "measurement_id": worker_all["measurement_id"],
        "image_filename": worker_all["image_filename"],
        "photo_url_large": metadata["photo_url_large"].astype(str),
        "photo_license": worker_all["photo_license"],
    })
    if tuple(worker_all.columns) != WORKER_FIELDS or tuple(acquisition_all.columns) != ACQUISITION_FIELDS:
        raise RuntimeError("validated blind worker/acquisition interface drifted")

    # The join key remains sealed until every terminal result is complete.
    metadata_join = metadata.drop(columns=["photo_url_large"]).copy()

    out = args.output_dir
    (out / "sealed_keys").mkdir(parents=True, exist_ok=True)
    metadata_join.to_csv(out / "sealed_keys" / "metadata_join_key.csv", index=False, lineterminator="\n")

    assignments = []
    batch_counts: dict[str, int] = {}
    nonempty = 0
    max_partition_rows = 0
    for batch in range(BLIND_BATCHES):
        keep = metadata["measurement_batch"].astype(int).eq(batch).to_numpy()
        worker = worker_all.loc[keep].reset_index(drop=True)
        acquisition = acquisition_all.loc[keep].reset_index(drop=True)
        batch_dir = out / f"batch_{batch}"
        (batch_dir / "worker_packet").mkdir(parents=True, exist_ok=True)
        (batch_dir / "sealed_keys").mkdir(parents=True, exist_ok=True)
        worker.to_csv(batch_dir / "worker_packet" / "measurement_manifest.csv", index=False, lineterminator="\n")
        acquisition.to_csv(batch_dir / "sealed_keys" / "acquisition_key.csv", index=False, lineterminator="\n")
        batch_counts[str(batch)] = int(len(worker))

        local = pd.DataFrame({
            "measurement_id": worker["measurement_id"].astype(str),
            "measurement_batch": batch,
            "semantic_shard": [semantic_shard(x, SEMANTIC_SHARDS) for x in worker["measurement_id"].astype(str)],
            "compute_partition": [compute_partition(x, COMPUTE_PARTITIONS) for x in worker["measurement_id"].astype(str)],
        })
        assignments.append(local)
        c = local.groupby(["semantic_shard", "compute_partition"], observed=True).size()
        nonempty += int(len(c))
        if len(c):
            max_partition_rows = max(max_partition_rows, int(c.max()))

    assignment = pd.concat(assignments, ignore_index=True).sort_values("measurement_id", kind="mergesort")
    assignment.to_csv(out / "partition_assignments.csv", index=False, lineterminator="\n")
    if sum(batch_counts.values()) != EXPECTED_ROWS:
        raise RuntimeError("blind batching changed denominator")

    result = {
        "analysis": "polymorphism_h2_p500_measurement_firewall",
        "status": "frozen_before_p500_pixels",
        "protocol": str(PROTOCOL.relative_to(ROOT)),
        "parent_metadata_result": str(PARENT_RESULT.relative_to(ROOT)),
        "lineage": {
            "P500_sha256": EXPECTED_P500_SHA,
            "metadata_sha256": EXPECTED_METADATA_SHA,
            "full100_species_sha256": EXPECTED_FULL100_SHA,
            "protocol_sha256": sha256_file(PROTOCOL),
        },
        "frozen_species": EXPECTED_SPECIES,
        "frozen_rows": EXPECTED_ROWS,
        "target_rows_per_species": TARGET_PER_SPECIES,
        "measurement_id_salt": BLIND_SALT,
        "blind_batches": BLIND_BATCHES,
        "rows_per_blind_batch": batch_counts,
        "semantic_shards_per_batch": SEMANTIC_SHARDS,
        "compute_partitions_per_semantic_shard": COMPUTE_PARTITIONS,
        "total_terminal_partitions": BLIND_BATCHES * SEMANTIC_SHARDS * COMPUTE_PARTITIONS,
        "nonempty_terminal_partitions": nonempty,
        "maximum_rows_in_terminal_partition": max_partition_rows,
        "worker_fields": list(worker_all.columns),
        "acquisition_fields": list(acquisition_all.columns),
        "worker_contains_species": False,
        "worker_contains_coordinates": False,
        "worker_contains_source_url": False,
        "coordinate_colour_join_opened": False,
        "candidate_pixels_opened": False,
        "H2_W_opened": False,
    }
    (out / "measurement_firewall_manifest.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
