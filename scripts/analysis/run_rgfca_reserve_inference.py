#!/usr/bin/env python3
"""Complete-census sharded inference for the fixed independent RGFCA reserve."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from fcp_pipeline.rgfca_reserve_inference import METRICS, PERMUTATIONS, evaluate_species, summarize
from fcp_pipeline.rgfca_reserve_replication import AUDIT, CONTRACT, ROOT, PREFIX, PROTOCOL, ID_FIELDS, id_digest, sha
from scripts.analysis.run_rgfca_reserve_replication import write_json

MEASURED = ROOT / f"data/derived/{PREFIX}_measured_photos_v1.csv"
MEASUREMENT_RESULT = ROOT / f"docs/supporting/{PREFIX}_measurement_result_v1.json"
SHARDS = 20
CODE = (
    "fcp_pipeline/rgfca_reserve_inference.py",
    "scripts/analysis/run_rgfca_reserve_inference.py",
    "scripts/analysis/run_global_rgfca_within_species_spatial_omnibus.py",
    "scripts/analysis/finalize_global_rgfca_background_control.py",
    "fcp_pipeline/global_g3.py",
)


def load_pool() -> tuple[pd.DataFrame, dict]:
    result = json.loads(MEASUREMENT_RESULT.read_text(encoding="utf-8"))
    audit = json.loads(AUDIT.read_text(encoding="utf-8"))
    if result.get("protocol") != PROTOCOL or result.get("status") != "complete_reserve_location_blind_measurement_and_join":
        raise ValueError("reserve measurement not complete")
    if result.get("terminal_partitions") != 256 or result.get("coordinate_colour_join_opened_after_complete_measurement") is not True:
        raise ValueError("reserve coordinate-colour join was not authorized by complete measurement")
    if result["files_sha256"][MEASURED.relative_to(ROOT).as_posix()] != sha(MEASURED.read_bytes()):
        raise ValueError("reserve measurement hash drift")
    frame = pd.read_csv(MEASURED)
    if len(frame) != 50000 or frame.inat_taxon_id.nunique() != 500 or not frame.groupby("inat_taxon_id").size().eq(100).all():
        raise ValueError("raw reserve denominator drift")
    if audit["ids_sha256"] != {f: id_digest(frame[f].unique()) for f in ID_FIELDS}:
        raise ValueError("raw reserve identities changed")
    if frame.photo_id.duplicated().any() or frame.observation_id.duplicated().any():
        raise ValueError("duplicate reserve measurements")
    if frame.groupby("species").inat_taxon_id.nunique().max() != 1 or frame.groupby("inat_taxon_id").species.nunique().max() != 1:
        raise ValueError("species label/taxon mapping is not one-to-one")
    if not frame.global_classifiable.astype(str).str.lower().isin(["true", "false"]).all():
        raise ValueError("invalid classification flag")
    pool = frame.loc[frame.global_classifiable.astype(str).str.lower().eq("true")].copy()
    counts = pool.groupby("inat_taxon_id").size()
    pool = pool.loc[pool.inat_taxon_id.isin(counts[counts >= 40].index)].copy()
    if result["evaluable_species"] != pool.inat_taxon_id.nunique():
        raise ValueError("measurement eligibility count differs from inference gate")
    if bool(pool.inat_taxon_id.nunique() >= 250) != result["measurement_gate_pass"]:
        raise ValueError("measurement eligibility decision drift")
    return pool.sort_values(["inat_taxon_id", "photo_id"], kind="stable").reset_index(drop=True), result


def lineage() -> dict:
    return {"contract_sha256": sha(CONTRACT.read_bytes()), "metadata_audit_sha256": sha(AUDIT.read_bytes()),
            "measured_sha256": sha(MEASURED.read_bytes()), "measurement_result_sha256": sha(MEASUREMENT_RESULT.read_bytes()),
            "code_sha256": {p: sha((ROOT / p).read_bytes()) for p in CODE}}


def run_shard(index: int, out: Path) -> None:
    if not 0 <= index < SHARDS:
        raise ValueError("invalid fixed species shard")
    pool, result = load_pool()
    out.mkdir(parents=True, exist_ok=True)
    receipt = {"protocol": PROTOCOL, "shard_index": index, "shards": SHARDS, "lineage": lineage(),
               "eligible_species": int(pool.inat_taxon_id.nunique()), "permutations": PERMUTATIONS,
               "metrics": list(METRICS), "github_run_id": os.environ.get("GITHUB_RUN_ID")}
    if result["measurement_gate_pass"] is not True:
        write_json(out / "receipt.json", {**receipt, "status": "not_evaluable_fewer_than_250_species"})
        return
    taxa = sorted(pool.inat_taxon_id.unique())[index::SHARDS]
    rows, distributions = [], []
    for taxon in taxa:
        g = pool.loc[pool.inat_taxon_id == taxon]
        try:
            row, null = evaluate_species(g)
        except (ValueError, RuntimeError) as exc:
            write_json(out / "receipt.json", {**receipt, "status": "not_evaluable_invalid_eligible_species",
                                              "failed_taxon": int(taxon), "reason": str(exc),
                                              "previously_computed_species_not_used": len(rows)})
            return
        rows.append(row)
        distributions.append(null)
    table = pd.DataFrame(rows)
    table.to_csv(out / "species.csv", index=False, lineterminator="\n")
    np.savez_compressed(out / "null.npz", taxon_ids=np.asarray(taxa, dtype="int64"), null=np.stack(distributions))
    write_json(out / "receipt.json", {**receipt, "status": "complete_reserve_species_shard",
                                      "taxa": list(map(int, taxa)),
                                      "files_sha256": {p: sha((out / p).read_bytes()) for p in ("species.csv", "null.npz")}})


def finalize(shards: Path, out: Path) -> None:
    pool, measured_result = load_pool()
    expected_lineage = lineage()
    directories = [shards / f"shard-{i}" for i in range(SHARDS)]
    if {p.name for p in shards.iterdir()} != {p.name for p in directories}:
        raise ValueError("incomplete or extraneous inference shards")
    receipts = [json.loads((p / "receipt.json").read_text(encoding="utf-8")) for p in directories]
    for index, receipt in enumerate(receipts):
        if (receipt.get("shard_index") != index or receipt.get("lineage") != expected_lineage
                or receipt.get("metrics") != list(METRICS) or receipt.get("permutations") != PERMUTATIONS
                or receipt.get("eligible_species") != pool.inat_taxon_id.nunique()):
            raise ValueError("inference shard provenance changed")
    out.mkdir(parents=True, exist_ok=True)
    common = {"protocol": PROTOCOL, "lineage": expected_lineage, "raw_species": 500, "raw_photos": 50000,
              "eligible_species": int(pool.inat_taxon_id.nunique()), "eligible_photos": len(pool),
              "github_run_id": os.environ.get("GITHUB_RUN_ID"),
              "measurement_run_id": measured_result["github_run_id"], "source_shard_receipts": receipts}
    failed = [r for r in receipts if r["status"] != "complete_reserve_species_shard"]
    if failed:
        write_json(out / "result.json", {**common, "status": "not_evaluable_reserve_replication",
                                        "failed_shards": [r["shard_index"] for r in failed],
                                        "primary_computed": False, "replacement_used": False})
        return
    species, tensors = [], []
    for index, (path, receipt) in enumerate(zip(directories, receipts)):
        for name, expected in receipt["files_sha256"].items():
            if sha((path / name).read_bytes()) != expected:
                raise ValueError("inference shard output hash mismatch")
        table = pd.read_csv(path / "species.csv")
        expected_taxa = list(map(int, sorted(pool.inat_taxon_id.unique())[index::SHARDS]))
        if table.inat_taxon_id.tolist() != expected_taxa or receipt["taxa"] != expected_taxa:
            raise ValueError("shard species membership mismatch")
        with np.load(path / "null.npz", allow_pickle=False) as stored:
            if stored["taxon_ids"].tolist() != expected_taxa or stored["null"].shape != (len(table), 4, 999):
                raise ValueError("null tensor census mismatch")
            tensors.append(stored["null"].copy())
        species.append(table)
    table = pd.concat(species, ignore_index=True)
    null = np.concatenate(tensors)
    if table.inat_taxon_id.nunique() != len(table) or set(table.inat_taxon_id) != set(pool.inat_taxon_id):
        raise ValueError("incomplete or duplicated final species census")
    counts = pool.groupby("inat_taxon_id").size()
    if not np.array_equal(table.photos.to_numpy(), counts.loc[table.inat_taxon_id].to_numpy()):
        raise ValueError("per-species photo counts changed")
    decision = summarize(table, null)
    table.sort_values("inat_taxon_id").to_csv(out / "species.csv", index=False, lineterminator="\n")
    global_null = pd.DataFrame(null.mean(axis=0).T, columns=METRICS)
    global_null.insert(0, "permutation_index", np.arange(999))
    global_null.to_csv(out / "global_null.csv", index=False, lineterminator="\n")
    write_json(out / "result.json", {**common, "status": "complete_reserve_replication_with_fixed_controls",
                                     "primary_computed": True, **decision,
                                     "maximum_direct_check_error": float(table.maximum_direct_check_error.max()),
                                     "direct_checks": int(table.direct_checks.sum()),
                                     "files_sha256": {p: sha((out / p).read_bytes()) for p in ("species.csv", "global_null.csv")},
                                     "claim_ceiling": "Species- and observation-disjoint replication in the same iNaturalist frame and measurement model. No causal, all-plant-prevalence or shared-boundary conclusion. Report every required control, including failures."})


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["shard", "finalize"])
    parser.add_argument("--shard-index", type=int)
    parser.add_argument("--shards-dir", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.stage == "shard":
        if args.shard_index is None:
            parser.error("shard requires --shard-index")
        run_shard(args.shard_index, args.output_dir)
    else:
        if args.shards_dir is None:
            parser.error("finalize requires --shards-dir")
        finalize(args.shards_dir, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
