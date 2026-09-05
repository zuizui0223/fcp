#!/usr/bin/env python3
"""Strictly reassemble all nine G1 shards and compute the one frozen G1 p-value."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from fcp_pipeline.global_repeated_atlas import consensus_field

ROOT = Path(__file__).resolve().parents[2]
EXECUTION = ROOT / "docs/supporting/global_monte_carlo_inference_execution_contract_v1.json"
OUT_OBSERVED = ROOT / "data/derived/global_rgfca_observed_outer_fields_v1.npz"
OUT_NULL = ROOT / "data/derived/global_rgfca_null_consensus_fields_v1.npz"
OUT_MANIFEST = ROOT / "docs/supporting/global_rgfca_g1_result_v1.json"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def array_sha256(*arrays: np.ndarray) -> str:
    h = hashlib.sha256()
    for array in arrays:
        a = np.ascontiguousarray(array)
        h.update(str(a.dtype).encode("ascii"))
        h.update(np.asarray(a.shape, dtype=np.int64).tobytes())
        h.update(a.tobytes())
    return h.hexdigest()


def field_correlation(left: np.ndarray, right: np.ndarray) -> float:
    a = np.asarray(left, dtype=float)
    b = np.asarray(right, dtype=float)
    keep = np.isfinite(a) & np.isfinite(b)
    if np.count_nonzero(keep) < 3:
        return float("nan")
    x = a[keep]
    y = b[keep]
    if np.ptp(x) <= 1e-15 or np.ptp(y) <= 1e-15:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--input-root", type=Path, required=True)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    execution = json.loads(EXECUTION.read_text(encoding="utf-8"))
    shard_count = int(execution["null_execution"]["deterministic_shards"])
    if shard_count != 9:
        raise RuntimeError("G1 reducer expects the frozen nine-shard execution")
    for path in (OUT_OBSERVED, OUT_NULL, OUT_MANIFEST):
        if path.exists():
            raise RuntimeError(f"refusing to overwrite frozen G1 output: {path}")

    manifests: list[dict[str, object]] = []
    npz_data: list[dict[str, np.ndarray]] = []
    observed_digests: set[str] = set()
    measured_hashes: set[str] = set()
    for shard in range(shard_count):
        json_matches = list(args.input_root.rglob(f"global_rgfca_g1_shard_{shard:02d}.json"))
        npz_matches = list(args.input_root.rglob(f"global_rgfca_g1_shard_{shard:02d}.npz"))
        if len(json_matches) != 1 or len(npz_matches) != 1:
            raise RuntimeError(f"G1 shard {shard} requires exactly one manifest and npz")
        manifest = json.loads(json_matches[0].read_text(encoding="utf-8"))
        if manifest.get("status") != "complete_global_rgfca_g1_null_shard":
            raise RuntimeError(f"G1 shard {shard} status invalid")
        if int(manifest.get("shard_index", -1)) != shard or int(manifest.get("shard_count", -1)) != shard_count:
            raise RuntimeError(f"G1 shard {shard} identity drift")
        if sha256_file(npz_matches[0]) != manifest["lineage"]["npz_sha256"]:
            raise RuntimeError(f"G1 shard {shard} npz SHA mismatch")
        with np.load(npz_matches[0], allow_pickle=False) as z:
            data = {name: z[name].copy() for name in z.files}
        digest = array_sha256(
            data["observed_outer_fields"],
            data["observed_outer_opportunities"],
            data["observed_consensus"],
            data["observed_aggregate_opportunity"],
            data["observed_concentration"].astype(np.float64),
        )
        if digest != manifest["observed_digest_sha256"]:
            raise RuntimeError(f"G1 shard {shard} observed digest does not match its manifest")
        if not np.array_equal(data["null_indices"].astype(int), np.asarray(manifest["null_indices"], dtype=int)):
            raise RuntimeError(f"G1 shard {shard} null index payload mismatch")
        observed_digests.add(digest)
        measured_hashes.add(str(manifest["lineage"]["measured_table_sha256"]))
        manifests.append(manifest)
        npz_data.append(data)

    if len(observed_digests) != 1:
        raise RuntimeError("observed G1 program is not bitwise identical across all nine null shards")
    if len(measured_hashes) != 1:
        raise RuntimeError("G1 shards do not share one measured-table lineage")

    all_indices = np.concatenate([data["null_indices"].astype(np.int64) for data in npz_data])
    if len(all_indices) != 999 or len(np.unique(all_indices)) != 999 or set(all_indices.tolist()) != set(range(999)):
        raise RuntimeError("G1 null shards do not contain indices 0..998 exactly once")
    all_fields = np.vstack([data["null_consensus_fields"] for data in npz_data])
    all_concentrations = np.concatenate([data["null_concentrations"].astype(float) for data in npz_data])
    order = np.argsort(all_indices, kind="mergesort")
    indices = all_indices[order]
    null_fields = all_fields[order]
    null_concentrations = all_concentrations[order]
    if not np.array_equal(indices, np.arange(999, dtype=np.int64)):
        raise RuntimeError("sorted G1 null index order drifted")
    if not np.isfinite(null_concentrations).all():
        raise RuntimeError("G1 null concentration contains non-finite values")

    first = npz_data[0]
    observed_fields = first["observed_outer_fields"].astype(float)
    observed_opportunities = first["observed_outer_opportunities"].astype(float)
    observed_consensus = first["observed_consensus"].astype(float)
    aggregate_opportunity = first["observed_aggregate_opportunity"].astype(float)
    observed_concentration = float(first["observed_concentration"][0])
    if observed_fields.shape[0] != int(execution["outer_schedule"]["observed_resamples"]):
        raise RuntimeError("observed G1 outer realization count drifted")

    exceed = int(np.count_nonzero(null_concentrations >= observed_concentration))
    p_upper = float((1 + exceed) / 1000.0)
    null_mean = float(np.mean(null_concentrations))
    alpha = float(execution["g1_primary"]["alpha"])
    supported = bool(observed_concentration > null_mean and p_upper < alpha)
    decision = "support_recurrent_global_flower_colour_barrier_concentration" if supported else "no_support_recurrent_global_flower_colour_barrier_concentration"

    checkpoints = {}
    final = consensus_field(observed_fields, observed_opportunities)
    for n in execution["g2_outputs_from_observed_outer_fields"]["running_consensus_checkpoints"]:
        partial = consensus_field(observed_fields[: int(n)], observed_opportunities[: int(n)])
        checkpoints[str(int(n))] = {
            "concentration": float(partial.concentration),
            "field_correlation_with_200": field_correlation(partial.field, final.field),
        }
    odd = consensus_field(observed_fields[0::2], observed_opportunities[0::2])
    even = consensus_field(observed_fields[1::2], observed_opportunities[1::2])
    odd_even_r = field_correlation(odd.field, even.field)

    OUT_OBSERVED.parent.mkdir(parents=True, exist_ok=True)
    OUT_MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        OUT_OBSERVED,
        observed_outer_fields=observed_fields,
        observed_outer_opportunities=observed_opportunities,
        observed_consensus=observed_consensus,
        observed_aggregate_opportunity=aggregate_opportunity,
        observed_concentration=np.asarray([observed_concentration], dtype=float),
    )
    np.savez_compressed(
        OUT_NULL,
        null_indices=indices,
        null_consensus_fields=null_fields,
        null_concentrations=null_concentrations,
    )
    manifest = {
        "protocol": execution["protocol"],
        "status": "complete_global_rgfca_g1_primary_inference",
        "decision": decision,
        "g1_supported": supported,
        "observed_concentration": observed_concentration,
        "null_mean_concentration": null_mean,
        "null_q025": float(np.quantile(null_concentrations, 0.025)),
        "null_q975": float(np.quantile(null_concentrations, 0.975)),
        "null_permutations": 999,
        "null_exceed_or_equal_count": exceed,
        "p_upper": p_upper,
        "alpha": alpha,
        "observed_above_null_mean": bool(observed_concentration > null_mean),
        "null_unit": manifests[0]["null_unit"],
        "eligible_species": int(manifests[0]["eligible_species"]),
        "classifiable_pool_rows": int(manifests[0]["classifiable_pool_rows"]),
        "schedule_audit": manifests[0]["schedule_audit"],
        "predeclared_observed_stability_diagnostics": {
            "running_consensus": checkpoints,
            "odd_even_consensus_field_pearson_r": odd_even_r,
            "scope": "partial G2 diagnostics only; leave-one-realm/family and spatial-support robustness remain separate and no persistent zone is frozen here",
        },
        "persistent_zone_extracted": False,
        "g4_overlay_run": False,
        "lineage": {
            "execution_contract_sha256": sha256_file(EXECUTION),
            "measured_table_sha256": next(iter(measured_hashes)),
            "observed_digest_sha256": next(iter(observed_digests)),
            "observed_outer_npz_sha256": sha256_file(OUT_OBSERVED),
            "null_consensus_npz_sha256": sha256_file(OUT_NULL),
            "shard_npz_sha256": [m["lineage"]["npz_sha256"] for m in manifests],
        },
        "files": {
            "observed_outer_fields": str(OUT_OBSERVED.relative_to(ROOT)),
            "null_consensus_fields": str(OUT_NULL.relative_to(ROOT)),
        },
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
