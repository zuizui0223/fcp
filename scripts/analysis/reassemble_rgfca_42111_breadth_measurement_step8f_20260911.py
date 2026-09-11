#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from fcp_pipeline.photo_first_measurement_execution import reassemble_complete_measurement

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/rgfca_42111_breadth_measurement_contract_v1.json"
BIOLOGICAL = ["white", "yellow_orange", "red_pink", "blue_purple"]
CHECKPOINTS = [500, 1000, 2000, 5000, 10000, 20000, 30000, 42111]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--results-dir", type=Path, required=True)
    p.add_argument("--worker-manifest", type=Path, required=True)
    p.add_argument("--metadata-join-key", type=Path, required=True)
    p.add_argument("--unresolved-terminal", type=Path, required=True)
    p.add_argument("--firewall-manifest", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def composition(frame: pd.DataFrame) -> dict[str, float | int | None]:
    classifiable = frame["breadth_classifiable"].astype(bool)
    n = int(classifiable.sum())
    out: dict[str, float | int | None] = {"rows": int(len(frame)), "classifiable": n}
    if n == 0:
        for morph in BIOLOGICAL:
            out[f"p_{morph}"] = None
        return out
    vc = frame.loc[classifiable, "morph"].value_counts()
    for morph in BIOLOGICAL:
        out[f"p_{morph}"] = float(vc.get(morph, 0) / n)
    return out


def vector(summary: dict[str, float | int | None]) -> np.ndarray:
    vals = [summary[f"p_{m}"] for m in BIOLOGICAL]
    return np.asarray([0.0 if v is None else float(v) for v in vals], dtype=float)


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    firewall = json.loads(args.firewall_manifest.read_text(encoding="utf-8"))
    if contract.get("status") != "authorize_exactly_one_breadth_measurement_after_transport_gate_before_pixels":
        raise RuntimeError("breadth measurement contract is not frozen")
    if firewall.get("status") != "breadth_measurement_firewall_frozen_before_pixels":
        raise RuntimeError("breadth firewall status mismatch")
    if int(firewall.get("resolved_measurement_rows", -1)) != 42110:
        raise RuntimeError("resolved measurement denominator drift")
    if int(firewall.get("terminal_partitions", -1)) != 128:
        raise RuntimeError("terminal partition count drift")

    worker = pd.read_csv(args.worker_manifest, dtype=str).fillna("")
    metadata_join = pd.read_csv(args.metadata_join_key, dtype={"measurement_id": str}).fillna("")
    unresolved = pd.read_csv(args.unresolved_terminal, dtype={"measurement_id": str}).fillna("")
    if len(worker) != 42110 or worker["measurement_id"].nunique() != 42110:
        raise RuntimeError("worker denominator incomplete")
    if len(metadata_join) != 42110 or metadata_join["measurement_id"].nunique() != 42110:
        raise RuntimeError("metadata join denominator incomplete")
    if len(unresolved) != 1:
        raise RuntimeError("unresolved terminal denominator must equal one")

    result_pattern = re.compile(r"^partition_s(\d{2})_p(\d{2})\.csv$")
    receipt_pattern = re.compile(r"^partition_s(\d{2})_p(\d{2})\.json$")
    result_files: list[Path] = []
    result_keys: set[tuple[int, int]] = set()
    receipt_keys: set[tuple[int, int]] = set()
    receipt_rows = 0
    for path in sorted(args.results_dir.iterdir()):
        m = result_pattern.match(path.name)
        if m:
            key = (int(m.group(1)), int(m.group(2)))
            if key in result_keys:
                raise RuntimeError(f"duplicate result partition {key}")
            result_keys.add(key)
            result_files.append(path)
            continue
        m = receipt_pattern.match(path.name)
        if m:
            key = (int(m.group(1)), int(m.group(2)))
            if key in receipt_keys:
                raise RuntimeError(f"duplicate receipt partition {key}")
            receipt_keys.add(key)
            receipt = json.loads(path.read_text(encoding="utf-8"))
            if receipt.get("status") != "complete_random_photo_first_terminal_partition":
                raise RuntimeError(f"bad partition receipt status {key}")
            if receipt.get("source_urls_present") is not False or receipt.get("species_present") is not False or receipt.get("coordinates_present") is not False:
                raise RuntimeError(f"metadata leak in terminal receipt {key}")
            receipt_rows += int(receipt.get("terminal_rows", 0))

    expected = {(s, p) for s in range(32) for p in range(4)}
    if result_keys != expected or receipt_keys != expected:
        raise RuntimeError(f"incomplete terminal partition set: missing_results={sorted(expected-result_keys)} missing_receipts={sorted(expected-receipt_keys)}")
    if receipt_rows != 42110:
        raise RuntimeError(f"terminal receipt rows {receipt_rows} != 42110")

    parts = [pd.read_csv(path, dtype={"measurement_id": str}).fillna("") for path in result_files]
    assembled = reassemble_complete_measurement(
        parts,
        worker,
        metadata_join,
        expected_partition_receipts=128,
    )
    measured = assembled.joined_photos.copy()
    if len(measured) != 42110 or measured["inat_taxon_id"].nunique() != 42110:
        raise RuntimeError("joined measurement denominator incomplete")

    # Append the one pre-pixel unresolved species only after complete terminal coverage.
    all_cols = sorted(set(measured.columns) | set(unresolved.columns))
    measured = measured.reindex(columns=all_cols)
    unresolved = unresolved.reindex(columns=all_cols)
    joined = pd.concat([measured, unresolved], ignore_index=True, sort=False)
    joined["inat_taxon_id"] = pd.to_numeric(joined["inat_taxon_id"], errors="raise").astype(int)
    joined["breadth_rank"] = pd.to_numeric(joined["breadth_rank"], errors="raise").astype(int)
    if len(joined) != 42111 or joined["inat_taxon_id"].nunique() != 42111:
        raise RuntimeError("final species denominator is not exactly 42,111")
    if sorted(joined["breadth_rank"].tolist()) != list(range(1, 42112)):
        raise RuntimeError("breadth rank denominator drift")

    joined["breadth_classifiable"] = (
        joined["morph"].astype(str).isin(BIOLOGICAL)
        & joined["measurement_status"].astype(str).eq("classified_four_state_morph")
    )
    joined = joined.sort_values("breadth_rank", kind="mergesort").reset_index(drop=True)

    full = composition(joined)
    full_vec = vector(full)
    saturation_rows: list[dict[str, object]] = []
    for n in CHECKPOINTS:
        q = joined.loc[joined["breadth_rank"] <= n].copy()
        s = composition(q)
        v = vector(s)
        denom = float(np.linalg.norm(v) * np.linalg.norm(full_vec))
        cosine = float(np.dot(v, full_vec) / denom) if denom > 0 else float("nan")
        l1 = float(np.abs(v - full_vec).sum())
        row: dict[str, object] = {
            "species_prefix": n,
            "rows": int(len(q)),
            "classifiable": int(s["classifiable"]),
            "classifiable_fraction": float(int(s["classifiable"]) / len(q)),
            "cosine_to_full": cosine,
            "l1_to_full": l1,
        }
        for morph in BIOLOGICAL:
            row[f"p_{morph}"] = s[f"p_{morph}"]
        saturation_rows.append(row)
    saturation = pd.DataFrame(saturation_rows)

    morph_counts = {str(k): int(v) for k, v in joined["morph"].astype(str).value_counts(dropna=False).to_dict().items()}
    status_counts = {str(k): int(v) for k, v in joined["measurement_status"].astype(str).value_counts(dropna=False).to_dict().items()}
    classifiable_n = int(joined["breadth_classifiable"].sum())
    result = {
        "protocol": contract["protocol"],
        "status": "complete_42111_species_breadth_measurement_and_postcomplete_join",
        "species_denominator": 42111,
        "image_measurement_rows": 42110,
        "preexisting_unresolved_terminal_rows": 1,
        "terminal_partition_receipts": 128,
        "terminal_receipt_rows": receipt_rows,
        "final_rows": int(len(joined)),
        "classifiable_rows": classifiable_n,
        "classifiable_fraction": float(classifiable_n / 42111),
        "mixed_or_unresolved_rows": int(42111 - classifiable_n),
        "morph_counts": morph_counts,
        "measurement_status_counts": status_counts,
        "species_metadata_join_opened_only_after_complete_measurement": True,
        "primary_species_equal_observed_colour_composition": {m: full[f"p_{m}"] for m in BIOLOGICAL},
        "estimand_boundary": contract["estimand_boundary"],
        "no_species_modal_colour_claim": True,
        "no_polymorphism_prevalence_claim": True,
        "no_Cstar_Sstar_claim": True,
        "lineage": {
            "outer_contract_sha256": sha256_file(CONTRACT),
            "firewall_manifest_sha256": sha256_file(args.firewall_manifest),
            "worker_manifest_sha256": sha256_file(args.worker_manifest),
            "metadata_join_sha256": sha256_file(args.metadata_join_key),
            "unresolved_terminal_sha256": sha256_file(args.unresolved_terminal),
        }
    }

    measured_path = args.output_dir / "rgfca_42111_species_breadth_measured.csv.gz"
    saturation_path = args.output_dir / "breadth_prefix_saturation.csv"
    joined.to_csv(measured_path, index=False, compression="gzip", lineterminator="\n")
    saturation.to_csv(saturation_path, index=False, lineterminator="\n")
    result["lineage"]["measured_table_sha256"] = sha256_file(measured_path)
    result["lineage"]["saturation_sha256"] = sha256_file(saturation_path)
    (args.output_dir / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    comp_text = ", ".join(
        f"{m}={result['primary_species_equal_observed_colour_composition'][m]:.4f}"
        if result['primary_species_equal_observed_colour_composition'][m] is not None else f"{m}=NA"
        for m in BIOLOGICAL
    )
    (args.output_dir / "RESULT.md").write_text(
        "# RGFCA Step 8F — 42,111-species breadth colour measurement\n\n"
        f"- final species denominator: **42,111**\n"
        f"- image measurement rows: **42,110**\n"
        f"- pre-pixel unresolved species retained: **1**\n"
        f"- terminal partitions: **128 / 128**\n"
        f"- classifiable anchors: **{classifiable_n:,} / 42,111 ({classifiable_n/42111:.3%})**\n"
        f"- species-equal observed-state composition among classifiable anchors: **{comp_text}**\n"
        "- interpretation ceiling: one frozen observed flower-colour state per discovered species; not modal species colour or polymorphism prevalence.\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
