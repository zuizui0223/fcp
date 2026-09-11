#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/rgfca_42111_taxon_cell_measurement_contract_v2.json"
FORBIDDEN_CONTEXT = {
    "species", "inat_taxon_id", "cell_id", "observation_id", "photo_id",
    "latitude", "longitude", "observer_id", "photo_url", "photo_url_large",
    "source_url", "discovery_source"
}


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
    p.add_argument("--firewall-manifest", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    firewall = json.loads(args.firewall_manifest.read_text(encoding="utf-8"))

    if contract.get("status") != "frozen_before_any_taxon_cell_specific_image_pixel_v2":
        raise RuntimeError("v2 contract is not frozen")
    if firewall.get("status") != "taxon_cell_new_only_firewall_frozen_before_pixels_v2":
        raise RuntimeError("new-only firewall status mismatch")
    if int(firewall.get("new_image_measurement_rows", -1)) != 58431:
        raise RuntimeError("new-only measurement denominator drift")
    if int(firewall.get("new_terminal_partitions", -1)) != 128:
        raise RuntimeError("terminal partition count drift")

    worker = pd.read_csv(args.worker_manifest, dtype=str).fillna("")
    if len(worker) != 58431 or worker["measurement_id"].nunique() != 58431:
        raise RuntimeError("worker denominator incomplete")
    if FORBIDDEN_CONTEXT & set(worker.columns):
        raise RuntimeError("worker manifest contains biological context")

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
        raise RuntimeError(
            f"incomplete terminal partition set: missing_results={sorted(expected-result_keys)} "
            f"missing_receipts={sorted(expected-receipt_keys)}"
        )
    if receipt_rows != 58431:
        raise RuntimeError(f"terminal receipt rows {receipt_rows} != 58431")

    parts = [pd.read_csv(path, dtype={"measurement_id": str}).fillna("") for path in result_files]
    terminal = pd.concat(parts, ignore_index=True, sort=False)
    if len(terminal) != 58431 or terminal["measurement_id"].nunique() != 58431:
        raise RuntimeError("terminal aggregate is not exactly one row per new measurement")
    if set(terminal["measurement_id"].astype(str)) != set(worker["measurement_id"].astype(str)):
        raise RuntimeError("terminal measurement IDs do not equal the frozen worker denominator")

    leaked = FORBIDDEN_CONTEXT & set(terminal.columns)
    if leaked:
        raise RuntimeError(f"blind terminal aggregate leaked biological context: {sorted(leaked)}")
    object_cols = [c for c in terminal.columns if terminal[c].dtype == object]
    for c in object_cols:
        vals = terminal[c].astype(str)
        if vals.str.contains(r"https?://", regex=True, na=False).any():
            raise RuntimeError(f"blind terminal aggregate leaked URL-like values in {c}")

    terminal = terminal.sort_values("measurement_id", kind="mergesort").reset_index(drop=True)
    out_table = args.output_dir / "rgfca_42111_taxon_cell_new_terminal_blind.csv.gz"
    terminal.to_csv(out_table, index=False, compression="gzip", lineterminator="\n")

    morph_counts = {}
    if "morph" in terminal.columns:
        morph_counts = {str(k): int(v) for k, v in terminal["morph"].astype(str).value_counts(dropna=False).to_dict().items()}
    status_counts = {}
    if "measurement_status" in terminal.columns:
        status_counts = {str(k): int(v) for k, v in terminal["measurement_status"].astype(str).value_counts(dropna=False).to_dict().items()}

    result = {
        "protocol": contract["protocol"],
        "status": "complete_58431_new_taxon_cell_blind_terminal_measurement_without_biological_join",
        "new_image_measurement_rows": 58431,
        "terminal_partition_receipts": 128,
        "terminal_receipt_rows": int(receipt_rows),
        "terminal_rows": int(len(terminal)),
        "measurement_ids_unique": int(terminal["measurement_id"].nunique()),
        "blind_morph_counts_not_biologically_joined": morph_counts,
        "blind_measurement_status_counts": status_counts,
        "species_context_present": False,
        "cell_context_present": False,
        "source_urls_present": False,
        "taxon_cell_colour_join_opened": False,
        "crosscell_discordance_opened": False,
        "lineage": {
            "contract_sha256": sha256_file(CONTRACT),
            "firewall_manifest_sha256": sha256_file(args.firewall_manifest),
            "worker_manifest_sha256": sha256_file(args.worker_manifest),
            "terminal_table_sha256": sha256_file(out_table),
        },
        "claim_boundary": "Colour outcomes exist only under opaque measurement IDs. No species, taxon, cell, or cross-cell biological inference is authorized from this aggregate."
    }
    (args.output_dir / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (args.output_dir / "RESULT.md").write_text(
        "# RGFCA Step 8G — new-only blind taxon-cell measurement\n\n"
        "- frozen new-photo denominator: **58,431**\n"
        "- terminal partitions: **128 / 128**\n"
        "- biological join: **not opened**\n"
        "- breadth-reuse join: **deferred**\n"
        "- cross-cell discordance: **not opened**\n\n"
        "This artifact is a location-blind terminal measurement aggregate keyed only by opaque measurement IDs.\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
