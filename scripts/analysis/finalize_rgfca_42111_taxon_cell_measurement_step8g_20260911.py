#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/rgfca_42111_taxon_cell_measurement_contract_v2.json"
BREADTH_RESULT = ROOT / "results/rgfca_42111_breadth_measurement_step8f_20260911/result.json"
BIOLOGICAL = ["white", "yellow_orange", "red_pink", "blue_purple"]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--new-terminal", type=Path, required=True)
    p.add_argument("--new-terminal-result", type=Path, required=True)
    p.add_argument("--new-metadata-join", type=Path, required=True)
    p.add_argument("--reuse-key", type=Path, required=True)
    p.add_argument("--unresolved-key", type=Path, required=True)
    p.add_argument("--breadth-table", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    new_result = json.loads(args.new_terminal_result.read_text(encoding="utf-8"))
    breadth_result = json.loads(BREADTH_RESULT.read_text(encoding="utf-8"))

    if contract.get("status") != "frozen_before_any_taxon_cell_specific_image_pixel_v2":
        raise RuntimeError("v2 taxon-cell contract is not frozen")
    if new_result.get("status") != "complete_58431_new_taxon_cell_blind_terminal_measurement_without_biological_join":
        raise RuntimeError("new-only blind measurement is incomplete")
    if int(new_result.get("terminal_partition_receipts", -1)) != 128 or int(new_result.get("terminal_rows", -1)) != 58431:
        raise RuntimeError("new-only terminal denominator incomplete")
    required_breadth_status = contract["execution_order"]["stage_b_final_taxon_cell_join"]["required_breadth_status"]
    if breadth_result.get("status") != required_breadth_status:
        raise RuntimeError("breadth measurement prerequisite is incomplete")
    if int(breadth_result.get("final_rows", -1)) != 42111 or int(breadth_result.get("terminal_partition_receipts", -1)) != 128:
        raise RuntimeError("breadth denominator incomplete")
    if breadth_result.get("species_metadata_join_opened_only_after_complete_measurement") is not True:
        raise RuntimeError("breadth result did not satisfy its complete-join firewall")

    new_terminal = pd.read_csv(args.new_terminal, dtype={"measurement_id": str}).fillna("")
    new_meta = pd.read_csv(args.new_metadata_join, dtype={"measurement_id": str}).fillna("")
    reuse_key = pd.read_csv(args.reuse_key).fillna("")
    unresolved = pd.read_csv(args.unresolved_key).fillna("")
    breadth = pd.read_csv(args.breadth_table).fillna("")

    if len(new_terminal) != 58431 or new_terminal["measurement_id"].nunique() != 58431:
        raise RuntimeError("new terminal denominator drift")
    if len(new_meta) != 58431 or new_meta["measurement_id"].nunique() != 58431:
        raise RuntimeError("new metadata denominator drift")
    if set(new_terminal["measurement_id"].astype(str)) != set(new_meta["measurement_id"].astype(str)):
        raise RuntimeError("new terminal IDs do not match sealed metadata IDs")
    if len(reuse_key) != 26904 or len(unresolved) != 2:
        raise RuntimeError("reuse/unresolved accounting drift")
    if len(breadth) != 42111 or breadth["inat_taxon_id"].nunique() != 42111:
        raise RuntimeError("breadth table denominator drift")

    new_joined = new_meta.merge(new_terminal, on="measurement_id", how="left", validate="one_to_one")
    if new_joined["measurement_status"].astype(str).eq("").all():
        raise RuntimeError("new terminal measurement columns failed to join")
    new_joined["measurement_origin"] = "new_taxon_cell_blind"

    breadth_cols = [
        "observation_id", "photo_id", "measurement_id", "morph", "measurement_status",
        "roi_status", "failure_reasons", "flower_effective_pixels"
    ]
    missing_breadth = [c for c in breadth_cols if c not in breadth.columns]
    if missing_breadth:
        raise RuntimeError(f"breadth table lacks required terminal columns: {missing_breadth}")
    reusable = breadth[breadth_cols].copy().rename(columns={"measurement_id": "breadth_measurement_id"})
    reusable["observation_id"] = pd.to_numeric(reusable["observation_id"], errors="coerce")
    reusable["photo_id"] = pd.to_numeric(reusable["photo_id"], errors="coerce")
    reuse_key["observation_id"] = pd.to_numeric(reuse_key["observation_id"], errors="raise")
    reuse_key["photo_id"] = pd.to_numeric(reuse_key["photo_id"], errors="raise")
    reuse_joined = reuse_key.merge(reusable, on=["observation_id", "photo_id"], how="left", validate="one_to_one")
    if reuse_joined["measurement_status"].astype(str).eq("").any():
        raise RuntimeError("one or more frozen reuse rows lack a breadth terminal measurement")
    reuse_joined["measurement_id"] = reuse_joined["breadth_measurement_id"].astype(str)
    reuse_joined["measurement_origin"] = "reused_breadth_terminal"

    unresolved = unresolved.copy()
    unresolved["measurement_id"] = ""
    unresolved["morph"] = ""
    unresolved["measurement_status"] = "transport_unresolved_pre_pixel"
    unresolved["roi_status"] = ""
    unresolved["failure_reasons"] = unresolved["resolution_status"].astype(str)
    unresolved["flower_effective_pixels"] = ""
    unresolved["measurement_origin"] = "transport_unresolved"

    all_cols = sorted(set(new_joined.columns) | set(reuse_joined.columns) | set(unresolved.columns))
    full = pd.concat(
        [
            new_joined.reindex(columns=all_cols),
            reuse_joined.reindex(columns=all_cols),
            unresolved.reindex(columns=all_cols),
        ],
        ignore_index=True,
        sort=False,
    )
    full["inat_taxon_id"] = pd.to_numeric(full["inat_taxon_id"], errors="raise").astype(int)
    full["cell_id"] = pd.to_numeric(full["cell_id"], errors="raise").astype(int)
    full["observation_id"] = pd.to_numeric(full["observation_id"], errors="coerce").astype("Int64")
    full["photo_id"] = pd.to_numeric(full["photo_id"], errors="coerce").astype("Int64")

    if len(full) != 85337:
        raise RuntimeError(f"final taxon-cell rows {len(full)} != 85337")
    if full[["inat_taxon_id", "cell_id"]].drop_duplicates().shape[0] != 85337:
        raise RuntimeError("final table is not exactly one row per frozen taxon-cell")
    origin_counts = full["measurement_origin"].value_counts().to_dict()
    if int(origin_counts.get("new_taxon_cell_blind", 0)) != 58431:
        raise RuntimeError("new-origin accounting drift")
    if int(origin_counts.get("reused_breadth_terminal", 0)) != 26904:
        raise RuntimeError("reuse-origin accounting drift")
    if int(origin_counts.get("transport_unresolved", 0)) != 2:
        raise RuntimeError("unresolved-origin accounting drift")

    full["taxon_cell_classifiable"] = (
        full["morph"].astype(str).isin(BIOLOGICAL)
        & full["measurement_status"].astype(str).eq("classified_four_state_morph")
    )
    full = full.sort_values(["cell_id", "inat_taxon_id"], kind="mergesort").reset_index(drop=True)

    cell_rows: list[dict[str, object]] = []
    for cell_id, g in full.groupby("cell_id", sort=True, observed=True):
        classifiable = g["taxon_cell_classifiable"].astype(bool)
        n_total = int(len(g))
        n_class = int(classifiable.sum())
        vc = g.loc[classifiable, "morph"].astype(str).value_counts()
        row: dict[str, object] = {
            "cell_id": int(cell_id),
            "taxon_denominator": n_total,
            "classifiable_taxa": n_class,
            "unclassifiable_taxa": int(n_total - n_class),
            "classifiable_fraction": float(n_class / n_total) if n_total else None,
        }
        for morph in BIOLOGICAL:
            row[f"p_{morph}"] = float(vc.get(morph, 0) / n_class) if n_class else None
        cell_rows.append(row)
    cells = pd.DataFrame(cell_rows)
    if len(cells) != 128:
        raise RuntimeError(f"occupied cell count {len(cells)} != 128")
    if int(cells["taxon_denominator"].sum()) != 85337:
        raise RuntimeError("cell denominators do not sum to the frozen taxon-cell denominator")

    full_path = args.output_dir / "taxon_cell_measured_85337.csv.gz"
    cell_path = args.output_dir / "cell_species_equal_colour_composition.csv"
    full.to_csv(full_path, index=False, compression="gzip", lineterminator="\n")
    cells.to_csv(cell_path, index=False, lineterminator="\n")

    classifiable_n = int(full["taxon_cell_classifiable"].sum())
    result = {
        "protocol": contract["protocol"],
        "status": "complete_85337_taxon_cell_measurement_and_postcomplete_join",
        "taxon_cell_rows": 85337,
        "species_universe": 42111,
        "occupied_cells": 128,
        "new_image_terminal_rows": 58431,
        "reused_breadth_terminal_rows": 26904,
        "transport_unresolved_rows": 2,
        "classifiable_taxon_cell_rows": classifiable_n,
        "classifiable_fraction": float(classifiable_n / 85337),
        "measurement_origin_counts": {str(k): int(v) for k, v in origin_counts.items()},
        "new_taxon_cell_metadata_join_opened_only_after_128_partitions": True,
        "breadth_reuse_join_opened_only_after_complete_breadth": True,
        "reused_image_pixels_reopened": False,
        "pool_all_85337_as_global_frequency": False,
        "global_species_equal_frequency_source": contract["geographic_estimand"]["global_species_equal_frequency_source"],
        "geographic_estimand": contract["geographic_estimand"]["primary"],
        "hard_nonclaims": contract["hard_nonclaims"],
        "lineage": {
            "contract_sha256": sha256_file(CONTRACT),
            "new_terminal_sha256": sha256_file(args.new_terminal),
            "new_terminal_result_sha256": sha256_file(args.new_terminal_result),
            "new_metadata_join_sha256": sha256_file(args.new_metadata_join),
            "reuse_key_sha256": sha256_file(args.reuse_key),
            "unresolved_key_sha256": sha256_file(args.unresolved_key),
            "breadth_result_sha256": sha256_file(BREADTH_RESULT),
            "breadth_table_sha256": sha256_file(args.breadth_table),
            "final_table_sha256": sha256_file(full_path),
            "cell_summary_sha256": sha256_file(cell_path),
        }
    }
    (args.output_dir / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (args.output_dir / "RESULT.md").write_text(
        "# RGFCA Step 8G — complete 42,111-species taxon×cell atlas\n\n"
        "- frozen taxon×cell denominator: **85,337**\n"
        "- new blind image measurements: **58,431**\n"
        "- exact breadth measurements reused without reopening pixels: **26,904**\n"
        "- transport-unresolved rows retained as missing: **2**\n"
        "- occupied equal-area cells: **128**\n"
        f"- classifiable taxon×cell anchors: **{classifiable_n:,} / 85,337 ({classifiable_n/85337:.3%})**\n"
        "- global pooled 85,337-row colour frequency: **not reported**; the global species-equal estimator remains the separate 42,111-species breadth layer.\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
