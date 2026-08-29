#!/usr/bin/env python3
"""Reconcile reviewer-1 and independently reblinded reviewer-2 ROI/condition decisions.

The reconciliation follows the protocol frozen before reviewer-2 decisions: direct
colour-calibration inclusion requires both reviewers to independently call ROI usable
and condition fresh. Any disagreement on those two binary gates is withheld for third
adjudication. No colour state is assigned here.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import pandas as pd

PROTOCOL = "jbi-ch1-r1-r2-roi-condition-reconciliation-v1"
ALLOWED_ROI = {"usable", "rescue_segment", "invalid", "ambiguous"}
ALLOWED_CONDITION = {"fresh", "senescent", "damaged", "mixed_or_ambiguous", "not_evaluable"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def reviewer1_decision(spec: dict, ordinal: int) -> tuple[str, str]:
    roi = "usable"
    if ordinal in {int(x) for x in spec.get("rescue_segment", [])}:
        roi = "rescue_segment"
    if ordinal in {int(x) for x in spec.get("invalid", [])}:
        roi = "invalid"
    if ordinal in {int(x) for x in spec.get("ambiguous", [])}:
        roi = "ambiguous"

    condition = "fresh"
    if ordinal in {int(x) for x in spec.get("senescent", [])}:
        condition = "senescent"
    if ordinal in {int(x) for x in spec.get("damaged", [])}:
        condition = "damaged"
    if ordinal in {int(x) for x in spec.get("mixed_or_ambiguous", [])}:
        condition = "mixed_or_ambiguous"
    if ordinal in {int(x) for x in spec.get("not_evaluable", [])}:
        condition = "not_evaluable"
    return roi, condition


def reviewer1_table(review: dict, features_path: Path) -> pd.DataFrame:
    rows = [json.loads(line) for line in features_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(rows) != 480:
        raise ValueError("expected 480 feature rows")
    records = []
    for species, spec in sorted(review["species"].items()):
        group = sorted([row for row in rows if str(row["species"]) == species], key=lambda row: str(row["blind_id"]))
        if len(group) != 80:
            raise ValueError(f"{species}: expected 80 rows")
        for ordinal, row in enumerate(group, start=1):
            roi, condition = reviewer1_decision(spec, ordinal)
            records.append({
                "blind_id": str(row["blind_id"]),
                "photo_id": str(row["photo_id"]),
                "species": species,
                "reviewer1_order_within_species": ordinal,
                "reviewer1_roi": roi,
                "reviewer1_condition": condition,
            })
    frame = pd.DataFrame(records)
    if len(frame) != 480 or frame["blind_id"].nunique() != 480:
        raise ValueError("reviewer1 table contract violated")
    return frame


def validate_reviewer2(review_csv: Path, mapping_csv: Path) -> pd.DataFrame:
    review = pd.read_csv(review_csv, keep_default_na=False)
    mapping = pd.read_csv(mapping_csv, keep_default_na=False)
    if len(review) != 480 or review["r2_id"].nunique() != 480:
        raise ValueError("reviewer2 queue must contain 480 unique r2 IDs")
    if len(mapping) != 480 or mapping["r2_id"].nunique() != 480 or mapping["blind_id"].nunique() != 480:
        raise ValueError("reviewer2 mapping contract violated")
    if set(review["target_roi_validity"]) - ALLOWED_ROI:
        raise ValueError("reviewer2 contains missing or invalid ROI decisions")
    if set(review["condition_review"]) - ALLOWED_CONDITION:
        raise ValueError("reviewer2 contains missing or invalid condition decisions")
    merged = mapping.merge(review, on=["review_order", "r2_id"], how="inner", validate="one_to_one")
    if len(merged) != 480:
        raise ValueError("reviewer2 mapping/review join lost records")
    return merged


def reconcile(r1: pd.DataFrame, r2: pd.DataFrame) -> pd.DataFrame:
    joined = r1.merge(
        r2[["r2_id", "blind_id", "target_roi_validity", "condition_review", "reviewer2_notes"]],
        on="blind_id",
        how="inner",
        validate="one_to_one",
    )
    if len(joined) != 480:
        raise ValueError("r1/r2 join lost records")
    joined = joined.rename(columns={
        "target_roi_validity": "reviewer2_roi",
        "condition_review": "reviewer2_condition",
    })
    joined["r1_roi_usable"] = joined["reviewer1_roi"].eq("usable")
    joined["r2_roi_usable"] = joined["reviewer2_roi"].eq("usable")
    joined["r1_fresh"] = joined["reviewer1_condition"].eq("fresh")
    joined["r2_fresh"] = joined["reviewer2_condition"].eq("fresh")
    joined["roi_usability_disagreement"] = joined["r1_roi_usable"].ne(joined["r2_roi_usable"])
    joined["fresh_condition_disagreement"] = joined["r1_fresh"].ne(joined["r2_fresh"])
    joined["third_adjudication_required"] = joined["roi_usability_disagreement"] | joined["fresh_condition_disagreement"]
    joined["direct_consensus_usable_fresh"] = (
        joined["r1_roi_usable"]
        & joined["r2_roi_usable"]
        & joined["r1_fresh"]
        & joined["r2_fresh"]
    )
    joined["consensus_excluded_without_disagreement"] = (~joined["third_adjudication_required"]) & (~joined["direct_consensus_usable_fresh"])
    if (joined["direct_consensus_usable_fresh"] & joined["third_adjudication_required"]).any():
        raise ValueError("consensus direct records cannot require adjudication")
    return joined


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reviewer1", type=Path, default=Path("docs/supporting/jbi_ch1_blind_roi_condition_review_r1_v1.json"))
    parser.add_argument("--reviewer2", type=Path, default=Path("data/calibration/jbi_ch1_reviewer2_reblind_queue_v1.csv"))
    parser.add_argument("--mapping", type=Path, default=Path("data/calibration/jbi_ch1_reviewer2_reblind_mapping_v1.csv"))
    parser.add_argument("--features", type=Path, default=Path("data/calibration/jbi_ch1_florence_calibration_features_v1.jsonl"))
    parser.add_argument("--frozen-protocol", type=Path, default=Path("docs/supporting/jbi_ch1_reviewer2_reblind_protocol_v1.json"))
    parser.add_argument("--output-csv", type=Path, default=Path("data/calibration/jbi_ch1_r1_r2_roi_condition_reconciliation_v1.csv"))
    parser.add_argument("--manifest", type=Path, default=Path("docs/supporting/jbi_ch1_r1_r2_roi_condition_reconciliation_v1.json"))
    args = parser.parse_args()

    frozen = json.loads(args.frozen_protocol.read_text(encoding="utf-8"))
    if frozen.get("status") != "frozen_before_reviewer2_decisions":
        raise ValueError("reviewer2 reconciliation protocol was not frozen before decisions")
    r1_json = json.loads(args.reviewer1.read_text(encoding="utf-8"))
    if r1_json.get("calibration_only") is not True or r1_json.get("evaluation_rows_opened") is not False:
        raise ValueError("reviewer1 firewall violation")
    r1 = reviewer1_table(r1_json, args.features)
    r2 = validate_reviewer2(args.reviewer2, args.mapping)
    out = reconcile(r1, r2)

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output_csv, index=False)
    by_species = {}
    for species, group in out.groupby("species", sort=True):
        by_species[species] = {
            "n": int(len(group)),
            "direct_consensus_usable_fresh": int(group["direct_consensus_usable_fresh"].sum()),
            "third_adjudication_required": int(group["third_adjudication_required"].sum()),
            "consensus_excluded_without_disagreement": int(group["consensus_excluded_without_disagreement"].sum()),
            "roi_usability_disagreement": int(group["roi_usability_disagreement"].sum()),
            "fresh_condition_disagreement": int(group["fresh_condition_disagreement"].sum()),
        }
    manifest = {
        "protocol": PROTOCOL,
        "status": "reviewer1_reviewer2_reconciled_pending_third_adjudication_if_needed",
        "calibration_only": True,
        "evaluation_rows_opened": False,
        "final_colour_label": False,
        "n_rows": 480,
        "n_direct_consensus_usable_fresh": int(out["direct_consensus_usable_fresh"].sum()),
        "n_third_adjudication_required": int(out["third_adjudication_required"].sum()),
        "n_consensus_excluded_without_disagreement": int(out["consensus_excluded_without_disagreement"].sum()),
        "per_species": by_species,
        "reviewer1_sha256": sha256(args.reviewer1),
        "reviewer2_sha256": sha256(args.reviewer2),
        "mapping_sha256": sha256(args.mapping),
        "frozen_protocol_sha256": sha256(args.frozen_protocol),
        "output_csv_sha256": sha256(args.output_csv),
        "colour_state_assigned": False,
        "next_gate": "complete third adjudication for disagreements, then rerun feature geometry on the final consensus fresh/evaluable calibration set"
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(manifest, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
