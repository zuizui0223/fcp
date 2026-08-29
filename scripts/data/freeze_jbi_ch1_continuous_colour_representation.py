#!/usr/bin/env python3
"""Freeze species-specific continuous colour vectors and scaling before evaluation opens."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

PROTOCOL = "jbi-ch1-continuous-colour-representation-v1"
EXPECTED_COUNTS = {
    "Antirrhinum majus": 60,
    "Dactylorhiza sambucina": 67,
    "Gentiana lutea": 23,
    "Ipomoea purpurea": 62,
    "Lysimachia arvensis": 66,
    "Raphanus sativus": 48,
}
FEATURE_SPECS = {
    "Antirrhinum majus": ("candidate_scores", ["magenta_pseudomajus_like", "yellow_striatum_like"]),
    "Dactylorhiza sambucina": ("candidate_scores", ["yellow", "purple"]),
    "Gentiana lutea": ("candidate_scores", ["yellow", "orange"]),
    "Ipomoea purpurea": ("candidate_scores", ["white", "pink", "blue_purple"]),
    "Lysimachia arvensis": ("candidate_scores", ["blue", "red"]),
    "Raphanus sativus": ("visual_colour_axes", ["anthocyanin_like_signal", "carotenoid_like_signal"]),
}
EXCLUSION_FIELDS = (
    "rescue_segment",
    "invalid",
    "ambiguous",
    "senescent",
    "damaged",
    "mixed_or_ambiguous",
    "not_evaluable",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_jsonl(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def usable_fresh_ordinals(spec: dict, n: int = 80) -> list[int]:
    excluded: set[int] = set()
    for field in EXCLUSION_FIELDS:
        excluded.update(int(x) for x in spec.get(field, []))
    return [i for i in range(1, n + 1) if i not in excluded]


def vector(row: dict, species: str) -> list[float]:
    source, names = FEATURE_SPECS[species]
    payload = row.get(source)
    if not isinstance(payload, dict):
        raise ValueError(f"{species}: missing {source}")
    out = [float(payload[name]) for name in names]
    if not np.isfinite(out).all():
        raise ValueError(f"{species}: non-finite feature vector")
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", type=Path, default=Path("data/calibration/jbi_ch1_florence_calibration_features_v1.jsonl"))
    parser.add_argument("--review", type=Path, default=Path("docs/supporting/jbi_ch1_blind_roi_condition_review_r1_v1.json"))
    parser.add_argument("--opening-contract", type=Path, default=Path("docs/supporting/jbi_ch1_evaluation_opening_contract_v1.json"))
    parser.add_argument("--output", type=Path, default=Path("docs/supporting/jbi_ch1_continuous_colour_representation_v1.json"))
    args = parser.parse_args()

    rows = load_jsonl(args.features)
    review = json.loads(args.review.read_text(encoding="utf-8"))
    opening = json.loads(args.opening_contract.read_text(encoding="utf-8"))
    if len(rows) != 480 or len({str(r["photo_id"]) for r in rows}) != 480:
        raise ValueError("expected exactly 480 unique calibration feature rows")
    if any(r.get("evaluation_row") is not False or r.get("final_label") is not False for r in rows):
        raise ValueError("calibration feature firewall violation")
    if opening.get("status") != "frozen_before_evaluation_image_feature_extraction":
        raise ValueError("evaluation opening contract is not frozen")
    if review.get("n_rows_reviewed") != 480 or review.get("final_label") is not False:
        raise ValueError("reviewer-1 manifest contract violated")

    per_species = {}
    for species in sorted(FEATURE_SPECS):
        group = sorted([r for r in rows if r["species"] == species], key=lambda r: str(r["blind_id"]))
        if len(group) != 80:
            raise ValueError(f"{species}: expected 80 calibration rows, found {len(group)}")
        ordinals = usable_fresh_ordinals(review["species"][species])
        if len(ordinals) != EXPECTED_COUNTS[species]:
            raise ValueError(f"{species}: expected {EXPECTED_COUNTS[species]} usable+fresh rows, found {len(ordinals)}")
        selected = [group[i - 1] for i in ordinals]
        if any(r.get("feature_status") != "ok" for r in selected):
            raise ValueError(f"{species}: selected calibration row lacks numeric feature")
        x = np.asarray([vector(r, species) for r in selected], dtype=float)
        mean = x.mean(axis=0)
        scale = x.std(axis=0, ddof=0)
        if np.any(~np.isfinite(mean)) or np.any(~np.isfinite(scale)) or np.any(scale <= 1e-12):
            raise ValueError(f"{species}: invalid calibration scaling")
        source, names = FEATURE_SPECS[species]
        per_species[species] = {
            "feature_source": source,
            "feature_names": names,
            "n_scaler_rows": len(selected),
            "mean": [float(v) for v in mean],
            "scale_population_sd": [float(v) for v in scale],
            "reviewer1_usable_fresh_ordinals_sha256": hashlib.sha256(
                ",".join(map(str, ordinals)).encode("utf-8")
            ).hexdigest(),
        }

    result = {
        "protocol": PROTOCOL,
        "status": "frozen_before_evaluation_values_inspected",
        "primary_representation": "species-specific continuous colour vector",
        "calibration_only_scaler_fit": True,
        "evaluation_values_used_to_define_axes": False,
        "evaluation_values_used_to_fit_scaling": False,
        "candidate_argmax_is_primary": False,
        "discrete_states_required_for_primary_analysis": False,
        "source_calibration_feature_sha256": sha256(args.features),
        "source_reviewer1_manifest_sha256": sha256(args.review),
        "source_evaluation_opening_contract_sha256": sha256(args.opening_contract),
        "n_scaler_rows_total": sum(v["n_scaler_rows"] for v in per_species.values()),
        "per_species": per_species,
        "standardization": "z=(raw-calibration_mean)/calibration_population_sd; parameters remain fixed for evaluation",
        "spatial_edge_metric": "RMS Euclidean difference across standardized features",
        "binary_special_case": "one-dimensional 0/1 vector reduces exactly to discrete mismatch indicator",
    }
    if result["n_scaler_rows_total"] != 326:
        raise ValueError("expected 326 total scaler rows")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
