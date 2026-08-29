#!/usr/bin/env python3
"""Re-run colour feature geometry on reviewer-1 usable+fresh calibration rows.

This is deliberately provisional: reviewer-1 is not an independent final adjudication.
The script tests whether obvious ROI/condition contamination changes the within-species
feature geometry before investing in a second blinded review. It never emits colour
labels and never opens evaluation rows.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from sklearn.preprocessing import StandardScaler

PROTOCOL = "jbi-ch1-reviewer1-fresh-feature-geometry-v1"
EXPECTED_COUNTS = {
    "Antirrhinum majus": 60,
    "Dactylorhiza sambucina": 67,
    "Gentiana lutea": 23,
    "Ipomoea purpurea": 62,
    "Lysimachia arvensis": 66,
    "Raphanus sativus": 48,
}


def load_geometry_module():
    path = Path(__file__).with_name("analyze_jbi_ch1_calibration_feature_geometry.py")
    spec = importlib.util.spec_from_file_location("base_geometry", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("could not load base feature geometry implementation")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def reviewer1_usable_fresh_ordinals(spec: dict, n: int = 80) -> list[int]:
    excluded = set()
    for field in (
        "rescue_segment",
        "invalid",
        "ambiguous",
        "senescent",
        "damaged",
        "mixed_or_ambiguous",
        "not_evaluable",
    ):
        excluded.update(int(x) for x in spec.get(field, []))
    return [ordinal for ordinal in range(1, n + 1) if ordinal not in excluded]


def select_rows(features: list[dict], review: dict) -> list[dict]:
    if review.get("calibration_only") is not True or review.get("evaluation_rows_opened") is not False:
        raise ValueError("reviewer-1 calibration/evaluation firewall violation")
    if review.get("final_label") is not False:
        raise ValueError("reviewer-1 unexpectedly contains final labels")
    if review.get("independent_second_review_completed") is not False:
        raise ValueError("this diagnostic expects reviewer-1-only state")

    selected = []
    review_species = review.get("species", {})
    if set(review_species) != set(EXPECTED_COUNTS):
        raise ValueError(f"reviewer-1 species set mismatch: {sorted(review_species)}")

    for species in sorted(EXPECTED_COUNTS):
        group = sorted(
            [row for row in features if str(row.get("species")) == species],
            key=lambda row: str(row["blind_id"]),
        )
        if len(group) != 80:
            raise ValueError(f"{species}: expected 80 rows, found {len(group)}")
        ordinals = reviewer1_usable_fresh_ordinals(review_species[species])
        manifest_n = int(review_species[species]["usable_fresh_n"])
        if len(ordinals) != manifest_n or len(ordinals) != EXPECTED_COUNTS[species]:
            raise ValueError(
                f"{species}: reviewer-1 usable+fresh count mismatch: derived={len(ordinals)} manifest={manifest_n} expected={EXPECTED_COUNTS[species]}"
            )
        for ordinal in ordinals:
            row = dict(group[ordinal - 1])
            if row.get("feature_status") != "ok":
                raise ValueError(f"{species} ordinal {ordinal}: reviewer-1 usable+fresh row lacks numeric feature")
            row["reviewer1_order_within_species"] = ordinal
            selected.append(row)

    expected_total = sum(EXPECTED_COUNTS.values())
    if expected_total != 326 or len(selected) != expected_total:
        raise ValueError(f"expected 326 reviewer-1 usable+fresh rows, found {len(selected)}")
    if int(review.get("overall", {}).get("usable_fresh_n", -1)) != expected_total:
        raise ValueError("reviewer-1 overall usable+fresh count does not equal 326")
    if len({str(row["blind_id"]) for row in selected}) != len(selected):
        raise ValueError("duplicate blind IDs in reviewer-1 filtered subset")
    if any(row.get("evaluation_row") is not False or row.get("final_label") is not False for row in selected):
        raise ValueError("evaluation/final-label firewall violation in filtered subset")
    return selected


def analyze_species(module, rows: list[dict], species: str) -> dict:
    spec = module.SPECIES_SPECS[species]
    group = [row for row in rows if str(row["species"]) == species]
    if len(group) != EXPECTED_COUNTS[species]:
        raise ValueError(f"{species}: filtered count mismatch")
    raw = np.asarray([module.feature_vector(row, species) for row in group], dtype=float)
    scaled = StandardScaler().fit_transform(raw)
    max_components = min(int(spec["max_components"]), max(1, len(group) // 10))
    grid = module.fit_bic_grid(scaled, max_components)
    return {
        "species": species,
        "status": "reviewer1_filtered_geometry_measured_not_morph_labels",
        "n_reviewer1_usable_fresh": len(group),
        "feature_source": spec["source"],
        "features": spec["features"],
        "literature_candidate_component_cap": spec["max_components"],
        "components_evaluated": list(range(1, max_components + 1)),
        "bic_grid": grid,
        "bic_selected_components": module.best_bic_components(grid),
        "bootstrap_bic_selected_component_frequency": module.bootstrap_component_support(
            scaled, max_components
        ),
        "raw_feature_quantiles": module.quantiles(raw, list(spec["features"])),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", type=Path, default=Path("data/calibration/jbi_ch1_florence_calibration_features_v1.jsonl"))
    parser.add_argument("--review", type=Path, default=Path("docs/supporting/jbi_ch1_blind_roi_condition_review_r1_v1.json"))
    parser.add_argument("--pre-filter-geometry", type=Path, default=Path("docs/supporting/jbi_ch1_calibration_feature_geometry_v1.json"))
    parser.add_argument("--output", type=Path, default=Path("docs/supporting/jbi_ch1_reviewer1_fresh_feature_geometry_v1.json"))
    args = parser.parse_args()

    module = load_geometry_module()
    features = module.load_rows(args.features)
    review = json.loads(args.review.read_text(encoding="utf-8"))
    selected = select_rows(features, review)
    pre = json.loads(args.pre_filter_geometry.read_text(encoding="utf-8"))
    pre_by_species = {row["species"]: row for row in pre["species"]}

    species_results = []
    for species in sorted(EXPECTED_COUNTS):
        result = analyze_species(module, selected, species)
        before = pre_by_species[species]
        result["pre_filter_bic_selected_components"] = before.get("bic_selected_components")
        result["pre_filter_bootstrap_bic_selected_component_frequency"] = before.get(
            "bootstrap_bic_selected_component_frequency"
        )
        species_results.append(result)

    output = {
        "protocol": PROTOCOL,
        "status": "reviewer1_filtered_provisional_geometry_complete_not_final_labels",
        "calibration_only": True,
        "evaluation_rows_opened": False,
        "final_label": False,
        "n_source_rows": 480,
        "n_reviewer1_usable_fresh_rows": len(selected),
        "reviewer1_filter_applied": True,
        "independent_second_review_completed": False,
        "reviewer1_decisions_are_final": False,
        "geography_used": False,
        "observer_used": False,
        "date_used": False,
        "environment_used": False,
        "colour_candidate_argmax_used_as_training_label": False,
        "gmm_components_are_biological_morph_labels": False,
        "bootstrap_replicates": module.BOOTSTRAPS,
        "species": species_results,
        "interpretation_limit": "This diagnostic measures sensitivity to reviewer-1 ROI/condition filtering only. It cannot freeze colour states until an independent second review is completed.",
        "next_gate": "independent blinded ROI/condition validation, then repeat geometry on the independently validated fresh/evaluable calibration subset"
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
