#!/usr/bin/env python3
"""Evaluate the prespecified Wave-0 duplicate-review calibration gate."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

THRESHOLDS = {
    "record_relevance": {"raw": 0.90, "kappa": 0.60},
    "natural_intraspecific_variation": {"raw": 0.90, "kappa": 0.60},
    "floral_display_colour": {"raw": 0.90, "kappa": 0.60},
    "full_text_required": {"raw": 0.85, "kappa": 0.60},
}
EXPECTED_N = 384


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--agreement", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    data = json.loads(Path(args.agreement).read_text(encoding="utf-8"))
    checks = {}
    ready = True
    passed = True
    for field, threshold in THRESHOLDS.items():
        metric = data.get("agreement", {}).get(field, {})
        n = int(metric.get("n_double_coded") or 0)
        raw = metric.get("raw_agreement")
        kap = metric.get("cohen_kappa")
        labels = metric.get("labels") or []
        complete = n == EXPECTED_N
        raw_ok = raw is not None and float(raw) >= threshold["raw"]
        # Kappa is mathematically undefined when both reviewers use a single category.
        # In that degenerate case we report it but do not fail solely because kappa is null.
        kappa_not_estimable_single_category = kap is None and len(labels) <= 1 and complete
        kappa_ok = kappa_not_estimable_single_category or (kap is not None and float(kap) >= threshold["kappa"])
        checks[field] = {
            "n_double_coded": n,
            "required_n": EXPECTED_N,
            "raw_agreement": raw,
            "raw_threshold": threshold["raw"],
            "cohen_kappa": kap,
            "kappa_threshold": threshold["kappa"],
            "kappa_not_estimable_single_category": kappa_not_estimable_single_category,
            "complete": complete,
            "raw_ok": raw_ok,
            "kappa_ok": kappa_ok,
            "pass": complete and raw_ok and kappa_ok,
        }
        ready = ready and complete
        passed = passed and checks[field]["pass"]

    if not ready:
        status = "not_ready"
        recommendation = "Complete independent duplicate coding before interpreting calibration performance."
    elif passed:
        status = "pass"
        recommendation = "Lock the codebook version and proceed to the full 12,064-record duplicate screen."
    else:
        status = "fail"
        recommendation = (
            "Inspect disagreement patterns, revise the codebook only for a documented systematic ambiguity, "
            "then run a new blinded calibration wave before full screening."
        )

    result = {
        "status": status,
        "prespecified_before_human_results": True,
        "expected_records": EXPECTED_N,
        "checks": checks,
        "focal_taxon_text_agreement": data.get("focal_taxon_text", {}),
        "focal_taxon_text_is_diagnostic_not_gate": True,
        "recommendation": recommendation,
        "semantic_guard": (
            "This gate measures coding reproducibility. Passing does not validate biological inclusion, "
            "taxonomic identity, or spatial state."
        ),
    }
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
