#!/usr/bin/env python3
"""Header-only audit of existing species-level covariates before Step 4 inference."""
from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "polymorphism_step4_covariate_audit_20260909"
OUT.mkdir(parents=True, exist_ok=True)
TERMS = {
    "taxonomy": ["family", "genus", "order", "taxon"],
    "life_form": ["life_form", "lifeform", "growth_form", "growthform", "woodiness", "annual", "perennial", "herb", "shrub", "tree"],
    "pollination": ["pollination", "pollinator", "pollinat", "syndrome", "breeding_system", "self_compat"],
    "range": ["span", "range", "radius", "extent", "area"],
    "latitude": ["latitude", "lat", "centroid"],
}


def header(path: Path):
    try:
        with path.open(newline="", encoding="utf-8") as fh:
            return next(csv.reader(fh))
    except Exception:
        return []


def main():
    matches = {k: [] for k in TERMS}
    scanned = 0
    for root in [ROOT / "data" / "frozen", ROOT / "data" / "derived"]:
        if not root.exists():
            continue
        for path in sorted(root.rglob("*.csv")):
            cols = header(path)
            if not cols:
                continue
            scanned += 1
            low = [c.lower() for c in cols]
            for group, terms in TERMS.items():
                hit = [cols[i] for i, c in enumerate(low) if any(t in c for t in terms)]
                if hit:
                    matches[group].append({"path": str(path.relative_to(ROOT)), "columns": hit})

    known = {
        "discovery_measurement": "data/derived/global_monte_carlo_measured_photos_v1.csv",
        "capacity_species": "data/frozen/global_monte_carlo_capacity_scan_selected_species_v3.csv",
        "candidate_species_audit": "data/frozen/global_monte_carlo_candidate_species_audit_v1.csv",
    }
    report = {
        "audit": "polymorphism_step4_covariate_availability",
        "date_jst": "2026-09-09",
        "outcomes_opened_by_this_audit": False,
        "csv_headers_scanned": scanned,
        "matches": matches,
        "known_sources": known,
        "preliminary_availability": {
            "genus": "derivable_from_binomial_species_name_without_external_lookup",
            "family": "available_only_if_a_matching tracked taxonomy column/file appears above; otherwise missing",
            "life_form": "available_only_if a matching tracked trait column/file appears above; otherwise missing",
            "pollination_mode": "available_only_if a matching tracked trait column/file appears above; otherwise missing",
            "range_size_proxy": "available_from frozen metadata span fields and/or raw coordinates; must freeze one definition before D analysis",
            "latitude_centroid": "derivable from frozen colour-blind sampled coordinates; must freeze raw-photo rather than classifiable-photo denominator before D analysis",
        },
    }
    (OUT / "audit.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    lines = [
        "# Step 4 covariate availability audit",
        "",
        "Header-only audit; no D/covariate association was computed.",
        "",
        f"CSV headers scanned: **{scanned}**",
        "",
    ]
    for group in ["taxonomy", "life_form", "pollination", "range", "latitude"]:
        lines.append(f"## {group}")
        hits = matches[group]
        if not hits:
            lines.append("- no tracked CSV header match")
        else:
            for h in hits[:50]:
                lines.append(f"- `{h['path']}`: {', '.join(h['columns'])}")
        lines.append("")
    lines += [
        "## Freeze implications",
        "",
        "- genus can be derived from the frozen binomial species name.",
        "- range-size and latitude-centroid proxies can be built from colour-blind frozen metadata; their exact estimators must be frozen before opening D associations.",
        "- family, life form, and pollination mode require a tracked source if none is listed above; they must not be back-filled after Step 4 outcomes are opened.",
    ]
    (OUT / "AUDIT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
