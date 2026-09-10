#!/usr/bin/env python3
"""Outcome-blind schema/status audit for reserve flower-colour measurement missingness."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "data" / "derived" / "rgfca_reserve_replication_measured_photos_v1.csv"
OUT = ROOT / "results" / "polymorphism_measurement_missingness_step8_preflight_20260910"
OUT.mkdir(parents=True, exist_ok=True)

TERMS = ("status", "class", "evalu", "morph", "palette", "roi", "flip", "fail", "error", "result")


def write_json(path: Path, obj: Any) -> None:
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def value_summary(s: pd.Series) -> dict[str, Any]:
    nonmissing = s.notna() & s.astype(str).ne("")
    vals = s.loc[nonmissing].astype(str)
    counts = vals.value_counts(dropna=False)
    n_unique = int(vals.nunique(dropna=True))
    limit = None if n_unique <= 50 else 50
    shown = counts if limit is None else counts.head(limit)
    return {
        "nonmissing": int(nonmissing.sum()),
        "unique_nonmissing": n_unique,
        "counts_truncated": bool(limit is not None),
        "value_counts": {str(k): int(v) for k, v in shown.items()},
    }


def main() -> None:
    header = pd.read_csv(SOURCE, nrows=0).columns.astype(str).tolist()
    forbidden = {"D", "D_unbiased", "spatial_observed_rho", "rho_primary", "log1p_span_primary", "reserve_log1p_span", "genus"}
    if forbidden.intersection(header):
        # The raw ledger may legitimately contain geometry but must not contain species-level outcome columns.
        bad = sorted(forbidden.intersection(header))
        raise RuntimeError(f"outcome firewall: forbidden species-level analysis columns present in raw ledger: {bad}")

    selected = [c for c in header if any(term in c.lower() for term in TERMS)]
    required = ["species", "inat_taxon_id"]
    usecols = list(dict.fromkeys(required + selected))
    df = pd.read_csv(SOURCE, usecols=usecols, low_memory=False)

    if df.empty:
        raise RuntimeError("reserve measurement ledger is empty")
    if df["species"].isna().any():
        raise RuntimeError("missing species in reserve measurement ledger")

    per_species = df.groupby("species", sort=False).size()
    exact_100 = bool(per_species.eq(100).all())

    fields = {c: value_summary(df[c]) for c in selected}
    result = {
        "analysis": "polymorphism_measurement_missingness_step8_preflight",
        "date_jst": "2026-09-10",
        "protocol": "docs/POLYMORPHISM_MEASUREMENT_MISSINGNESS_STEP8_PREFLIGHT_20260910.md",
        "D_association_computed": False,
        "species_stratified_status_analysis_computed": False,
        "row_count": int(len(df)),
        "unique_species": int(df["species"].nunique()),
        "rows_per_species_min": int(per_species.min()),
        "rows_per_species_median": float(per_species.median()),
        "rows_per_species_max": int(per_species.max()),
        "every_species_exactly_100_rows": exact_100,
        "header": header,
        "candidate_status_columns": selected,
        "candidate_status_fields": fields,
    }
    write_json(OUT / "result.json", result)

    lines = [
        "# Polymorphism measurement missingness — Step 8 preflight result",
        "",
        f"- rows: **{len(df)}**",
        f"- unique species: **{df['species'].nunique()}**",
        f"- rows/species min–median–max: **{per_species.min()} – {per_species.median():.0f} – {per_species.max()}**",
        f"- every species exactly 100 rows: **{exact_100}**",
        f"- candidate status fields: **{len(selected)}**",
        "- D association computed: **False**",
        "",
        "## Candidate status fields",
        "",
    ]
    for c in selected:
        x = fields[c]
        lines.append(f"### `{c}`")
        lines.append("")
        lines.append(f"nonmissing={x['nonmissing']}; unique={x['unique_nonmissing']}; truncated={x['counts_truncated']}")
        lines.append("")
        for k, v in x["value_counts"].items():
            lines.append(f"- `{k}`: {v}")
        lines.append("")
    (OUT / "RESULT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(json.dumps({
        "row_count": result["row_count"],
        "unique_species": result["unique_species"],
        "every_species_exactly_100_rows": exact_100,
        "candidate_status_columns": selected,
    }, indent=2))


if __name__ == "__main__":
    main()
