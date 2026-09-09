#!/usr/bin/env python3
"""Pre-inference provenance audit for polymorphism directionality Step 1."""
from __future__ import annotations

import csv
import json
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "polymorphism_directionality_step1_input_audit_20260909"
OUT.mkdir(parents=True, exist_ok=True)
TABLES = {
    "discovery": ROOT / "data" / "derived" / "global_monte_carlo_measured_photos_v1.csv",
    "reserve": ROOT / "data" / "derived" / "rgfca_reserve_replication_measured_photos_v1.csv",
}
MORPHS = ("white", "yellow_orange", "red_pink", "blue_purple")


def as_bool(x: str) -> bool:
    return str(x).strip().lower() in {"true", "1", "yes"}


def git_grep(pattern: str) -> list[str]:
    proc = subprocess.run(
        ["git", "grep", "-n", "-I", "-E", pattern, "--", ":(exclude)data/**"],
        cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
    )
    return [x for x in proc.stdout.splitlines() if x.strip()][:500]


def fmt(value, digits: int = 4) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def summarize_table(path: Path) -> dict:
    if not path.exists():
        return {"exists": False, "path": str(path.relative_to(ROOT))}
    per_species_morph: dict[str, Counter] = defaultdict(Counter)
    per_species_terminal: dict[str, Counter] = defaultdict(Counter)
    n_rows = n_classifiable = 0
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        columns = reader.fieldnames or []
        morph_field = "global_morph" if "global_morph" in columns else ("morph" if "morph" in columns else None)
        for row in reader:
            n_rows += 1
            species = (row.get("species") or "").strip()
            morph = (row.get(morph_field) or "").strip() if morph_field else ""
            if not species:
                continue
            classifiable = as_bool(row.get("global_classifiable", "")) if "global_classifiable" in columns else (
                morph in MORPHS and (not row.get("measurement_status") or row.get("measurement_status") == "classified_four_state_morph")
            )
            per_species_terminal[species][morph or (row.get("measurement_status") or "blank")] += 1
            if classifiable and morph in MORPHS:
                per_species_morph[species][morph] += 1
                n_classifiable += 1

    palette_count = [c for c in columns if c.startswith("palette_count_")]
    background_palette_count = [c for c in columns if c.startswith("background_palette_count_")]
    flower_fraction = [c for c in columns if c.startswith("flower_fraction_")]
    background_fraction = [c for c in columns if c.startswith("background_fraction_")]
    colour4 = [c for c in columns if c.startswith("colour_")]

    eligible = {}
    for sp, cnt in per_species_morph.items():
        n = sum(cnt.values())
        if n < 40:
            continue
        p = sorted((cnt.get(m, 0) / n for m in MORPHS), reverse=True)
        eligible[sp] = {
            "n_four_state": n,
            "D": 1.0 - sum(x * x for x in p),
            "second_fraction": p[1],
            "counts": dict(cnt),
        }
    dvals = [v["D"] for v in eligible.values()]
    seconds = [v["second_fraction"] for v in eligible.values()]
    return {
        "exists": True,
        "path": str(path.relative_to(ROOT)),
        "rows": n_rows,
        "classifiable_rows_reconstructed": n_classifiable,
        "columns": columns,
        "morph_field": morph_field,
        "has_global_classifiable": "global_classifiable" in columns,
        "palette_count_columns": palette_count,
        "background_palette_count_columns": background_palette_count,
        "flower_fraction_columns": flower_fraction,
        "background_fraction_columns": background_fraction,
        "colour4_columns": colour4,
        "n_species_any_four_state": len(per_species_morph),
        "n_species_ge40_four_state": len(eligible),
        "D_min": min(dvals) if dvals else None,
        "D_max": max(dvals) if dvals else None,
        "second_ge_0_10_fraction": sum(x >= 0.10 for x in seconds) / len(seconds) if seconds else None,
        "second_ge_0_20_fraction": sum(x >= 0.20 for x in seconds) / len(seconds) if seconds else None,
        "eligible_species": eligible,
    }


def fingerprint_distance(x: dict) -> float | None:
    vals = (x.get("n_species_ge40_four_state"), x.get("D_max"), x.get("second_ge_0_10_fraction"), x.get("second_ge_0_20_fraction"))
    if any(v is None for v in vals):
        return None
    n, dmax, s10, s20 = vals
    return abs(n - 369) / 369 + abs(dmax - 0.71) + abs(s10 - 0.466) + abs(s20 - 0.266)


def main() -> None:
    report = {
        "audit": "polymorphism_directionality_step1_input_provenance",
        "date_jst": "2026-09-09",
        "inference_opened": False,
        "reported_fingerprint": {
            "n_species": 369, "D_range": [0.0, 0.71],
            "second_ge_0_10_fraction": 0.466, "second_ge_0_20_fraction": 0.266,
            "reported_correlations": {"D_SD_M1": 0.332, "D_SD_M2": 0.299, "D_mixed_rate": 0.537},
        },
        "tables": {},
    }
    for name, path in TABLES.items():
        report["tables"][name] = summarize_table(path)
        report["tables"][name]["fingerprint_distance"] = fingerprint_distance(report["tables"][name])
    report["repository_grep"] = {
        "exact_reported_numbers": git_grep(r"0\.332|0\.299|0\.537"),
        "simpson_polymorphism_terms": git_grep(r"[Ss]impson|polymorph|second[_ -]?morph|morph[_ -]?divers"),
        "mode_sd_terms": git_grep(r"SD\(M[123]\)|sd_m[123]|mode.*sd|species.*mode.*std"),
    }
    candidates = sorted((x["fingerprint_distance"], name) for name, x in report["tables"].items() if x.get("fingerprint_distance") is not None)
    report["fingerprint_best_table"] = candidates[0][1] if candidates else None
    report["fingerprint_best_distance"] = candidates[0][0] if candidates else None
    (OUT / "input_audit.json").write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")

    lines = [
        "# Polymorphism directionality Step 1 — input provenance audit", "",
        "This audit was run before the Step 1 inferential statistic. `inference_opened = false`.", "",
        "## Fingerprint reconstruction", "",
        "| table | rows | classifiable | morph field | species >=40 | D max | second >=10% | second >=20% | flower 12 | background 12 | colour4 |",
        "|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for name in ("discovery", "reserve"):
        x = report["tables"][name]
        lines.append(
            f"| {name} | {x.get('rows','NA')} | {x.get('classifiable_rows_reconstructed','NA')} | {x.get('morph_field','NA')} | "
            f"{x.get('n_species_ge40_four_state','NA')} | {fmt(x.get('D_max'),6)} | {fmt(x.get('second_ge_0_10_fraction'))} | "
            f"{fmt(x.get('second_ge_0_20_fraction'))} | {len(x.get('palette_count_columns',[]))} | "
            f"{len(x.get('background_palette_count_columns',[]))} | {len(x.get('colour4_columns',[]))} |"
        )
    lines += ["", f"Fingerprint-best table: **{report['fingerprint_best_table']}** (distance {fmt(report['fingerprint_best_distance'],6)}).", ""]
    for key, hits in report["repository_grep"].items():
        lines += [f"### {key}"] + ([f"- `{h}`" for h in hits[:80]] if hits else ["- no tracked-text hit"]) + [""]
    lines += ["## Decision boundary", "", "No Step 1 effect estimate or p-value is computed by this audit. Any source-cohort correction must be committed before the inferential runner is executed.", ""]
    (OUT / "input_audit.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines[:40]))


if __name__ == "__main__":
    main()
