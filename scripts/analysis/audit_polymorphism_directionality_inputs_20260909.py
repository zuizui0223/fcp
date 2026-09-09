#!/usr/bin/env python3
"""Pre-inference provenance audit for polymorphism directionality Step 1.

This script deliberately does not run the Step 1 inferential test. It establishes
which frozen table reconstructs the reported discrete-polymorphism fingerprint,
which continuous palette fields are available on that table, and where the
already-reported exploratory correlations were generated in the repository.
"""
from __future__ import annotations

import csv
import json
import math
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


def git_grep(pattern: str) -> list[str]:
    proc = subprocess.run(
        ["git", "grep", "-n", "-I", "-E", pattern, "--", ":(exclude)data/**"],
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    return [x for x in proc.stdout.splitlines() if x.strip()][:500]


def summarize_table(path: Path) -> dict:
    if not path.exists():
        return {"exists": False, "path": str(path.relative_to(ROOT))}
    per_species_morph: dict[str, Counter] = defaultdict(Counter)
    per_species_status: dict[str, Counter] = defaultdict(Counter)
    n_rows = 0
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        columns = reader.fieldnames or []
        for row in reader:
            n_rows += 1
            species = (row.get("species") or "").strip()
            morph = (row.get("global_morph") or "").strip()
            if species:
                if morph:
                    per_species_status[species][morph] += 1
                if morph in MORPHS:
                    per_species_morph[species][morph] += 1

    palette_count = [c for c in columns if c.startswith("palette_count_")]
    background_palette_count = [c for c in columns if c.startswith("background_palette_count_")]
    flower_fraction = [c for c in columns if c.startswith("flower_fraction_")]
    background_fraction = [c for c in columns if c.startswith("background_fraction_")]

    eligible = {}
    for sp, cnt in per_species_morph.items():
        n = sum(cnt.values())
        if n < 40:
            continue
        p = sorted((cnt.get(m, 0) / n for m in MORPHS), reverse=True)
        d = 1.0 - sum(x * x for x in p)
        eligible[sp] = {
            "n_four_state": n,
            "D": d,
            "second_fraction": p[1] if len(p) >= 2 else 0.0,
            "counts": dict(cnt),
        }

    dvals = [v["D"] for v in eligible.values()]
    seconds = [v["second_fraction"] for v in eligible.values()]
    mixed_rates = []
    for sp in eligible:
        st = per_species_status[sp]
        denom = sum(st.values())
        if denom:
            mixed_rates.append(st.get("mixed_uncertain", 0) / denom)

    return {
        "exists": True,
        "path": str(path.relative_to(ROOT)),
        "rows": n_rows,
        "columns": columns,
        "has_global_morph": "global_morph" in columns,
        "palette_count_columns": palette_count,
        "background_palette_count_columns": background_palette_count,
        "flower_fraction_columns": flower_fraction,
        "background_fraction_columns": background_fraction,
        "n_species_any_four_state": len(per_species_morph),
        "n_species_ge40_four_state": len(eligible),
        "D_min": min(dvals) if dvals else None,
        "D_max": max(dvals) if dvals else None,
        "second_ge_0_10_fraction": sum(x >= 0.10 for x in seconds) / len(seconds) if seconds else None,
        "second_ge_0_20_fraction": sum(x >= 0.20 for x in seconds) / len(seconds) if seconds else None,
        "mixed_rate_mean_among_ge40_species": sum(mixed_rates) / len(mixed_rates) if mixed_rates else None,
        "eligible_species": eligible,
    }


def fingerprint_distance(x: dict) -> float | None:
    required = (x.get("n_species_ge40_four_state"), x.get("D_max"), x.get("second_ge_0_10_fraction"), x.get("second_ge_0_20_fraction"))
    if any(v is None for v in required):
        return None
    n, dmax, s10, s20 = required
    # Scale count so a one-species discrepancy has comparable influence to ~0.0027 in fractions.
    return abs(n - 369) / 369 + abs(dmax - 0.71) + abs(s10 - 0.466) + abs(s20 - 0.266)


def main() -> None:
    report = {
        "audit": "polymorphism_directionality_step1_input_provenance",
        "date_jst": "2026-09-09",
        "inference_opened": False,
        "reported_fingerprint": {
            "n_species": 369,
            "D_range": [0.0, 0.71],
            "second_ge_0_10_fraction": 0.466,
            "second_ge_0_20_fraction": 0.266,
            "reported_correlations": {"D_SD_M1": 0.332, "D_SD_M2": 0.299, "D_mixed_rate": 0.537},
        },
        "tables": {},
        "repository_grep": {},
    }
    for name, path in TABLES.items():
        report["tables"][name] = summarize_table(path)
        report["tables"][name]["fingerprint_distance"] = fingerprint_distance(report["tables"][name])

    patterns = {
        "exact_reported_numbers": r"0\.332|0\.299|0\.537",
        "simpson_polymorphism_terms": r"[Ss]impson|polymorph|second[_ -]?morph|morph[_ -]?divers",
        "mode_sd_terms": r"SD\(M[123]\)|sd_m[123]|mode.*sd|species.*mode.*std",
    }
    report["repository_grep"] = {k: git_grep(v) for k, v in patterns.items()}

    # Pick a provenance candidate only from the fingerprint, never from inferential outcomes.
    candidates = []
    for name, x in report["tables"].items():
        dist = x.get("fingerprint_distance")
        if dist is not None:
            candidates.append((dist, name))
    candidates.sort()
    report["fingerprint_best_table"] = candidates[0][1] if candidates else None
    report["fingerprint_best_distance"] = candidates[0][0] if candidates else None

    json_path = OUT / "input_audit.json"
    json_path.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")

    lines = [
        "# Polymorphism directionality Step 1 — input provenance audit",
        "",
        "This audit was run before the Step 1 inferential statistic. `inference_opened = false`.",
        "",
        "## Fingerprint reconstruction",
        "",
        "| table | rows | species >=40 four-state | D max | second >=10% | second >=20% | 12 flower counts | 12 background counts |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for name in ("discovery", "reserve"):
        x = report["tables"][name]
        if not x.get("exists"):
            lines.append(f"| {name} | missing | | | | | | |")
            continue
        lines.append(
            f"| {name} | {x['rows']} | {x['n_species_ge40_four_state']} | "
            f"{x['D_max']:.6f} | {x['second_ge_0_10_fraction']:.4f} | {x['second_ge_0_20_fraction']:.4f} | "
            f"{len(x['palette_count_columns'])} | {len(x['background_palette_count_columns'])} |"
        )
    lines += [
        "",
        f"Fingerprint-best table: **{report['fingerprint_best_table']}** (distance {report['fingerprint_best_distance']}).",
        "",
        "## Exact-number / provenance grep",
        "",
    ]
    for key, hits in report["repository_grep"].items():
        lines.append(f"### {key}")
        if hits:
            lines.extend(f"- `{h}`" for h in hits[:80])
        else:
            lines.append("- no tracked-text hit")
        lines.append("")
    lines += [
        "## Decision boundary",
        "",
        "No Step 1 effect estimate or p-value is computed by this audit. Any source-cohort correction must be committed before the inferential runner is executed.",
        "",
    ]
    (OUT / "input_audit.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines[:30]))


if __name__ == "__main__":
    main()
