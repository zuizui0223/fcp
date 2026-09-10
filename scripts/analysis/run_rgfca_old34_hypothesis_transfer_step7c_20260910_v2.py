#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
SRC = HERE / "run_rgfca_old34_hypothesis_transfer_step7c_20260910.py"
spec = importlib.util.spec_from_file_location("step7c_base", SRC)
base = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(base)


def fmt(v, spec: str) -> str:
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return "not estimable"
    return format(v, spec)


def main() -> None:
    states = pd.read_csv(base.STATES)
    states["C_star"] = base.bseries(states["C_star"])
    states["S_star"] = base.bseries(states["S_star"])
    dm, ds = base.run_one(base.DISC, "discovery", states, base.SEED + 1)
    rm, rs = base.run_one(base.RES, "reserve", states, base.SEED + 2)
    metrics = pd.concat([dm, rm], ignore_index=True)
    metrics.to_csv(base.OUT / "species_transfer_metrics.csv", index=False)

    recurrent = {
        "H2": bool(ds["H2_primary_pass"] and rs["H2_primary_pass"]),
        "H9": bool(ds["H9_primary_pass"] and rs["H9_primary_pass"]),
    }
    result = {
        "analysis": "rgfca_old34_hypothesis_transfer_step7c_photo_only",
        "protocol": "docs/RGFCA_OLD34_HYPOTHESIS_TRANSFER_STEP7C_PROTOCOL_20260910.md",
        "discovery": ds,
        "reserve": rs,
        "two_tranche_recurrent": recurrent,
        "environmental_family_status": "NOT_OPENED_REQUIRES_SEPARATE_SOURCE_FREEZE",
        "claim_boundary": "H2 is sampled-coordinate fragmentation, not biological range fragmentation. H9 is photo-date partitioning, not direct flowering phenology. C*/S* are photo-derived states, so this is a transferred test rather than a confirmatory replication of the historical 34-species response.",
    }
    (base.OUT / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    lines = ["# RGFCA Step 7C — historical hypothesis transfer, photo-only family", ""]
    for tr, s in [("Discovery", ds), ("Reserve", rs)]:
        h2 = s["H2_fragmentation"]
        h2s = s["H2_span_adjusted_robustness"]
        h9 = s["H9_photo_date_partitioning"]
        lines += [
            f"## {tr}", "",
            f"- eligible C*/S* species: **{s['species']}**",
            f"- H2 fragmentation S*=1 minus S*=0 median contrast: **{fmt(h2['delta'], '.6f')}**, raw p **{fmt(h2['p'], '.6g')}**, Holm p **{fmt(h2['holm_p'], '.6g')}**",
            f"- H2 span-adjusted robustness contrast: **{fmt(h2s['delta'], '.6f')}**, p **{fmt(h2s['p'], '.6g')}**",
            f"- H9 eligible species: **{h9['eligible_species']}**",
            f"- H9 C*=1 minus C*=0 phenology-delta contrast: **{fmt(h9['delta'], '.3f')} days**, raw p **{fmt(h9['p'], '.6g')}**, Holm p **{fmt(h9['holm_p'], '.6g')}**",
            f"- H9 rho(phenology proxy, D): **{fmt(h9['rho_with_D'], '.4f')}**; rho(proxy, S*): **{fmt(h9['rho_with_S_star'], '.4f')}**", "",
        ]
    lines += [
        "## Two-tranche recurrence", "",
        f"- H2: **{recurrent['H2']}**",
        f"- H9: **{recurrent['H9']}**", "",
        "H0/H1a/H1b/H3a/H3b/H3c/H8 remain unopened until the RGFCA environmental source is frozen. H2 is sampled-coordinate fragmentation; H9 is photo-date partitioning rather than direct flowering phenology.",
    ]
    (base.OUT / "RESULT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
