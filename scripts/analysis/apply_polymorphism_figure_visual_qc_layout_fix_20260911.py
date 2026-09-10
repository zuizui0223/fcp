#!/usr/bin/env python3
"""Apply layout-only fixes discovered by journal-size visual QC, fail closed."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PATH = ROOT / "scripts/analysis/render_polymorphism_paper_v0_1_figures_20260910.py"

REPLACEMENTS = [
    (
        '(0.03, 0.62, 0.26, 0.22, "100 fixed\\nphotos / species"),',
        '(0.03, 0.62, 0.26, 0.22, "100 photos\\nper species"),',
    ),
    (
        '(0.37, 0.62, 0.26, 0.22, "measurement +\\nstatus partition"),',
        '(0.37, 0.62, 0.26, 0.22, "image measurement\\n+ status gate"),',
    ),
    (
        'ax.text(x + w/2, y + h/2, label, transform=ax.transAxes, ha="center", va="center")',
        'ax.text(x + w/2, y + h/2, label, transform=ax.transAxes, ha="center", va="center", fontsize=7.5)',
    ),
    (
        'ax.set(xlabel="D", ylabel="Second-largest morph fraction", title="Threshold summaries preserve a continuous gradient")',
        'ax.set(xlabel="D", ylabel="Second-largest morph fraction", title="Threshold summaries preserve\\na continuous gradient")',
    ),
    (
        'ax.set(xlabel="Species in admitted discovery frame (%)", xlim=(0, 55), title="Non-trivial secondary morphs are common in-frame")',
        'ax.set(xlabel="Species in admitted discovery frame (%)", xlim=(0, 55), title="Non-trivial secondary morphs\\nare common in-frame")',
    ),
    (
        'ax.text(0, 1.35, "Not an angiosperm-wide prevalence estimate", fontsize=7, color=MID)',
        'ax.text(0.98, 0.04, "not angiosperm-wide", transform=ax.transAxes, ha="right", va="bottom", fontsize=6.6, color=MID)',
    ),
    (
        'ax.set(ylabel="Mean species spatial ρ ± bootstrap 95% CI", title="Polymorphic subset is more spatially organized")',
        'ax.set(ylabel="Mean species spatial ρ ± bootstrap 95% CI", title="Polymorphic subset has stronger\\nspatial organization")',
    ),
    (
        'ax.set(xlabel="Equal-species mean spatial ρ (≥10% subset)", title="Reserve replication survives nuisance controls")',
        'ax.set(xlabel="Equal-species mean spatial ρ (≥10% subset)", title="Reserve replication survives\\nnuisance controls")',
    ),
    (
        'ax.set(ylabel="Partial Spearman ρ\\n| sampled span + technical failure", title="Spatial signal survives both ambiguity endpoints")',
        'ax.set(ylabel="Partial Spearman ρ\\n| sampled span + technical failure", title="Spatial signal survives\\nboth ambiguity endpoints")',
    ),
    (
        'ax.set(xlabel="Reserve species ranked by observed D", ylabel="D", title="Observed ordering is stable across broad completion bounds")',
        'ax.set(xlabel="Reserve species ranked by observed D", ylabel="D", title="Observed ordering remains stable\\nacross completion bounds")',
    ),
    (
        'ax.set(xlabel="Genus clustering gain", title="Taxonomic clustering is replicated but ambiguity-sensitive")',
        'ax.set(xlabel="Genus clustering gain", title="Genus clustering replicates,\\nbut is ambiguity-sensitive")',
    ),
    (
        'ax.set(xlim=lim, ylim=lim, xlabel="Discovery genus mean D", ylabel="Reserve genus mean D", title="23 repeated genera show cross-tranche concordance")',
        'ax.set(xlim=lim, ylim=lim, xlabel="Discovery genus mean D", ylabel="Reserve genus mean D", title="Repeated genera show\\ncross-tranche concordance")',
    ),
]


def main() -> int:
    text = PATH.read_text(encoding="utf-8")
    original = text
    for old, new in REPLACEMENTS:
        count = text.count(old)
        if count != 1:
            raise RuntimeError(f"visual-QC patch expected exactly one target, found {count}: {old[:100]!r}")
        text = text.replace(old, new, 1)
    if text == original:
        raise RuntimeError("visual-QC patch made no change")
    PATH.write_text(text, encoding="utf-8", newline="\n")
    print(f"patched {PATH.relative_to(ROOT)} with {len(REPLACEMENTS)} layout-only replacements")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
