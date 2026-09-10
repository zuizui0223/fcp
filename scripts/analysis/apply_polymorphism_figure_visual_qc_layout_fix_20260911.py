#!/usr/bin/env python3
"""Apply second-pass layout-only fixes discovered by journal-size visual QC, fail closed."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PATH = ROOT / "scripts/analysis/render_polymorphism_paper_v0_1_figures_20260910.py"

REPLACEMENTS = [
    (
        '(0.37, 0.62, 0.26, 0.22, "image measurement\\n+ status gate"),',
        '(0.37, 0.62, 0.26, 0.22, "measurement\\n+ status gate"),',
    ),
    (
        '(0.71, 0.62, 0.26, 0.22, "4 admitted\\ncolour states"),',
        '(0.71, 0.62, 0.26, 0.22, "4 colour\\nstates"),',
    ),
    (
        'ax.text(x + w/2, y + h/2, label, transform=ax.transAxes, ha="center", va="center", fontsize=7.5)',
        'ax.text(x + w/2, y + h/2, label, transform=ax.transAxes, ha="center", va="center", fontsize=6.5)',
    ),
    (
        'labels = ["Disc. <10%", "Disc. ≥10%", "Res. <10%", "Res. ≥10%"]',
        'labels = ["Disc.\\n<10%", "Disc.\\n≥10%", "Reserve\\n<10%", "Reserve\\n≥10%"]',
    ),
    (
        'ax.legend(frameon=False, ncol=3, loc="upper left")',
        'ax.legend(frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 0.04))',
    ),
    (
        'ax.text(0.02, 0.03, "numbers above points are geometry-preserving p", transform=ax.transAxes, fontsize=6.5, color=MID)',
        '# Geometry-preserving P-value provenance is stated in the figure caption.',
    ),
    (
        'ax.legend(frameon=False, loc="upper left")',
        'ax.legend(frameon=False, loc="lower right")',
    ),
    (
        'fig = plt.figure(figsize=(7.2, 3.3))',
        'fig = plt.figure(figsize=(7.2, 3.8))',
    ),
    (
        'gs = fig.add_gridspec(1, 2, wspace=0.42)',
        'gs = fig.add_gridspec(1, 2, wspace=0.42, top=0.78)',
    ),
    (
        'fig.suptitle("Figure 4 | Genus-level structure is a secondary, qualified result", fontsize=11, y=1.01)',
        'fig.suptitle("Figure 4 | Genus-level structure is a secondary, qualified result", fontsize=11, y=0.98)',
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
    print(f"patched {PATH.relative_to(ROOT)} with {len(REPLACEMENTS)} second-pass layout-only replacements")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
