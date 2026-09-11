#!/usr/bin/env python3
"""Apply final label-placement-only fixes after 168-mm visual inspection, fail closed.

This patch changes no data, statistics, axis limits, panel geometry, or typography floor.
It only shortens two Figure 2c tranche tick labels and separates the six deliberately
labelled genera in Figure 4b using fixed point offsets.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PATH = ROOT / "scripts/analysis/render_polymorphism_paper_v0_1_figures_20260910.py"

OLD_TICKS = 'labels = ["Disc.\\n<10%", "Disc.\\n≥10%", "Reserve\\n<10%", "Reserve\\n≥10%"]'
NEW_TICKS = 'labels = ["Disc.\\n<10%", "Disc.\\n≥10%", "Res.\\n<10%", "Res.\\n≥10%"]'

OLD_LABEL_LOOP = '''    for i in sorted(label_idx):
        row = shared.iloc[i]
        ax.text(row.discovery_mean_D + 0.007, row.reserve_mean_D + 0.007, row.genus, fontsize=9.5)
'''

NEW_LABEL_LOOP = '''    label_offsets = {
        "Allium": (5, 4, "left"),
        "Calochortus": (5, 5, "left"),
        "Erythranthe": (5, 4, "left"),
        "Ipomoea": (6, 8, "left"),
        "Rosa": (5, 8, "left"),
        "Salvia": (-6, 5, "right"),
    }
    for i in sorted(label_idx):
        row = shared.iloc[i]
        dx, dy, ha = label_offsets.get(row.genus, (5, 5, "left"))
        ax.annotate(
            row.genus,
            (row.discovery_mean_D, row.reserve_mean_D),
            xytext=(dx, dy),
            textcoords="offset points",
            ha=ha,
            va="center",
            fontsize=9.5,
        )
'''


def replace_once(text: str, old: str, new: str, label: str) -> str:
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"final-label QC expected exactly one {label} target, found {count}")
    return text.replace(old, new, 1)


def main() -> int:
    text = PATH.read_text(encoding="utf-8")
    original = text
    text = replace_once(text, OLD_TICKS, NEW_TICKS, "Figure 2c tick-label")
    text = replace_once(text, OLD_LABEL_LOOP, NEW_LABEL_LOOP, "Figure 4b genus-label loop")
    for token in [
        '"font.size": 9.5,',
        '"xtick.labelsize": 9.5,',
        '"ytick.labelsize": 9.5,',
        '"legend.fontsize": 9.5,',
        'f"({letter.lower()})"',
        '(0.16, 0.78, 0.68, 0.14, "100 photos per species")',
    ]:
        if token not in text:
            raise RuntimeError(f"required final figure contract token missing: {token}")
    if 'ax.text(x, y + 0.006' in text:
        raise RuntimeError("Figure 3c pointwise p-value labels unexpectedly reintroduced")
    if text == original:
        raise RuntimeError("final-label QC patch made no change")
    PATH.write_text(text, encoding="utf-8", newline="\n")
    print("patched Figure 2c tranche abbreviations and Figure 4b genus-label offsets")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
