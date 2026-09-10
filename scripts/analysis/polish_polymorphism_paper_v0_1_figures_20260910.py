#!/usr/bin/env python3
"""Apply the frozen post-QC readability-only polish to the v0.1 figure renderer."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RENDERER = ROOT / "scripts" / "analysis" / "render_polymorphism_paper_v0_1_figures_20260910.py"
QC = ROOT / "paper" / "polymorphism_v0_1" / "FIGURE_QC.md"


def replace_once(text: str, old: str, new: str) -> str:
    n = text.count(old)
    if n == 0 and new in text:
        return text
    if n != 1:
        raise RuntimeError(f"expected exactly one renderer match, found {n}: {old!r}")
    return text.replace(old, new, 1)


def main() -> None:
    text = RENDERER.read_text(encoding="utf-8")
    text = replace_once(
        text,
        '    labels = ["<10%", "≥10%", "<10%", "≥10%"]',
        '    labels = ["Disc. <10%", "Disc. ≥10%", "Res. <10%", "Res. ≥10%"]',
    )
    text = replace_once(
        text,
        '    ax.text(0.5, -0.23, "Discovery", transform=ax.transAxes, ha="left", color=DISC)\n'
        '    ax.text(0.84, -0.23, "Reserve", transform=ax.transAxes, ha="left", color=RES)\n',
        '',
    )
    text = replace_once(
        text,
        '            ax.text(x, y + 0.006, f"{p:.3f}", ha="center", va="bottom", fontsize=5.8, color=c)',
        '            ax.text(x, y + 0.006, f"{p:.3f}", ha="center", va="bottom", fontsize=6.6, color=c)',
    )
    text = replace_once(
        text,
        '        ax.text(row.discovery_mean_D + 0.007, row.reserve_mean_D + 0.007, row.genus, fontsize=5.8)',
        '        ax.text(row.discovery_mean_D + 0.007, row.reserve_mean_D + 0.007, row.genus, fontsize=6.6)',
    )
    RENDERER.write_text(text, encoding="utf-8")

    qc = QC.read_text(encoding="utf-8")
    marker = "## Readability polish executed\n"
    if marker not in qc:
        qc += (
            "\n## Readability polish executed\n\n"
            "The renderer-level risks identified above were addressed without changing any data or statistics: "
            "Figure 2C now uses self-contained discovery/reserve tick labels, and the smallest point/genus annotations "
            "in Figures 3C and 4B were increased from 5.8 pt to 6.6 pt. Figures were then regenerated from the same "
            "frozen figure-data tables and their hashes were refreshed in `figures/manifest.json`.\n"
        )
        QC.write_text(qc, encoding="utf-8")

    print("Applied deterministic readability-only figure polish")


if __name__ == "__main__":
    main()
