#!/usr/bin/env python3
"""Apply Journal of Biogeography figure-format compliance fixes, fail closed.

Layout/data/statistics are unchanged. The patch only raises the renderer typography
floor for 168-mm reduction, changes panel tags to lowercase parenthesized labels,
and synchronizes the generated captions to the same panel convention.
"""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PATH = ROOT / "scripts/analysis/render_polymorphism_paper_v0_1_figures_20260910.py"
MIN_SOURCE_FONT_PT = 9.5

RC_REPLACEMENTS = {
    '"font.size": 8.5,': '"font.size": 9.5,',
    '"axes.titlesize": 10,': '"axes.titlesize": 10.5,',
    '"axes.labelsize": 9,': '"axes.labelsize": 9.5,',
    '"xtick.labelsize": 7.5,': '"xtick.labelsize": 9.5,',
    '"ytick.labelsize": 7.5,': '"ytick.labelsize": 9.5,',
    '"legend.fontsize": 7.5,': '"legend.fontsize": 9.5,',
}

EXPECTED_SMALL_FONT_COUNTS = {
    6.5: 1,
    6.6: 3,
    6.7: 1,
    7.0: 5,
    7.2: 1,
    7.5: 4,
}

PANEL_OLD = 'ax.text(-0.12, 1.06, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")'
PANEL_NEW = 'ax.text(-0.12, 1.06, f"({letter.lower()})", transform=ax.transAxes, fontsize=11.5, fontweight="bold", va="top")'

CAPTION_COUNTS = {
    "(A–B)": 1,
    "(A)": 3,
    "(B)": 3,
    "(C)": 3,
    "(D)": 3,
}
CAPTION_REPLACEMENTS = {
    "(A–B)": "(a–b)",
    "(A)": "(a)",
    "(B)": "(b)",
    "(C)": "(c)",
    "(D)": "(d)",
}


def fail_closed_replace(text: str, old: str, new: str, *, expected: int = 1) -> str:
    count = text.count(old)
    if count != expected:
        raise RuntimeError(f"JBI figure patch expected {expected} occurrence(s), found {count}: {old!r}")
    return text.replace(old, new)


def main() -> int:
    text = PATH.read_text(encoding="utf-8")
    original = text

    for old, new in RC_REPLACEMENTS.items():
        text = fail_closed_replace(text, old, new)

    text = fail_closed_replace(text, PANEL_OLD, PANEL_NEW)

    values = [float(v) for v in re.findall(r"fontsize=(\d+(?:\.\d+)?)", text)]
    observed_small = {value: values.count(value) for value in sorted(set(values)) if value < MIN_SOURCE_FONT_PT}
    if observed_small != EXPECTED_SMALL_FONT_COUNTS:
        raise RuntimeError(
            "unexpected sub-floor explicit font inventory before compliance patch: "
            f"observed={observed_small}, expected={EXPECTED_SMALL_FONT_COUNTS}"
        )

    def raise_font(match: re.Match[str]) -> str:
        value = float(match.group(1))
        if value < MIN_SOURCE_FONT_PT:
            return f"fontsize={MIN_SOURCE_FONT_PT:g}"
        return match.group(0)

    text = re.sub(r"fontsize=(\d+(?:\.\d+)?)", raise_font, text)

    for old, expected in CAPTION_COUNTS.items():
        count = text.count(old)
        if count != expected:
            raise RuntimeError(f"caption panel token drift for {old}: {count} != {expected}")
    for old, new in CAPTION_REPLACEMENTS.items():
        text = text.replace(old, new)

    remaining_explicit = [float(v) for v in re.findall(r"fontsize=(\d+(?:\.\d+)?)", text)]
    below = sorted({v for v in remaining_explicit if v < MIN_SOURCE_FONT_PT})
    if below:
        raise RuntimeError(f"explicit figure font remains below {MIN_SOURCE_FONT_PT:g} pt: {below}")
    if any(token in text for token in CAPTION_COUNTS):
        raise RuntimeError("uppercase panel labels remain in generated caption text")
    if 'f"({letter.lower()})"' not in text:
        raise RuntimeError("lowercase parenthesized panel renderer not installed")

    if text == original:
        raise RuntimeError("JBI figure compliance patch made no change")
    PATH.write_text(text, encoding="utf-8", newline="\n")
    print(
        f"patched {PATH.relative_to(ROOT)}: typography floor={MIN_SOURCE_FONT_PT:g} pt; "
        "panel tags=(a),(b),...; generated captions synchronized"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
