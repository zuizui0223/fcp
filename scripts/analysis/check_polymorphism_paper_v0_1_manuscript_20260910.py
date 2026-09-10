#!/usr/bin/env python3
"""Fail if the FCP polymorphism v0.1 manuscript drifts from the frozen paper ledger."""
from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "paper" / "polymorphism_v0_1"
MANUSCRIPT = PAPER / "MANUSCRIPT.md"
NUMBERS = PAPER / "paper_numbers.json"
CAPTIONS = PAPER / "FIGURE_CAPTIONS.md"
QC = PAPER / "FIGURE_QC.md"


def require(text: str, token: str) -> None:
    if token not in text:
        raise RuntimeError(f"required manuscript token missing: {token}")


def forbid(text: str, token: str) -> None:
    if token.lower() in text.lower():
        raise RuntimeError(f"forbidden overclaim present: {token}")


def main() -> None:
    for p in [MANUSCRIPT, NUMBERS, CAPTIONS, QC]:
        if not p.exists():
            raise RuntimeError(f"missing paper surface: {p.relative_to(ROOT)}")

    m = MANUSCRIPT.read_text(encoding="utf-8")
    n = json.loads(NUMBERS.read_text(encoding="utf-8"))

    first = m.splitlines()[0].removeprefix("# ").strip()
    if first != n["title"]:
        raise RuntimeError(f"title drift: {first!r} != {n['title']!r}")

    # Section surface contract.
    for heading in [
        "## Abstract",
        "## Introduction",
        "## Materials and Methods",
        "## Results",
        "## Discussion",
        "## Data availability and reproducibility",
        "## References",
    ]:
        require(m, heading)

    # Key numbers are intentionally checked in rounded paper-facing form.
    required_numbers = [
        "369 species",
        "363 species",
        "0.707645",
        "46.61%",
        "26.56%",
        "0.089213",
        "0.101601",
        "0.099288",
        "0.116241",
        "0.097078",
        "0.125286",
        "0.116299",
        "0.132785",
        "0.202963",
        "0.153654",
        "0.530632",
        "0.06165",
    ]
    for x in required_numbers:
        require(m, x)

    # Inferential boundaries must remain visible.
    for phrase in [
        "species-disjoint reserve",
        "not random samples of all angiosperms",
        "not a formal phylogenetic comparative test",
        "not an estimate of biological range size",
        "uniform endpoint stress tests",
        "does not identify selection",
    ]:
        require(m, phrase)

    for phrase in [
        "polymorphism is adaptive",
        "range size causes polymorphism",
        "phylogenetic signal is significant",
        "pollinators drive global polymorphism",
        "universal climate threshold",
        "ambiguous observations are true morphs",
    ]:
        forbid(m, phrase)

    # Main-text word count: Introduction through end of Discussion, excluding headings.
    start = m.index("## Introduction")
    end = m.index("## Data availability and reproducibility")
    main_text = m[start:end]
    main_text = re.sub(r"^#+.*$", " ", main_text, flags=re.MULTILINE)
    main_text = re.sub(r"`[^`]+`", " ", main_text)
    words = re.findall(r"\b[\w–-]+\b", main_text)
    if not (2500 <= len(words) <= 7500):
        raise RuntimeError(f"unexpected main-text word count: {len(words)}")

    # Figure references and core literature surface.
    for fig in ["Figure 1", "Figure 2", "Figure 3", "Figure 4"]:
        require(m, fig)
    for ref in ["Delph & Kelly, 2014", "Narbona et al., 2018", "Trunschke et al., 2021", "Dalrymple et al., 2020"]:
        require(m, ref)

    receipt = {
        "paper_version": n["paper_version"],
        "title": first,
        "main_text_word_count": len(words),
        "required_number_tokens_checked": len(required_numbers),
        "forbidden_overclaims_checked": 6,
        "status": "PASS",
    }
    out = PAPER / "MANUSCRIPT_CHECK.json"
    out.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
