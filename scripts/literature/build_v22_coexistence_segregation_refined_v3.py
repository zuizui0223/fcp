#!/usr/bin/env python3
"""Third refinement: combine title+abstract binomial candidates before GBIF validation."""
from __future__ import annotations

import re

import build_v22_coexistence_segregation_refined_v2 as rules

base = rules.base
BINOMIAL_RE = re.compile(r"\b([A-Z][a-z][A-Za-z-]{1,})\s+([a-z][a-z-]{2,})\b")


def extract_candidates(text: str) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for genus, epithet in BINOMIAL_RE.findall(text):
        variants = [f"{genus} {epithet}"]
        # OpenAlex occasionally duplicates the terminal `is` in indexed abstracts
        # (e.g. arvensis -> arvensisis). Keep the literal string and add a conservative
        # correction as a second GBIF candidate rather than silently replacing it.
        if epithet.endswith("isis") and len(epithet) > 6:
            variants.append(f"{genus} {epithet[:-2]}")
        for name in variants:
            if name not in seen:
                seen.add(name)
                out.append(name)
    return out


def ranked_candidates_v3(hidden: dict[str, str], title: str, abstract: str, historical_taxon: str) -> list[str]:
    if historical_taxon:
        return [historical_taxon]

    combined: list[str] = []
    seen: set[str] = set()
    for name in extract_candidates(title) + extract_candidates(abstract[:5000]):
        if name not in seen:
            seen.add(name)
            combined.append(name)
    if combined:
        return combined[:20]

    fallback = [base.clean(x) for x in base.clean(hidden.get("detected_binomial_strings")).split(";") if base.clean(x)]
    return fallback[:8]


base.ranked_candidates = ranked_candidates_v3

if __name__ == "__main__":
    base.main()
