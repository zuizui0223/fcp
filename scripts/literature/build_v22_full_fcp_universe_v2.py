#!/usr/bin/env python3
"""Refine the all-species FCP universe with hierarchical focal-taxon attribution.

The v1 universe correctly made membership independent of C/S positivity, but inherited a
taxon candidate function that pooled title and abstract binomials. A source mentioning
one focal plant in the title and comparison plants in the abstract could therefore be
marked taxonomically ambiguous. This wrapper uses source-structure priority instead:

1. historical exact-source taxon rescue (taxonomy only; no old state labels),
2. source-text binomials from the title,
3. source-text binomials from the abstract,
4. archived navigation hints as a final fallback.

Only one tier is sent to GBIF. Thus lower-priority comparison taxa cannot override a
clear higher-priority focal taxon. No C/S, climate or historical state information is
used to choose the tier.
"""
from __future__ import annotations

import re

import build_v22_full_fcp_universe as v1

base = v1.base

# More tolerant than the old word-boundary extractor: OpenAlex occasionally concatenates
# preceding words with a genus (e.g. dimorphicBuddleja delavayi). Starting at the capital
# letter remains conservative because GBIF validates every candidate as a plant species.
BINOMIAL_RE = re.compile(r"([A-Z][a-z][A-Za-z-]{1,})\s+([a-z][a-z-]{2,})")


def extract_candidates(text: str) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for genus, epithet in BINOMIAL_RE.findall(text or ""):
        variants = [f"{genus} {epithet}"]
        if epithet.endswith("isis") and len(epithet) > 6:
            variants.append(f"{genus} {epithet[:-2]}")
        for name in variants:
            if name not in seen:
                seen.add(name)
                out.append(name)
    return out


def ranked_candidates_hierarchical(hidden: dict[str, str], title: str, abstract: str, historical_taxon: str) -> list[str]:
    if historical_taxon:
        return [historical_taxon]

    title_candidates = extract_candidates(title)
    if title_candidates:
        return title_candidates[:12]

    abstract_candidates = extract_candidates(abstract[:5000])
    if abstract_candidates:
        return abstract_candidates[:20]

    fallback = [base.clean(x) for x in base.clean(hidden.get("detected_binomial_strings")).split(";") if base.clean(x)]
    return fallback[:10]


base.ranked_candidates = ranked_candidates_hierarchical

if __name__ == "__main__":
    v1.main()
