#!/usr/bin/env python3
"""Second refinement of the C/S strict audit.

Adds population-resolved coexistence evidence (e.g. `polymorphic population`) while
removing a broad S pattern that could confuse between-population experimental methods
with spatial segregation. Positive-source taxon extraction is rebuilt directly from
source text before GBIF validation.
"""
from __future__ import annotations

import re

import build_v22_coexistence_segregation_refined as base

COLORS = base.COLORS

# A population explicitly described as flower-colour polymorphic is a positive local
# coexistence statement. This remains distinct from generic species-level `polymorphic`.
C_ADDITIONS = [
    re.compile(r"(?is)\bpolymorphic\s+populations?\b"),
    re.compile(r"(?is)\bwithin[- ]population\s+(?:flower|floral)\s+colou?r\s+polymorph\w*\b"),
    re.compile(r"(?is)\b(?:flower|floral)\s+colou?r\s+polymorph\w*\b.{0,90}\b(?:found|observed|recorded)\s+in\s+(?:only\s+)?(?:one|a|the)\s+population\b"),
    re.compile(r"(?is)\bvariation\s+in\s+(?:flower|floral)\s+colou?r\b.{0,90}\bamong\s+individuals\s+within\s+(?:plant\s+)?populations?\b"),
    re.compile(r"(?is)\b(?:mixed|polymorphic)\s*,?\s*(?:colour|color)?[- ]?polymorphic\s+populations?\b"),
]
base.C_PATTERNS = C_ADDITIONS + base.C_PATTERNS

# Remove the broad inverse `among populations ... colour` pattern because it can match
# crossing methods or another response variable. Replace it with direct colour/morph
# spatial contrasts only.
S_ADDITIONS = [
    re.compile(rf"(?is)\beach\s+population\b.{{0,90}}\bpredominantly\s+{COLORS}\b.{{0,70}}\b(?:or|and)\b.{{0,70}}\b{COLORS}\b"),
    re.compile(rf"(?is)\bmost\s+populations\b.{{0,90}}\b{COLORS}\b.{{0,120}}\bsome\b.{{0,90}}\b{COLORS}\b"),
    re.compile(r"(?is)\bpopulations?\b.{0,70}\b(?:differ|differed|differing|vary|varied)\b.{0,70}\b(?:flower|floral)\s+colou?r\b"),
    re.compile(r"(?is)\b(?:flower|floral)\s+colou?r\b.{0,70}\b(?:differ|differed|differing|vary|varied)\b.{0,70}\b(?:among|between|across)\s+populations?\b"),
    re.compile(r"(?is)\b(?:proportion|proportions|frequency|frequencies)\s+of\b.{0,70}\b(?:morphs?|colou?r\s+morphs?)\b.{0,100}\b(?:significantly\s+)?different\s+between\s+populations?\b"),
    re.compile(r"(?is)\bmonomorphic[- ]?[a-z]*\s+populations?\b.{0,150}\bpolymorphic\s+populations?\b"),
    re.compile(r"(?is)\bpolymorphic\s+populations?\b.{0,150}\bmonomorphic[- ]?[a-z]*\s+populations?\b"),
    re.compile(r"(?is)\bpolymorphic\s+populations?\s+(?:are|were)\s+restricted\s+to\b.{0,150}\bmonomorphic\s+populations?\b"),
    re.compile(r"(?is)\b(?:flower|floral)?\s*colou?r\s+morphs?\b.{0,90}\b(?:distinct|different)\s+geographic\s+regions?\b"),
]
base.S_PATTERNS = S_ADDITIONS + [pattern for index, pattern in enumerate(base.S_PATTERNS) if index != 1]

# Articles labelled as perspectives are not primary evidence even when OpenAlex calls
# them `article`. Recompile explicitly instead of concatenating an inline-flag pattern.
base.NONPRIMARY_TEXT_RE = re.compile(
    r"\b(?:systematic\s+review|meta[- ]analysis|review\s+article|we\s+review|we\s+also\s+review|"
    r"here\s+we\s+summari[sz]e|this\s+review|review\s+the\s+incidence|PERSPECTIVE)\b",
    re.IGNORECASE | re.DOTALL,
)

# Fresh source-text binomials are used for C/S-positive records because the archived
# navigation extractor can be dominated by false two-word phrases before the focal
# species name appears. GBIF still decides whether candidates are accepted plant taxa.
BINOMIAL_RE = re.compile(r"\b([A-Z][a-z][A-Za-z-]{1,})\s+([a-z][a-z-]{2,})\b")


def _extract_binomials(text: str) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for genus, epithet in BINOMIAL_RE.findall(text):
        name = f"{genus} {epithet}"
        if name not in seen:
            seen.add(name)
            out.append(name)
    return out


def ranked_candidates_v2(hidden: dict[str, str], title: str, abstract: str, historical_taxon: str) -> list[str]:
    if historical_taxon:
        return [historical_taxon]

    title_candidates = _extract_binomials(title)
    if title_candidates:
        return title_candidates[:8]

    abstract_candidates = _extract_binomials(abstract[:4000])
    if abstract_candidates:
        return abstract_candidates[:12]

    # Navigation hints are a final fallback only; they never override source text.
    fallback = [base.clean(x) for x in base.clean(hidden.get("detected_binomial_strings")).split(";") if base.clean(x)]
    return fallback[:8]


base.ranked_candidates = ranked_candidates_v2

if __name__ == "__main__":
    base.main()
