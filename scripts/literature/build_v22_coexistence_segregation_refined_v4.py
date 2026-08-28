#!/usr/bin/env python3
"""Fourth refinement of the coexistence (C) / segregation (S) audit.

Adds generic population-resolved evidence patterns that were still missed after v3:
- explicit within-population polymorphism / mixed populations;
- natural colour forms described as intermixed in populations;
- population-level morph-frequency ranges where at least two morphs have positive
  lower bounds, which directly implies local coexistence;
- explicit restriction of flower-colour polymorphism to one/some populations;
- monomorphic-vs-polymorphic population contrasts and altitudinal/spatial changes in
  polymorphism or morph frequencies.

These are positive-evidence rules. They do not infer C=0 or S=0 from silence.
"""
from __future__ import annotations

import re

import build_v22_coexistence_segregation_refined_v3 as v3

base = v3.base
COLORS = base.COLORS

# Positive percentage/range whose lower bound is >0. Used only in a population-level
# morph-frequency context. Two such morphs imply at least two colour morphs coexist in
# the described population set.
POS_PCT = r"(?:[1-9]\d*(?:\.\d+)?)(?:\s*[-–]\s*\d+(?:\.\d+)?)?\s*%"

C_POPULATION_EVIDENCE = [
    # Direct terminology: source explicitly says polymorphism is within populations.
    re.compile(r"(?is)\bwithin[- ]population\s+(?:flower[- ]?colou?r\s+)?polymorph\w*\b"),
    re.compile(r"(?is)\b(?:some|several|many|the|these|those)?\s*populations?\s+(?:show(?:ing|ed|s)?|exhibit(?:ing|ed|s)?|contain(?:ing|ed|s)?|have|had)\b.{0,90}\b(?:flower|floral)?\s*colou?r\s+polymorph\w*\b"),
    re.compile(r"(?is)\b(?:mixed[- ]colou?r|mixed\s+colou?r)\s+populations?\b"),
    re.compile(r"(?is)\bin\s+(?:a|the|one|some|several|many)\s+mixed[- ]?colou?r\s+populations?\b"),
    re.compile(r"(?is)\b(?:colour|color)\s+(?:forms?|morphs?)\b.{0,80}\b(?:seen|found|observed|occurring?)\s+intermixed\b.{0,90}\bpopulations?\b"),
    re.compile(r"(?is)\bintermixed\b.{0,80}\b(?:colour|color)\s+(?:forms?|morphs?)\b.{0,90}\bpopulations?\b"),
    # A flower-colour polymorphism explicitly located in a named/numbered population.
    re.compile(r"(?is)\b(?:flower|floral)\s+colou?r\s+polymorph\w*\b.{0,90}\b(?:found|observed|recorded|present|occurs?|occurring)\b.{0,45}\bin\s+(?:only\s+)?(?:one|a|the|some|several)\s+populations?\b"),
    # Population-level morph-frequency ranges with >=2 morphs whose lower frequency is >0.
    # This captures tables/results such as violet 20-43%, pink 38-69%, white 11-30%.
    re.compile(
        rf"(?is)\bpopulations?\b.{{0,520}}\bmorph\s+frequenc(?:y|ies)\b.{{0,220}}"
        rf"\b{COLORS}\b\s*{POS_PCT}.{{0,130}}\b{COLORS}\b\s*{POS_PCT}"
    ),
    re.compile(
        rf"(?is)\bmorph\s+frequenc(?:y|ies)\b.{{0,220}}"
        rf"\b{COLORS}\b\s*{POS_PCT}.{{0,130}}\b{COLORS}\b\s*{POS_PCT}"
        rf".{{0,420}}\bpopulations?\b"
    ),
]

S_POPULATION_EVIDENCE = [
    # Flower-colour polymorphism explicitly restricted to one/some populations.
    re.compile(r"(?is)\b(?:flower|floral)\s+colou?r\s+polymorph\w*\b.{0,100}\b(?:found|observed|recorded|present|restricted)\b.{0,55}\b(?:only\s+)?(?:one|a|some|several)\s+populations?\b"),
    re.compile(r"(?is)\b(?:only|some|several)\s+populations?\b.{0,100}\b(?:flower|floral)\s+colou?r\s+polymorph\w*\b"),
    # Direct monomorphic-versus-polymorphic population contrasts.
    re.compile(r"(?is)\bmonomorphic\b.{0,70}\bpopulations?\b.{0,240}\bpolymorph\w*\b.{0,70}\bpopulations?\b"),
    re.compile(r"(?is)\bpolymorph\w*\b.{0,70}\bpopulations?\b.{0,240}\bmonomorphic\b.{0,70}\bpopulations?\b"),
    re.compile(r"(?is)\bpopulations?\b.{0,90}\b(?:uniformly\s+colou?red|monomorphic)\b.{0,180}\b(?:other|some|several|four|three|two|one|\d+)\b.{0,60}\b(?:mixed|polymorph\w*)\b"),
    # Spatial/altitudinal distribution of the polymorphism itself.
    re.compile(r"(?is)\bdistribution\s+of\s+(?:the\s+)?(?:flower[- ]?colou?r\s+)?polymorph\w*\b.{0,120}\b(?:altitude|elevation|geographic|spatial|region|latitude)\b"),
    re.compile(r"(?is)\b(?:degree|extent)\s+of\s+(?:flower[- ]?colou?r\s+)?polymorph\w*\b.{0,130}\b(?:altitude|elevation|geographic|spatial|region|latitude)\b"),
    re.compile(r"(?is)\bmorph\s+frequenc(?:y|ies)\s+among\b.{0,80}\bpopulations?\b.{0,140}\b(?:altitude|elevation|geographic|spatial|region|latitude)\b"),
    # Below/above or low/high geographic units changing from monomorphic to polymorphic.
    re.compile(r"(?is)\bpopulations?\s+(?:below|above|at|from)\b.{0,130}\bmonomorphic\b.{0,220}\b(?:above|below|higher|lower)\b.{0,140}\bpolymorph\w*\b"),
    # Explicit regional frequency contrasts in the same study.
    re.compile(
        rf"(?is)\b(?:in|within)\s+[A-Z][A-Za-z-]+\b.{{0,180}}\b{COLORS}\b\s*{POS_PCT}.{{0,260}}"
        rf"\b(?:in|within)\s+[A-Z][A-Za-z-]+\b.{{0,180}}\b{COLORS}\b\s*{POS_PCT}"
    ),
]

# Prepend high-specificity population evidence to the existing v3 rule set.
base.C_PATTERNS = C_POPULATION_EVIDENCE + base.C_PATTERNS
base.S_PATTERNS = S_POPULATION_EVIDENCE + base.S_PATTERNS

if __name__ == "__main__":
    base.main()
