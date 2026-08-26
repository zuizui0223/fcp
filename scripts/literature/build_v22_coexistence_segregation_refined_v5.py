#!/usr/bin/env python3
"""Final strict refinement of the coexistence (C) / segregation (S) audit.

Adds a general non-display exclusion for studies whose colour polymorphism target is
androecium/anther/pollen colour rather than floral display surfaces. This removes sexual-
organ colour cases without changing the positive C/S evidence rules established in v4.
"""
from __future__ import annotations

import re

import build_v22_coexistence_segregation_refined_v4 as v4

base = v4.base

SEXUAL_ORGAN_COLOUR_ONLY = (
    r"(?:androecium|androecial|anther(?:/pollen)?|pollen)\s+colou?r"
    r"|colou?r\s+(?:of\s+)?(?:the\s+)?(?:androecium|androecial|anther(?:s)?|pollen)"
    r"|(?:androecium|anther|pollen)\s+colou?r\s+polymorph\w*"
)

# Treat explicit sexual-organ colour targets as non-display conflicts. The base builder
# still requires display/discrete context for C/S, but this guard prevents incidental
# mentions of corolla/floral traits from admitting anther/pollen colour studies.
base.HARD_CONFLICT_RE = re.compile(
    rf"(?:{base.HARD_CONFLICT_RE.pattern})|(?:{SEXUAL_ORGAN_COLOUR_ONLY})",
    flags=re.IGNORECASE | re.DOTALL,
)

if __name__ == "__main__":
    base.main()
