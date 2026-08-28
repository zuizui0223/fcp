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

# The inherited pattern starts with inline (?is) flags. Remove only that leading flag
# token before combining it with the new alternative, then compile with equivalent
# explicit flags. The biological rule itself is unchanged.
_inherited_conflict = base.HARD_CONFLICT_RE.pattern
if _inherited_conflict.startswith("(?is)"):
    _inherited_conflict = _inherited_conflict[5:]

base.HARD_CONFLICT_RE = re.compile(
    rf"(?:{_inherited_conflict})|(?:{SEXUAL_ORGAN_COLOUR_ONLY})",
    flags=re.IGNORECASE | re.DOTALL,
)

if __name__ == "__main__":
    base.main()
