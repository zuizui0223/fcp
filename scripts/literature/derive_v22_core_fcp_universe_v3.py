#!/usr/bin/env python3
"""Tighten the primary display-colour FCP core without using C/S outcomes.

This wrapper removes `population` as a sufficient display context for generic colour-morph
wording and broadens the sexual-organ-only exclusion. Expanded-universe membership is
unchanged; this only sharpens the primary display-colour core.
"""
from __future__ import annotations

import re

import derive_v22_core_fcp_universe as base

# Replace the two generic context rules: population alone is not a floral-display structure.
strict_patterns=[]
for name, pattern in base.DIRECT_PATTERNS:
    if name in {"colour_morphs_with_display_context", "display_context_with_colour_morphs"}:
        continue
    strict_patterns.append((name, pattern))
strict_patterns.extend([
    ("colour_morphs_with_strict_display_context", re.compile(
        r"(?is)\bcolou?r\s+morphs?\b.{0,120}\b(?:flower|floral|petal|corolla|perianth|tepal|labellum|bract|inflorescence)\b"
    )),
    ("strict_display_context_with_colour_morphs", re.compile(
        r"(?is)\b(?:flower|floral|petal|corolla|perianth|tepal|labellum|bract|inflorescence)\b.{0,120}\bcolou?r\s+morphs?\b"
    )),
])
base.DIRECT_PATTERNS = strict_patterns

base.NONDISPLAY_ONLY_RE = re.compile(
    r"(?is)\b(?:stigma|stigmatic|gynoecium|gynoecial|anther|androecium|androecial|pollen)\b.{0,35}\bcolou?r\b"
    r"|\bcolou?r\b.{0,35}\b(?:stigma|stigmatic|gynoecium|gynoecial|anther|androecium|androecial|pollen)\b"
)
base.DISPLAY_STRUCTURE_RE = re.compile(
    r"(?is)\b(?:petal|corolla|perianth|tepal|labellum|bract|inflorescence)\s+colou?r\b|"
    r"\b(?:flower|floral)[- ]?colou?r\s+(?:morphs?|forms?|variants?|polymorph\w*)\b"
)

if __name__ == "__main__":
    base.main()
