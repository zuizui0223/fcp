#!/usr/bin/env python3
"""Tighten the primary display-colour FCP core without using C/S outcomes.

This wrapper removes `population` as a sufficient display context for generic colour-morph
wording and prevents sexual-organ-target studies from re-entering via generic background
phrases such as `flower color polymorphism` in the abstract. Expanded-universe membership
is unchanged; this only sharpens the primary display-colour core.
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
TARGET_NONDISPLAY_TITLE_RE = re.compile(
    r"(?is)\b(?:stigma|stigmatic|gynoecium|gynoecial|anther|androecium|androecial|pollen)\b.{0,50}\bcolou?r\b"
    r"|\bcolou?r\b.{0,50}\b(?:stigma|stigmatic|gynoecium|gynoecial|anther|androecium|androecial|pollen)\b"
)
TITLE_DISPLAY_RE = re.compile(
    r"(?is)\b(?:petal|corolla|perianth|tepal|labellum|bract|inflorescence|flower|floral)\b.{0,35}\bcolou?r\b"
    r"|\bcolou?r\b.{0,35}\b(?:petal|corolla|perianth|tepal|labellum|bract|inflorescence|flower|floral)\b"
)

_original_direct_reason = base.direct_reason

def direct_reason(text: str, title: str) -> str:
    # If the paper title explicitly defines the target as a sexual-organ colour trait,
    # generic background language in the abstract cannot turn it into a display-colour study.
    if TARGET_NONDISPLAY_TITLE_RE.search(title or "") and not TITLE_DISPLAY_RE.search(title or ""):
        return ""
    return _original_direct_reason(text, title)

base.direct_reason = direct_reason

if __name__ == "__main__":
    base.main()
