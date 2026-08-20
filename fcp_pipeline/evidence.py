"""Shared evidence-classification helpers for the frozen paper pipeline."""
from __future__ import annotations

import re

WITHIN = re.compile(
    r"\b(within[- ]populations?|same populations?|coexist|co-occurr|morph frequenc|"
    r"frequency[- ]dependent|polymorphic populations?|multiple (?:colou?r|flower) morphs|"
    r"colour morphs|color morphs)\b",
    re.I,
)
GEOGRAPHIC = re.compile(
    r"\b(geographic|spatial variation|among populations|between populations|population differentiation|"
    r"cline|hybrid zone|range edge|local adaptation|translocation|regional variation)\b",
    re.I,
)
MAINTENANCE = re.compile(
    r"\b(maintenance|balancing selection|frequency[- ]dependent|pollinator[- ]mediated|"
    r"assortative mating|negative frequency)\b",
    re.I,
)

ALLOWED_SPATIAL_LABELS = {"within_population", "among_population", "mixed", "unclear"}


def rule_label(text: str) -> tuple[str, int, int]:
    """Return deterministic spatial label plus within/geographic signal flags."""
    text = str(text or "")
    within = int(bool(WITHIN.search(text)))
    geographic = int(bool(GEOGRAPHIC.search(text)))
    if within and not geographic:
        label = "within_population"
    elif geographic and not within:
        label = "among_population"
    elif within and geographic:
        label = "mixed"
    else:
        label = "unclear"
    return label, within, geographic


def normalize_doi(value: str) -> str:
    x = str(value or "").strip().lower()
    return re.sub(r"^https?://(?:dx\.)?doi\.org/", "", x)


def clean_excerpt(value: str, limit: int = 1800) -> str:
    return re.sub(r"\s+", " ", str(value or "")).strip()[:limit]


def source_match(source_id: str, doi: str, openalex_id: str) -> str:
    source = str(source_id or "").strip()
    if not source:
        return "classification_source_missing"
    if normalize_doi(source) and normalize_doi(source) == normalize_doi(doi):
        return "matches_queue_best_doi"
    if source.rstrip("/").lower() == str(openalex_id or "").strip().rstrip("/").lower():
        return "matches_queue_best_openalex"
    return "classification_source_differs_from_queue_best"
