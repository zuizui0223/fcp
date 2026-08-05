#!/usr/bin/env python3
"""Merged spatial-scale extraction: systematic-search breadth with baseline decision rules.

Path A (34-species baseline) and path B (107-species systematic pipeline) each hold half
of a defensible extraction. This module keeps A's decision layer verbatim -- the tiered
evidence thresholds, the richer directional vocabulary, the conflict gate and the freeze --
and applies it to whichever corpus is supplied, including B's systematic retrieval.

Reused unchanged from the repository:
  - evidence_class() / classify()      tiered evidence thresholds        (build_evidence_review_queue)
  - colour_linked_geographic_signal()  colour/geography proximity test   (enrich_spatial_scale_evidence)
  - WITHIN / GEOGRAPHIC / COLOUR / DIRECT / NEGATIVE regexes             (enrich_spatial_scale_evidence)
  - BINOMIAL / norm()                  binomial recognition              (build_provisional_spatial_classification)

New here, and deliberately small:
  1. a corpus adapter so one decision layer can read either retrieval;
  2. the A-union-B directional vocabulary;
  3. a symmetric proximity test, so `within` is judged as strictly as `among`;
  4. species attribution with a dilution cap on multi-species abstracts;
  5. aggregation into the schema evidence_class()/classify() already expect;
  6. the conflict gate and freeze, ported from analysis_evidence_spatial_scale_enriched.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

from build_evidence_review_queue import (
    ARTIFICIAL_TERMS,
    DIRECT_POLYMORPHISM_RE,
    NATURAL_TERMS,
    classify,
    count_terms,
    evidence_class,
    to_int,
)
from build_provisional_spatial_classification import BINOMIAL, norm
from enrich_spatial_scale_evidence import (
    COLOUR,
    DIRECT,
    GEOGRAPHIC,
    NEGATIVE,
    WITHIN,
    colour_linked_geographic_signal,
)

PROXIMITY_WINDOW = 180
DEFAULT_SPECIES_DILUTION_CAP = 3


# --- 2. directional vocabulary: A's regexes unioned with B's configured terms ----------

def union_pattern(base: re.Pattern[str], extra_terms: list[str]) -> re.Pattern[str]:
    """Widen a path-A regex with path-B's configured vocabulary.

    B stores plain substrings; they are escaped and allowed to hyphenate so that
    `within population` also matches `within-population`.
    """
    parts = [base.pattern]
    for term in extra_terms:
        parts.append(re.escape(term.strip()).replace(r"\ ", r"[- ]").replace(r"\-", r"[- ]"))
    return re.compile("|".join(p for p in parts if p), re.I)


def load_vocabulary(config_path: Path | None) -> tuple[re.Pattern[str], re.Pattern[str]]:
    if config_path is None or not config_path.exists():
        return WITHIN, GEOGRAPHIC
    config = json.loads(config_path.read_text(encoding="utf-8"))
    terms = config.get("screening_terms") or config.get("terms") or {}
    return (
        union_pattern(WITHIN, list(terms.get("within") or [])),
        union_pattern(GEOGRAPHIC, list(terms.get("among") or [])),
    )


# --- 3. symmetric proximity -----------------------------------------------------------

def colour_linked(text: str, pattern: re.Pattern[str], window: int = PROXIMITY_WINDOW) -> bool:
    """Generalisation of colour_linked_geographic_signal to either direction.

    Path A applies the proximity requirement only to `among`, leaving `within` on a bare
    bag-of-words match; path B applies it to neither. Judging the two asymmetrically is
    what lets unrelated population language classify a record.
    """
    colour_hits = list(COLOUR.finditer(text))
    direction_hits = list(pattern.finditer(text))
    return any(
        max(c.start(), d.start()) - min(c.end(), d.end()) <= window
        for c in colour_hits
        for d in direction_hits
    )


def self_test(texts: list[str]) -> None:
    """colour_linked must reproduce the original exactly on the geographic direction."""
    mismatches = [
        t for t in texts
        if colour_linked(t, GEOGRAPHIC) != colour_linked_geographic_signal(t)
    ]
    if mismatches:
        raise SystemExit(
            f"colour_linked diverges from colour_linked_geographic_signal on "
            f"{len(mismatches)}/{len(texts)} texts"
        )
    print(f"self-test ok: generalised proximity matches the original on {len(texts)} texts")


# --- 1. corpus adapters ---------------------------------------------------------------

def _split(value: str) -> list[str]:
    return [norm(v) for v in str(value or "").split(";") if norm(v)]


def _evidence_snippet(value: str) -> str:
    """match_evidence is `species|score|reason_tokens|snippet`; keep the snippet."""
    parts = str(value or "").split("|")
    return parts[-1] if len(parts) >= 4 else str(value or "")


def read_global_works(path: Path) -> list[dict]:
    """Path-A discovery corpus (data/global_flower_colour_works.csv)."""
    records = []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            title = norm(row.get("title"))
            snippet = _evidence_snippet(row.get("match_evidence"))
            title_species = _split(row.get("title_matched_species"))
            context_species = _split(row.get("context_matched_species"))
            records.append(
                {
                    "record_id": norm(row.get("openalex_id")),
                    "title": title,
                    "text": f"{title} {snippet}",
                    "doi": norm(row.get("doi")),
                    "family_hint": (_split(row.get("matched_families")) or [""])[0],
                    "title_species": title_species,
                    "context_species": [s for s in context_species if s not in title_species],
                    "n_species_mentioned": max(
                        to_int(row.get("raw_species_mentions")),
                        len(set(title_species) | set(context_species)),
                    ),
                    "snippet": snippet,
                }
            )
    return records


SCREENING_QUEUE_REQUIRED = {"screen_priority", "title", "abstract", "candidate_species"}


def read_screening_queue(path: Path) -> list[dict]:
    """Path-B systematic screening queue (itv_fcp_human_screening_queue.csv).

    `abstract` is load-bearing: the proximity test needs more than a title, and a queue
    without it would silently classify everything from headline words alone.
    """
    records = []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        missing = sorted(SCREENING_QUEUE_REQUIRED - set(reader.fieldnames or []))
        if missing:
            raise SystemExit(
                f"screening queue is missing required columns: {missing}; "
                f"present columns: {sorted(reader.fieldnames or [])}"
            )
        for row in reader:
            if row.get("screen_priority") not in {"P1_high_natural_itv", "P1_high_population_itv"}:
                continue
            title = norm(row.get("title"))
            abstract = norm(row.get("abstract"))[:1200]
            names = [n for n in _split(row.get("candidate_species")) if BINOMIAL.fullmatch(n)]
            lowered_title = title.lower()
            title_species = [n for n in names if n.lower() in lowered_title]
            lowered_abstract = abstract.lower()
            context_species = [
                n for n in names if n not in title_species and n.lower() in lowered_abstract
            ]
            records.append(
                {
                    "record_id": norm(row.get("record_id")) or norm(row.get("dedup_key")),
                    "title": title,
                    "text": f"{title} {abstract}",
                    "doi": norm(row.get("doi")),
                    "family_hint": "",
                    "title_species": title_species,
                    "context_species": context_species,
                    "n_species_mentioned": len(names),
                    "snippet": abstract,
                    "work_type": norm(row.get("work_type")).lower(),
                }
            )
    return records


CORPUS_READERS = {"global_works": read_global_works, "screening_queue": read_screening_queue}

# Codebook: review_discovery_seed and secondary_synthesis are `discovery_only`; such
# sources may surface a candidate species but must not supply directional evidence.
DISCOVERY_ONLY_WORK_TYPES = {"review", "book", "book-chapter", "editorial", "letter", "erratum"}


# --- 4/5. attribution and aggregation -------------------------------------------------

def record_signals(record: dict, within_re: re.Pattern[str], among_re: re.Pattern[str]) -> dict:
    text = record["text"]
    direct = bool(DIRECT.search(text)) or bool(DIRECT_POLYMORPHISM_RE.search(text))
    artificial = bool(NEGATIVE.search(text)) or count_terms(text, ARTIFICIAL_TERMS) > 0
    discovery_only = record.get("work_type", "") in DISCOVERY_ONLY_WORK_TYPES
    within = (not discovery_only) and colour_linked(text, within_re)
    among = (not discovery_only) and colour_linked(text, among_re)
    return {
        "direct": direct,
        "artificial": artificial,
        "within": within,
        "among": among,
        "natural": count_terms(text, NATURAL_TERMS) > 0,
        "discovery_only": discovery_only,
    }


def score_attribution(in_title: bool, signals: dict) -> int:
    """Scoring shape reused from enrich_spatial_scale_evidence.main()."""
    return (
        (12 if in_title else 5)
        + (10 if signals["direct"] else 0)
        + (8 if signals["within"] else 0)
        + (8 if signals["among"] else 0)
        - (10 if signals["artificial"] else 0)
    )


def attribute(records: list[dict], within_re, among_re, dilution_cap: int) -> tuple[dict, Counter]:
    """Assign each record's evidence to species, capping multi-species abstracts.

    Path B hands a record's directional signal to every binomial in the title or the first
    1200 characters of the abstract, with no check that the species is the study subject.
    Title attributions are kept unconditionally; abstract-only attributions are dropped
    once a record mentions more species than the cap, which is where comparative reviews
    and congener asides live.
    """
    per_species: dict[str, list[dict]] = defaultdict(list)
    diag: Counter = Counter()
    for record in records:
        diag["records_read"] += 1
        signals = record_signals(record, within_re, among_re)
        if signals["discovery_only"]:
            diag["records_discovery_only"] += 1
        if not (signals["direct"] or signals["within"] or signals["among"]):
            diag["records_without_usable_signal"] += 1
            continue
        diag["records_retained"] += 1
        diluted = record["n_species_mentioned"] > dilution_cap
        if diluted:
            diag["records_over_dilution_cap"] += 1
        for name in record["title_species"]:
            per_species[name].append({"record": record, "signals": signals, "in_title": True})
            diag["attributions_title"] += 1
        if diluted:
            diag["attributions_context_suppressed"] += len(record["context_species"])
            continue
        for name in record["context_species"]:
            per_species[name].append({"record": record, "signals": signals, "in_title": False})
            diag["attributions_context"] += 1
    return per_species, diag


def aggregate(name: str, hits: list[dict]) -> dict:
    """Build exactly the row shape evidence_class()/classify() already consume."""
    scored = sorted(
        hits, key=lambda h: -score_attribution(h["in_title"], h["signals"])
    )
    best = scored[0]
    scores = [score_attribution(h["in_title"], h["signals"]) for h in hits]
    return {
        "canonical_name": name,
        "family": next((h["record"]["family_hint"] for h in hits if h["record"]["family_hint"]), ""),
        "n_works": len({h["record"]["record_id"] for h in hits}),
        "n_title_matches": sum(1 for h in hits if h["in_title"]),
        "n_context_matches": sum(1 for h in hits if not h["in_title"]),
        "max_score": max(scores),
        "total_score": sum(scores),
        "best_title": best["record"]["title"],
        "best_doi": best["record"]["doi"],
        "best_openalex_id": best["record"]["record_id"],
        "best_match_evidence": best["record"]["snippet"],
        "followup_evidence_count": len(hits),
        "followup_direct_count": sum(1 for h in hits if h["signals"]["direct"]),
        "followup_natural_count": sum(1 for h in hits if h["signals"]["natural"]),
        "followup_artificial_count": sum(1 for h in hits if h["signals"]["artificial"]),
        "n_within_records": sum(1 for h in hits if h["signals"]["within"]),
        "n_among_records": sum(1 for h in hits if h["signals"]["among"]),
        "n_title_within_records": sum(1 for h in hits if h["in_title"] and h["signals"]["within"]),
        "n_title_among_records": sum(1 for h in hits if h["in_title"] and h["signals"]["among"]),
    }


# --- 6. spatial classification with the conflict gate ---------------------------------

def spatial_class(row: dict) -> tuple[str, str]:
    """Port of analysis_evidence_spatial_scale_enriched.integrate().

    The baseline call is made on records naming the species in the title; the remaining
    records act as the independent corroboration path A checked enrichment works against.
    A contradiction demotes the species to `mixed` rather than forcing a binary class.
    """
    baseline_within = row["n_title_within_records"] > 0
    baseline_among = row["n_title_among_records"] > 0
    other_within = row["n_within_records"] - row["n_title_within_records"] > 0
    other_among = row["n_among_records"] - row["n_title_among_records"] > 0

    if baseline_within and not baseline_among:
        return ("mixed", "baseline_enrichment_conflict") if other_among else (
            "within_population",
            "baseline_unambiguous",
        )
    if baseline_among and not baseline_within:
        return ("mixed", "baseline_enrichment_conflict") if other_within else (
            "among_population",
            "baseline_unambiguous",
        )
    if baseline_within and baseline_among:
        return "mixed", "baseline_mixed"
    if other_within and not other_among:
        return "within_population", "high_confidence_enrichment"
    if other_among and not other_within:
        return "among_population", "high_confidence_enrichment"
    if other_within and other_among:
        return "mixed", "high_confidence_enrichment_mixed"
    return "unclear", "unresolved"


# --- 7. freeze ------------------------------------------------------------------------

def write_frozen_manifest(rows: list[dict], path: Path) -> str:
    keep = sorted(
        (r for r in rows if r["classification_source"] == "baseline_unambiguous"
         and r["spatial_scale"] in {"within_population", "among_population"}),
        key=lambda r: r["canonical_name"],
    )
    fields = ["canonical_name", "family", "spatial_scale", "review_priority", "evidence_class",
              "n_works", "n_title_matches", "max_score", "best_doi", "best_openalex_id"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(keep)
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--corpus", required=True)
    parser.add_argument("--corpus-format", choices=sorted(CORPUS_READERS), default="global_works")
    parser.add_argument("--search-config", default="literature/itv_fcp_search_config.json")
    parser.add_argument("--dilution-cap", type=int, default=DEFAULT_SPECIES_DILUTION_CAP)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    within_re, among_re = load_vocabulary(Path(args.search_config))
    records = CORPUS_READERS[args.corpus_format](Path(args.corpus))
    if args.self_test:
        self_test([r["text"] for r in records])
    per_species, diag = attribute(records, within_re, among_re, args.dilution_cap)

    rows = []
    for name, hits in per_species.items():
        row = aggregate(name, hits)
        klass, natural, artificial, colour_change = evidence_class(row)
        priority, reason = classify(row, klass)
        scale, source = spatial_class(row)
        row.update(
            {
                "evidence_class": klass,
                "natural_signal_count": natural,
                "artificial_signal_count": artificial,
                "colour_change_signal_count": colour_change,
                "review_priority": priority,
                "review_reason": reason,
                "spatial_scale": scale,
                "classification_source": source,
                "review_status": "unreviewed",
            }
        )
        rows.append(row)

    rows.sort(key=lambda r: (r["review_priority"], -r["max_score"], r["canonical_name"]))
    classification_path = outdir / "merged_spatial_classification.csv"
    with classification_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else ["canonical_name"])
        writer.writeheader()
        writer.writerows(rows)

    digest = write_frozen_manifest(rows, outdir / "merged_baseline_unambiguous_manifest.csv")

    eligible = [
        r for r in rows
        if r["review_priority"] in {"P0", "P1", "P2"}
        and r["evidence_class"] == "natural_polymorphism"
        and r["spatial_scale"] in {"within_population", "among_population"}
    ]
    binary = Counter(r["spatial_scale"] for r in eligible)
    qc = {
        "status": "complete" if len(eligible) >= 20 and len(binary) == 2 else "insufficient",
        "corpus": args.corpus,
        "corpus_format": args.corpus_format,
        "dilution_cap": args.dilution_cap,
        "record_diagnostics": dict(diag),
        "species_total": len(rows),
        "review_priority": dict(Counter(r["review_priority"] for r in rows)),
        "evidence_class": dict(Counter(r["evidence_class"] for r in rows)),
        "spatial_scale_all": dict(Counter(r["spatial_scale"] for r in rows)),
        "classification_source": dict(Counter(r["classification_source"] for r in rows)),
        "model_eligible_species": len(eligible),
        "model_eligible_within": binary.get("within_population", 0),
        "model_eligible_among": binary.get("among_population", 0),
        "frozen_manifest_sha256": digest,
        "semantic_guard": (
            "Deterministic screening labels. review_status remains `unreviewed`; the codebook "
            "requires `adjudicated` before any set is eligible_for_freeze."
        ),
    }
    (outdir / "merged_classification_qc.json").write_text(
        json.dumps(qc, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(qc, indent=2))


if __name__ == "__main__":
    main()
