#!/usr/bin/env python3
"""Adjudication-ready spatial-scale extraction built on PR #12.

This stricter layer keeps the systematic corpus and PR #12 decision vocabulary, but
adds the reviewer-facing safeguards needed before ecological modelling:

* explicit source_role gating (secondary sources discover candidates but cannot vote);
* hard exclusion of artificial/induced/developmental records;
* title attribution without dilution, abstract attribution only when K <= threshold;
* species + direct flower-colour + direction proximity within one local window;
* species-level P0--P3 evidence hierarchy and mixed/conflict retention;
* immutable evidence ledger, reviewer-1/reviewer-2/adjudication template;
* SHA-256 pre-freeze manifest; no dataset is labelled frozen before adjudication.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

from build_evidence_review_queue import (
    ARTIFICIAL_TERMS,
    DIRECT_POLYMORPHISM_RE,
    NATURAL_TERMS,
    classify,
    count_terms,
    evidence_class,
)
from build_merged_spatial_classification import (
    COLOUR,
    DIRECT,
    NEGATIVE,
    PROXIMITY_WINDOW,
    load_vocabulary,
    norm,
    score_attribution,
    spatial_class,
)
from build_provisional_spatial_classification import BINOMIAL

REQUIRED = {"screen_priority", "title", "abstract", "candidate_species"}
PRIMARY_ROLES = {"primary_evidence", "primary_study", "primary", "empirical_primary"}
DISCOVERY_ROLES = {"review_discovery_seed", "secondary_synthesis", "discovery_only", "secondary"}
DISCOVERY_WORK_TYPES = {
    "review", "book", "book-chapter", "book chapter", "editorial", "letter", "erratum",
    "meta-analysis", "systematic-review",
}
VALID_LABELS = {"within_population", "among_population", "mixed", "unclear"}


def split_names(value: str) -> list[str]:
    return [norm(v) for v in str(value or "").split(";") if norm(v) and BINOMIAL.fullmatch(norm(v))]


def normalise_role(row: dict) -> tuple[str, str]:
    explicit = norm(row.get("source_role")).lower().replace("-", "_").replace(" ", "_")
    if explicit:
        if explicit in DISCOVERY_ROLES:
            return "discovery_only", "explicit_source_role"
        if explicit in PRIMARY_ROLES:
            return "primary_evidence", "explicit_source_role"
        return explicit, "explicit_source_role_unrecognised"
    work_type = norm(row.get("work_type")).lower()
    if work_type in DISCOVERY_WORK_TYPES:
        return "discovery_only", "inferred_from_work_type"
    return "primary_evidence", "inferred_from_work_type"


def local_windows(text: str, species: str, radius: int) -> list[str]:
    hits = list(re.finditer(re.escape(species), text, flags=re.I))
    return [text[max(0, h.start() - radius): min(len(text), h.end() + radius)] for h in hits]


def directional_signal(
    text: str,
    species: str,
    pattern: re.Pattern[str],
    radius: int,
) -> bool:
    """Require species, direct flower-colour evidence and direction in one local window."""
    for window in local_windows(text, species, radius):
        direct = bool(DIRECT.search(window)) or bool(DIRECT_POLYMORPHISM_RE.search(window))
        if direct and COLOUR.search(window) and pattern.search(window):
            return True
    return False


def record_is_artificial(text: str) -> bool:
    return bool(NEGATIVE.search(text)) or count_terms(text, ARTIFICIAL_TERMS) > 0


def read_queue(path: Path) -> tuple[list[dict], list[str]]:
    records: list[dict] = []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or [])
        missing = sorted(REQUIRED - set(fields))
        if missing:
            raise SystemExit(f"screening queue missing required columns: {missing}; present={fields}")
        for row in reader:
            if row.get("screen_priority") not in {"P1_high_natural_itv", "P1_high_population_itv"}:
                continue
            title = norm(row.get("title"))
            abstract = norm(row.get("abstract"))
            names = split_names(row.get("candidate_species"))
            role, role_basis = normalise_role(row)
            records.append({
                "record_id": norm(row.get("record_id")) or norm(row.get("dedup_key")),
                "doi": norm(row.get("doi")),
                "title": title,
                "abstract": abstract,
                "text": f"{title} {abstract}",
                "work_type": norm(row.get("work_type")).lower(),
                "source_role": role,
                "source_role_basis": role_basis,
                "names": names,
                "title_names": [n for n in names if re.search(re.escape(n), title, re.I)],
            })
    return records, fields


def collect_evidence(records: list[dict], within_re, among_re, cap: int, radius: int):
    per_species: dict[str, list[dict]] = defaultdict(list)
    ledger: list[dict] = []
    diag: Counter = Counter()

    for record in records:
        diag["records_read"] += 1
        artificial = record_is_artificial(record["text"])
        if artificial:
            diag["records_hard_excluded_artificial"] += 1
        discovery_only = record["source_role"] == "discovery_only"
        if discovery_only:
            diag["records_discovery_only"] += 1

        n_names = len(set(record["names"]))
        for species in record["names"]:
            in_title = species in record["title_names"]
            attribution_allowed = in_title or n_names <= cap
            within_raw = directional_signal(record["text"], species, within_re, radius)
            among_raw = directional_signal(record["text"], species, among_re, radius)
            direction_allowed = attribution_allowed and not artificial and not discovery_only
            within = within_raw and direction_allowed
            among = among_raw and direction_allowed
            direct_local = any(
                (DIRECT.search(w) or DIRECT_POLYMORPHISM_RE.search(w)) and COLOUR.search(w)
                for w in local_windows(record["text"], species, radius)
            )
            natural = count_terms(record["text"], NATURAL_TERMS) > 0

            exclusion_reason = ""
            if not attribution_allowed:
                exclusion_reason = "abstract_species_dilution"
                diag["attributions_suppressed_dilution"] += 1
            elif artificial:
                exclusion_reason = "hard_exclusion_artificial_or_induced"
            elif discovery_only:
                exclusion_reason = "source_role_discovery_only"
            elif not direct_local:
                exclusion_reason = "no_local_direct_colour_signal"
            elif not (within_raw or among_raw):
                exclusion_reason = "no_local_direction_signal"

            evidence = {
                "record_id": record["record_id"],
                "doi": record["doi"],
                "title": record["title"],
                "canonical_name": species,
                "in_title": int(in_title),
                "n_binomials_in_record": n_names,
                "dilution_cap": cap,
                "attribution_allowed": int(attribution_allowed),
                "source_role": record["source_role"],
                "source_role_basis": record["source_role_basis"],
                "work_type": record["work_type"],
                "artificial_signal": int(artificial),
                "direct_colour_signal_local": int(bool(direct_local)),
                "within_signal_raw": int(within_raw),
                "among_signal_raw": int(among_raw),
                "within_vote": int(within),
                "among_vote": int(among),
                "exclusion_reason": exclusion_reason,
            }
            ledger.append(evidence)

            # Discovery-only records still preserve candidate provenance but do not vote.
            if attribution_allowed and not artificial and direct_local:
                signals = {
                    "direct": bool(direct_local),
                    "artificial": artificial,
                    "within": within,
                    "among": among,
                    "natural": natural,
                    "discovery_only": discovery_only,
                }
                per_species[species].append({
                    "record": {
                        "record_id": record["record_id"],
                        "doi": record["doi"],
                        "title": record["title"],
                        "snippet": record["abstract"],
                        "family_hint": "",
                    },
                    "signals": signals,
                    "in_title": in_title,
                    "source_role": record["source_role"],
                })
                diag["attributions_retained"] += 1
    return per_species, ledger, diag


def aggregate_species(name: str, hits: list[dict]) -> dict:
    scored = sorted(hits, key=lambda h: -score_attribution(h["in_title"], h["signals"]))
    best = scored[0]
    scores = [score_attribution(h["in_title"], h["signals"]) for h in hits]
    primary = [h for h in hits if h["source_role"] != "discovery_only"]
    return {
        "canonical_name": name,
        "family": "",
        "n_works": len({h["record"]["record_id"] for h in primary}),
        "n_candidate_discovery_works": len({h["record"]["record_id"] for h in hits if h["source_role"] == "discovery_only"}),
        "n_title_matches": sum(h["in_title"] for h in primary),
        "n_context_matches": sum(not h["in_title"] for h in primary),
        "max_score": max(scores),
        "total_score": sum(scores),
        "best_title": best["record"]["title"],
        "best_doi": best["record"]["doi"],
        "best_openalex_id": best["record"]["record_id"],
        "best_match_evidence": best["record"]["snippet"],
        "followup_evidence_count": len(primary),
        "followup_direct_count": sum(h["signals"]["direct"] for h in primary),
        "followup_natural_count": sum(h["signals"]["natural"] for h in primary),
        "followup_artificial_count": sum(h["signals"]["artificial"] for h in primary),
        "n_within_records": sum(h["signals"]["within"] for h in primary),
        "n_among_records": sum(h["signals"]["among"] for h in primary),
        "n_title_within_records": sum(h["in_title"] and h["signals"]["within"] for h in primary),
        "n_title_among_records": sum(h["in_title"] and h["signals"]["among"] for h in primary),
    }


def write_csv(path: Path, rows: list[dict], fields: Iterable[str] | None = None) -> None:
    fieldnames = list(fields or (list(rows[0]) if rows else ["canonical_name"]))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def review_template(rows: list[dict]) -> list[dict]:
    return [{
        "canonical_name": r["canonical_name"],
        "automated_evidence_class": r["evidence_class"],
        "automated_spatial_scale": r["spatial_scale"],
        "review_priority": r["review_priority"],
        "best_doi": r["best_doi"],
        "reviewer_1_evidence_class": "",
        "reviewer_1_spatial_scale": "",
        "reviewer_1_notes": "",
        "reviewer_2_evidence_class": "",
        "reviewer_2_spatial_scale": "",
        "reviewer_2_notes": "",
        "adjudicated_evidence_class": "",
        "adjudicated_spatial_scale": "",
        "adjudicator": "",
        "adjudication_reason": "",
        "review_status": "unreviewed",
        "eligible_for_freeze": "false",
    } for r in rows]


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--corpus", required=True)
    p.add_argument("--search-config", default="literature/itv_fcp_search_config.json")
    p.add_argument("--dilution-cap", type=int, default=3)
    p.add_argument("--proximity-window", type=int, default=PROXIMITY_WINDOW)
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    within_re, among_re = load_vocabulary(Path(args.search_config))
    records, input_fields = read_queue(Path(args.corpus))
    per_species, ledger, diag = collect_evidence(
        records, within_re, among_re, args.dilution_cap, args.proximity_window
    )

    rows: list[dict] = []
    for name, hits in per_species.items():
        row = aggregate_species(name, hits)
        klass, natural, artificial, colour_change = evidence_class(row)
        priority, reason = classify(row, klass)
        scale, source = spatial_class(row)
        row.update({
            "evidence_class": klass,
            "natural_signal_count": natural,
            "artificial_signal_count": artificial,
            "colour_change_signal_count": colour_change,
            "review_priority": priority,
            "review_reason": reason,
            "spatial_scale": scale,
            "classification_source": source,
            "review_status": "unreviewed",
            "eligible_for_freeze": "false",
        })
        rows.append(row)

    rows.sort(key=lambda r: (r["review_priority"], -r["max_score"], r["canonical_name"]))
    ledger.sort(key=lambda r: (r["canonical_name"], r["record_id"]))
    class_path = outdir / "audited_spatial_classification_unreviewed.csv"
    ledger_path = outdir / "record_species_evidence_ledger.csv"
    review_path = outdir / "double_review_adjudication_template.csv"
    correction_path = outdir / "correction_log_template.csv"
    write_csv(class_path, rows)
    write_csv(ledger_path, ledger)
    write_csv(review_path, review_template(rows))
    write_csv(correction_path, [], [
        "canonical_name", "old_label", "new_label", "reason", "reviewer", "date", "source_doi"
    ])

    prefreeze = {
        "status": "prefreeze_only",
        "eligible_for_freeze": False,
        "reason": "All species remain unreviewed; codebook requires adjudicated status.",
        "input_corpus": str(args.corpus),
        "input_fields": input_fields,
        "parameters": {
            "dilution_cap": args.dilution_cap,
            "proximity_window": args.proximity_window,
            "source_role_policy": "explicit source_role, otherwise deterministic work_type fallback",
        },
        "counts": {
            "records": len(records),
            "species": len(rows),
            "spatial_scale": dict(Counter(r["spatial_scale"] for r in rows)),
            "review_priority": dict(Counter(r["review_priority"] for r in rows)),
            "source_role": dict(Counter(r["source_role"] for r in ledger)),
            "diagnostics": dict(diag),
        },
        "sha256": {
            class_path.name: sha256(class_path),
            ledger_path.name: sha256(ledger_path),
            review_path.name: sha256(review_path),
            correction_path.name: sha256(correction_path),
        },
    }
    (outdir / "prefreeze_manifest.json").write_text(json.dumps(prefreeze, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(prefreeze, indent=2))


if __name__ == "__main__":
    main()
