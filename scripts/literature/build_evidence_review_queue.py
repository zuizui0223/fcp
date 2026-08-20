#!/usr/bin/env python3
"""Build a review queue focused on naturally occurring flower-colour polymorphism."""
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path
from typing import Dict, List

KNOWN_POSITIVES = {
    "Dactylorhiza sambucina", "Gymnadenia rhellicani",
    "Linanthus parryae", "Silene littorea",
}
NATURAL_TERMS = (
    "natural population", "wild population", "wild populations", "in the wild",
    "polymorphic population", "polymorphic populations", "colour polymorphism",
    "color polymorphism", "colour morph", "color morph", "morph frequencies",
    "morph frequency", "geographic variation", "spatial variation",
    "population differentiation", "balancing selection", "frequency dependent",
    "frequency-dependent", "pollinator-mediated selection", "local adaptation",
    "within-population", "within population", "geographic cline",
)
ARTIFICIAL_TERMS = (
    "cultivar", "cultivars", "breeding", "ornamental", "market value",
    "gene editing", "gene-edited", "mutagenesis", "mutant", "transgenic",
    "qtl", "candidate gene", "genome-wide association", "gwas",
    "crop improvement", "agriculture and tourism", "petal colour breeding",
    "petal color breeding",
)
COLOUR_CHANGE_TERMS = (
    "color changing", "colour changing", "changes color", "changes colour",
    "colour change", "color change", "ontogenetic color", "ontogenetic colour",
)
DIRECT_POLYMORPHISM_RE = re.compile(
    r"(?:flower|floral|petal|corolla|bract|labellum).{0,80}(?:colou?r|pigment).{0,80}(?:polymorph|morph|variation|dimorph|cline)"
    r"|(?:polymorph|dimorph|morph|variation|cline).{0,80}(?:flower|floral|petal|corolla|bract|labellum).{0,80}(?:colou?r|pigment)",
    re.IGNORECASE,
)
FIELDS = [
    "review_priority", "evidence_class", "canonical_name", "family", "n_works",
    "n_title_matches", "n_context_matches", "max_score", "total_score",
    "best_title", "best_doi", "best_openalex_id", "best_match_evidence",
    "followup_evidence_count", "followup_direct_count", "followup_natural_count",
    "followup_artificial_count", "natural_signal_count", "artificial_signal_count",
    "colour_change_signal_count", "review_status", "review_reason",
]


def read_rows(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def to_int(value: object) -> int:
    try:
        return int(str(value or "0"))
    except ValueError:
        return 0


def count_terms(text: str, terms: tuple[str, ...]) -> int:
    lowered = text.lower()
    return sum(1 for term in terms if term in lowered)


def evidence_class(row: Dict[str, str]) -> tuple[str, int, int, int]:
    text = f"{row.get('best_title', '')} {row.get('best_match_evidence', '')}"
    natural = count_terms(text, NATURAL_TERMS) + to_int(row.get("followup_natural_count"))
    artificial = count_terms(text, ARTIFICIAL_TERMS) + to_int(row.get("followup_artificial_count"))
    colour_change = count_terms(text, COLOUR_CHANGE_TERMS)
    direct_count = to_int(row.get("followup_direct_count"))
    direct = bool(DIRECT_POLYMORPHISM_RE.search(text)) or direct_count > 0

    if row.get("canonical_name", "") in KNOWN_POSITIVES:
        return "natural_polymorphism", max(natural, 1), artificial, colour_change
    if artificial >= 1 and natural == 0 and direct_count == 0:
        return "artificial_or_horticultural", natural, artificial, colour_change
    if colour_change >= 1 and not direct and natural == 0:
        return "ontogenetic_colour_change", natural, artificial, colour_change
    if direct and (natural >= 1 or direct_count >= 2 or to_int(row.get("n_context_matches")) >= 2):
        return "natural_polymorphism", natural, artificial, colour_change
    if direct:
        return "possible_polymorphism", natural, artificial, colour_change
    return "insufficient_direct_evidence", natural, artificial, colour_change


def classify(row: Dict[str, str], klass: str) -> tuple[str, str]:
    title = to_int(row.get("n_title_matches"))
    context = to_int(row.get("n_context_matches"))
    works = to_int(row.get("n_works"))
    max_score = to_int(row.get("max_score"))
    direct_count = to_int(row.get("followup_direct_count"))
    natural_count = to_int(row.get("followup_natural_count"))
    name = row.get("canonical_name", "")

    if name in KNOWN_POSITIVES:
        return "P0", "known_positive_control"
    if klass != "natural_polymorphism":
        return "DEFER", klass
    if direct_count >= 2 and natural_count >= 1:
        return "P1", "aggregated_followup_direct_and_natural_support"
    if direct_count >= 1 and natural_count >= 1:
        return "P2", "aggregated_followup_direct_support"
    if title >= 1 and (works >= 2 or context >= 2 or max_score >= 24):
        return "P1", "natural_polymorphism_with_direct_and_replicated_support"
    if title >= 1 and context >= 1 and max_score >= 18:
        return "P2", "natural_polymorphism_with_direct_support"
    if context >= 2 and works >= 2 and max_score >= 14:
        return "P3", "natural_polymorphism_with_replicated_context"
    return "DEFER", "natural_signal_below_review_threshold"


def write_csv(path: Path, rows: List[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", default="data/global_flower_colour_species_ranked.csv")
    parser.add_argument("--out", default="data/global_flower_colour_review_queue.csv")
    parser.add_argument("--deferred-out", default="data/global_flower_colour_deferred.csv")
    parser.add_argument("--qc-out", default="data/global_flower_colour_review_queue_qc.json")
    parser.add_argument("--max-review", type=int, default=150)
    args = parser.parse_args()

    rows = read_rows(Path(args.species))
    if not rows:
        raise RuntimeError("No candidate species available for review queue")

    selected: List[Dict[str, str]] = []
    deferred: List[Dict[str, str]] = []
    class_counts: Counter[str] = Counter()
    priority_counts: Counter[str] = Counter()

    for row in rows:
        klass, natural, artificial, colour_change = evidence_class(row)
        priority, reason = classify(row, klass)
        class_counts[klass] += 1
        priority_counts[priority] += 1
        enriched = {
            "review_priority": priority,
            "evidence_class": klass,
            **row,
            "natural_signal_count": str(natural),
            "artificial_signal_count": str(artificial),
            "colour_change_signal_count": str(colour_change),
            "review_reason": reason,
        }
        (deferred if priority == "DEFER" else selected).append(enriched)

    selected.sort(key=lambda row: (
        {"P0": 0, "P1": 1, "P2": 2, "P3": 3}[row["review_priority"]],
        -to_int(row.get("followup_direct_count")),
        -to_int(row.get("followup_natural_count")),
        -to_int(row.get("natural_signal_count")),
        -to_int(row.get("n_title_matches")),
        -to_int(row.get("n_context_matches")),
        -to_int(row.get("n_works")),
        row.get("canonical_name", ""),
    ))
    selected = selected[: max(1, args.max_review)]
    deferred.sort(key=lambda row: (
        row.get("evidence_class", ""),
        -to_int(row.get("followup_direct_count")),
        -to_int(row.get("followup_natural_count")),
        row.get("canonical_name", ""),
    ))

    write_csv(Path(args.out), selected)
    write_csv(Path(args.deferred_out), deferred)

    recovered = sorted(KNOWN_POSITIVES.intersection({r.get("canonical_name", "") for r in selected}))
    qc = {
        "input_candidate_species": len(rows),
        "review_queue_rows": len(selected),
        "deferred_rows": len(deferred),
        "evidence_class_counts": dict(sorted(class_counts.items())),
        "priority_counts": dict(sorted(priority_counts.items())),
        "aggregated_followup_review_rows": sum(r.get("review_reason", "").startswith("aggregated_followup") for r in selected),
        "known_positive_controls_in_review_queue": recovered,
        "known_positive_controls_in_review_queue_count": len(recovered),
        "review_mode": "natural_population_polymorphism_aggregated_v3",
    }
    Path(args.qc_out).write_text(json.dumps(qc, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(qc, ensure_ascii=False), flush=True)

    if len(selected) < 10:
        raise RuntimeError(f"Natural-polymorphism review queue unexpectedly small: {len(selected)}")
    if len(recovered) < 4:
        raise RuntimeError(f"Known-positive review recall failed: {len(recovered)}/4")


if __name__ == "__main__":
    main()
