#!/usr/bin/env python3
"""Validate the ITV/FCP systematic-search protocol and configuration."""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


def require(condition: bool, message: str, failures: list[str]) -> None:
    if not condition:
        failures.append(message)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--protocol", required=True)
    parser.add_argument("--script", required=True)
    args = parser.parse_args()

    config_path = Path(args.config)
    protocol_path = Path(args.protocol)
    script_path = Path(args.script)
    failures: list[str] = []

    for path in (config_path, protocol_path, script_path):
        require(path.is_file() and path.stat().st_size > 0, f"Missing or empty: {path}", failures)

    config = json.loads(config_path.read_text(encoding="utf-8"))
    protocol = protocol_path.read_text(encoding="utf-8")
    script = script_path.read_text(encoding="utf-8")

    queries = config.get("queries", [])
    require(len(queries) >= 12, "At least 12 complementary query blocks are required", failures)
    query_ids = [q.get("id") for q in queries]
    require(len(query_ids) == len(set(query_ids)), "Query IDs must be unique", failures)
    for query in queries:
        for field in ("id", "concept", "openalex", "crossref", "europe_pmc"):
            require(bool(str(query.get(field, "")).strip()), f"Query missing {field}: {query}", failures)

    seeds = config.get("seed_reviews", [])
    require(len(seeds) >= 5, "At least five prespecified seed reviews are required", failures)
    for seed in seeds:
        require(bool(re.match(r"^10\.\S+/\S+$", str(seed.get("doi", "")))), f"Invalid seed DOI: {seed}", failures)

    required_term_groups = {
        "display_colour", "variation", "within", "among", "natural", "cultivated",
        "induced", "ontogenetic", "non_display_floral", "continuous", "discrete",
    }
    terms = config.get("screening_terms", {})
    require(required_term_groups.issubset(terms), f"Missing term groups: {sorted(required_term_groups - set(terms))}", failures)
    for name in required_term_groups:
        require(len(terms.get(name, [])) >= 2, f"Term group too small: {name}", failures)

    required_protocol_phrases = [
        "systematic map", "ROSES", "PRISMA", "Intraspecific floral-colour variation",
        "Flower-colour polymorphism", "within_population", "among_population",
        "mixed_within_among", "cultivated_only_exclude", "duplicate independent review",
        "Cohen's kappa", "model results must not be consulted",
    ]
    for phrase in required_protocol_phrases:
        require(phrase.lower() in protocol.lower(), f"Protocol phrase missing: {phrase}", failures)

    required_script_phrases = [
        "openalex_search", "crossref_search", "europe_pmc_search", "citation_chase",
        "deduplicate", "human_review_required", "review_status", "adjudicated_decision",
    ]
    for phrase in required_script_phrases:
        require(phrase in script, f"Script capability missing: {phrase}", failures)

    report = {
        "status": "pass" if not failures else "fail",
        "query_blocks": len(queries),
        "seed_reviews": len(seeds),
        "term_groups": sorted(terms),
        "failures": failures,
    }
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
