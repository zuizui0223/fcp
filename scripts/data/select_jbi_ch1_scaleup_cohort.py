#!/usr/bin/env python3
"""Select the prospective JBI Chapter 1 scale-up cohort.

Selection uses only the frozen literature ranking, a human evidence ledger and the
predeclared family cap. Candidate photographs, colour measurements, Stage-A effects,
Stage-B surfaces and environmental layers are not read by this program.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
from typing import Iterable


PROTOCOL = "jbi-ch1-scaleup-cohort-v1"

RANKING_REQUIRED = {
    "rank",
    "canonical_name",
    "family",
    "best_doi",
}
LEDGER_REQUIRED = {
    "rank",
    "canonical_name",
    "family",
    "best_doi",
    "decision",
    "evidence_class",
    "photo_state_risk",
    "decision_basis",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path: Path, required: set[str], label: str) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or [])
        missing = sorted(required - fields)
        if missing:
            raise ValueError(f"{label} missing required columns: {missing}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
    if not rows:
        raise ValueError(f"{label} is empty")
    return rows


def parse_rank(value: str, *, label: str) -> int:
    try:
        rank = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} has invalid rank {value!r}") from exc
    if rank <= 0:
        raise ValueError(f"{label} rank must be positive: {rank}")
    return rank


def ranking_index(rows: Iterable[dict[str, str]]) -> dict[str, dict[str, str]]:
    by_name: dict[str, dict[str, str]] = {}
    ranks: set[int] = set()
    for row in rows:
        name = row["canonical_name"]
        if not name:
            raise ValueError("ranking contains a blank canonical_name")
        rank = parse_rank(row["rank"], label=f"ranking {name}")
        if name in by_name:
            raise ValueError(f"ranking contains duplicate species: {name}")
        if rank in ranks:
            raise ValueError(f"ranking contains duplicate rank: {rank}")
        ranks.add(rank)
        by_name[name] = row
    return by_name


def validate_ledger(
    ledger_rows: list[dict[str, str]],
    ranking_by_name: dict[str, dict[str, str]],
    contract: dict,
) -> list[dict[str, object]]:
    allowed = set(contract["allowed_decisions"])
    evidence_allowed = set(contract["eligible_evidence_classes"])
    eligible_risk = str(contract["eligible_photo_state_risk"])
    development = set(contract["parent_result"]["development_species"])

    seen_names: set[str] = set()
    seen_ranks: set[int] = set()
    validated: list[dict[str, object]] = []

    for raw in ledger_rows:
        name = raw["canonical_name"]
        if not name:
            raise ValueError("ledger contains a blank canonical_name")
        if name in seen_names:
            raise ValueError(f"ledger contains duplicate species: {name}")
        seen_names.add(name)

        rank = parse_rank(raw["rank"], label=f"ledger {name}")
        if rank in seen_ranks:
            raise ValueError(f"ledger contains duplicate rank: {rank}")
        seen_ranks.add(rank)

        source = ranking_by_name.get(name)
        if source is None:
            raise ValueError(f"ledger species absent from frozen ranking: {name}")
        source_rank = parse_rank(source["rank"], label=f"ranking {name}")
        if rank != source_rank:
            raise ValueError(f"rank mismatch for {name}: ledger={rank}, ranking={source_rank}")
        if raw["family"] != source["family"]:
            raise ValueError(
                f"family mismatch for {name}: ledger={raw['family']!r}, ranking={source['family']!r}"
            )
        if raw["best_doi"] != source["best_doi"]:
            raise ValueError(
                f"DOI mismatch for {name}: ledger={raw['best_doi']!r}, ranking={source['best_doi']!r}"
            )

        decision = raw["decision"]
        if decision not in allowed:
            raise ValueError(f"invalid decision for {name}: {decision!r}")
        if name in development and decision == "eligible":
            raise ValueError(f"completed development species cannot re-enter scale-up: {name}")
        if decision == "eligible":
            if raw["evidence_class"] not in evidence_allowed:
                raise ValueError(
                    f"eligible species {name} has unsupported evidence_class={raw['evidence_class']!r}"
                )
            if raw["photo_state_risk"] != eligible_risk:
                raise ValueError(
                    f"eligible species {name} must have photo_state_risk={eligible_risk!r}"
                )
            if not raw["decision_basis"]:
                raise ValueError(f"eligible species {name} lacks a decision basis")

        validated.append({**raw, "rank_int": rank})

    return validated


def select_cohort(
    ranking_rows: list[dict[str, str]],
    ledger_rows: list[dict[str, str]],
    contract: dict,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    if contract.get("protocol") != PROTOCOL:
        raise ValueError(f"unexpected contract protocol: {contract.get('protocol')!r}")
    if contract.get("status") != "prospective_before_scaleup_photo_acquisition":
        raise ValueError("scale-up contract is not in its prospective pre-acquisition state")
    if any(
        contract.get(key) is not False
        for key in (
            "selection_used_candidate_images",
            "selection_used_flower_colour_measurements",
            "selection_used_stage_a_effects",
            "selection_used_stage_b_surfaces",
            "selection_used_environmental_layers",
        )
    ):
        raise ValueError("selection contract permits a forbidden post-colour input")

    ranking_by_name = ranking_index(ranking_rows)
    validated = validate_ledger(ledger_rows, ranking_by_name, contract)
    eligible = sorted(
        (row for row in validated if row["decision"] == "eligible"),
        key=lambda row: (int(row["rank_int"]), str(row["canonical_name"])),
    )

    target = int(contract["target_new_species"])
    family_cap = int(contract["maximum_species_per_family"])
    if target <= 0 or family_cap <= 0:
        raise ValueError("target_new_species and maximum_species_per_family must be positive")

    selected: list[dict[str, object]] = []
    skipped_for_family_cap: list[dict[str, object]] = []
    family_counts: Counter[str] = Counter()

    for row in eligible:
        family = str(row["family"])
        if family_counts[family] >= family_cap:
            skipped_for_family_cap.append(row)
            continue
        family_counts[family] += 1
        selected.append(
            {
                "cohort_order": len(selected) + 1,
                "rank": int(row["rank_int"]),
                "canonical_name": str(row["canonical_name"]),
                "family": family,
                "best_doi": str(row["best_doi"]),
                "evidence_class": str(row["evidence_class"]),
                "photo_state_risk": str(row["photo_state_risk"]),
                "decision_basis": str(row["decision_basis"]),
            }
        )
        if len(selected) == target:
            break

    if len(selected) != target:
        raise ValueError(
            f"only {len(selected)} eligible species satisfy the frozen family cap; target={target}"
        )
    if len({row["canonical_name"] for row in selected}) != target:
        raise ValueError("selected cohort contains duplicate species")
    if max(Counter(str(row["family"]) for row in selected).values()) > family_cap:
        raise ValueError("selected cohort violates the family cap")

    return selected, skipped_for_family_cap


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError("refusing to write an empty cohort")
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--ranking",
        type=Path,
        default=Path("data/global_flower_colour_species_ranked.csv"),
    )
    parser.add_argument(
        "--ledger",
        type=Path,
        default=Path("docs/supporting/jbi_ch1_scaleup_species_ledger_v1.csv"),
    )
    parser.add_argument(
        "--contract",
        type=Path,
        default=Path("docs/supporting/jbi_ch1_scaleup_contract_v1.json"),
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path("docs/supporting/jbi_ch1_scaleup_cohort_v1.csv"),
    )
    parser.add_argument(
        "--manifest-json",
        type=Path,
        default=Path("docs/supporting/jbi_ch1_scaleup_cohort_manifest_v1.json"),
    )
    args = parser.parse_args()

    ranking_rows = read_csv(args.ranking, RANKING_REQUIRED, "ranking")
    ledger_rows = read_csv(args.ledger, LEDGER_REQUIRED, "ledger")
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    selected, skipped_for_family_cap = select_cohort(ranking_rows, ledger_rows, contract)
    write_csv(args.output_csv, selected)

    decisions = Counter(row["decision"] for row in ledger_rows)
    manifest = {
        "protocol": PROTOCOL,
        "status": "scaleup_cohort_selected_before_photo_acquisition",
        "contract_sha256": sha256(args.contract),
        "ranking_sha256": sha256(args.ranking),
        "ledger_sha256": sha256(args.ledger),
        "cohort_csv_sha256": sha256(args.output_csv),
        "target_new_species": int(contract["target_new_species"]),
        "selected_new_species": len(selected),
        "selected_species": [row["canonical_name"] for row in selected],
        "selected_families": dict(sorted(Counter(row["family"] for row in selected).items())),
        "maximum_selected_rank": max(int(row["rank"]) for row in selected),
        "ledger_decision_counts": dict(sorted(decisions.items())),
        "skipped_for_family_cap": [row["canonical_name"] for row in skipped_for_family_cap],
        "selection_used_candidate_images": False,
        "selection_used_flower_colour_measurements": False,
        "selection_used_stage_a_effects": False,
        "selection_used_stage_b_surfaces": False,
        "selection_used_environmental_layers": False,
        "next_gate": "iNaturalist metadata-only feasibility and deterministic 200-photo acquisition per selected species",
    }
    args.manifest_json.parent.mkdir(parents=True, exist_ok=True)
    args.manifest_json.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
