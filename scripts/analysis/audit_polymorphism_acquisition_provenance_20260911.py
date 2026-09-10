#!/usr/bin/env python3
"""Fail-closed acquisition-provenance audit for polymorphism paper v0.1.

This audit links the current reserve constructor to the exact historical frozen
candidate pool and to the pre-pixel acquisition contract that created that pool.
It is deliberately metadata/code only: no image pixels or flower-colour outcomes
are opened.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
SOURCE_COMMIT = "29584f3ad7ae0cd99a1d8f43459252af38f615da"
CONTRACT_PATH = "docs/supporting/global_monte_carlo_candidate_acquisition_contract_v1.json"
MANIFEST_PATH = "docs/supporting/global_monte_carlo_candidate_acquisition_manifest_v1.json"
QUERY_MODULE = "fcp_pipeline/global_candidate_acquisition.py"
PARSER_MODULE = "fcp_pipeline/random_photo_h9_pool.py"
PRODUCER = "scripts/acquisition/run_global_monte_carlo_candidate_acquisition_shard.py"
RESERVE_MODULE = ROOT / "fcp_pipeline/rgfca_reserve_replication.py"
OUT_DIR = ROOT / "results/polymorphism_acquisition_provenance_20260911"
EXPECTED_CANDIDATE_SHA256 = "f1319461d8883f3094cec8ac0e5fc247ff902b464f34146af275575b94edc9d2"
EXPECTED_QUERY = {
    "quality_grade": "research",
    "photos": True,
    "geo": True,
    "rank": "species",
    "flowering_term_id": 12,
    "flowering_term_value_id": 13,
    "maximum_positional_accuracy_m": 5000,
    "obscuration": "none",
    "allowed_photo_licenses": ["cc0", "cc-by", "cc-by-sa", "cc-by-nc", "cc-by-nc-sa"],
}


def git_bytes(path: str) -> bytes:
    return subprocess.check_output(["git", "show", f"{SOURCE_COMMIT}:{path}"], cwd=ROOT)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def literal_assignments(source: str) -> dict[str, Any]:
    tree = ast.parse(source)
    out: dict[str, Any] = {}
    for node in tree.body:
        if not isinstance(node, ast.Assign) or len(node.targets) != 1 or not isinstance(node.targets[0], ast.Name):
            continue
        try:
            out[node.targets[0].id] = ast.literal_eval(node.value)
        except Exception:
            pass
    return out


def stable_query_literal_keys(source: str) -> set[str]:
    tree = ast.parse(source)
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == "stable_candidate_query":
            for inner in ast.walk(node):
                if isinstance(inner, ast.Return) and isinstance(inner.value, ast.Dict):
                    keys: set[str] = set()
                    for key in inner.value.keys:
                        if isinstance(key, ast.Constant) and isinstance(key.value, str):
                            keys.add(key.value)
                    return keys
    raise RuntimeError("stable_candidate_query return dictionary not found")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def build_result() -> dict[str, Any]:
    contract_raw = git_bytes(CONTRACT_PATH)
    manifest_raw = git_bytes(MANIFEST_PATH)
    query_raw = git_bytes(QUERY_MODULE)
    parser_raw = git_bytes(PARSER_MODULE)
    producer_raw = git_bytes(PRODUCER)
    contract = json.loads(contract_raw)
    manifest = json.loads(manifest_raw)

    require(contract["protocol"] == "global-monte-carlo-candidate-acquisition-v1", "candidate contract protocol drift")
    require(manifest["protocol"] == contract["protocol"], "manifest/contract protocol mismatch")
    require(manifest["status"] == "complete_global_candidate_acquisition_premeasurement_gate_passed", "candidate premeasurement gate not passed")
    require(manifest["candidate_rows"] == 100000 and manifest["full_target_species"] == 1000, "candidate pool denominator drift")
    require(manifest["capacity_selected_raw_photo_target"] == 100, "candidate raw-photo target drift")
    require(manifest["candidate_image_pixels_opened"] is False and manifest["flower_colour_used"] is False, "candidate outcome firewall drift")
    require(manifest["lineage"]["candidate_contract_sha256"] == sha256(contract_raw), "manifest contract hash mismatch")
    require(manifest["lineage"]["candidate_photos_sha256"] == EXPECTED_CANDIDATE_SHA256, "candidate-pool hash drift")

    q = contract["query"]
    for key, value in EXPECTED_QUERY.items():
        require(q.get(key) == value, f"frozen query field drift: {key}")
    require(contract["selection"]["observer_cap_per_species"] == 2, "observer cap drift")
    require(contract["selection"]["final_photo_selection"].startswith("deterministic geographic maximin"), "candidate final selection drift")
    require(contract["outcome_firewall"]["flower_colour_used_for_page_selection"] is False, "colour used in page selection")
    require(contract["outcome_firewall"]["flower_colour_used_for_species_selection"] is False, "colour used in species selection")

    query_source = query_raw.decode("utf-8")
    query_keys = stable_query_literal_keys(query_source)
    required_runtime_keys = {"taxon_id", "quality_grade", "photos", "geo", "rank", "term_id", "term_value_id", "acc_below", "obscuration", "photo_license", "order_by", "order", "per_page", "page"}
    require(required_runtime_keys <= query_keys, "historical stable_candidate_query is missing frozen runtime keys")
    require(not ({"captive", "wild", "native"} & query_keys), "unexpected captive/wild/native query key present")

    parser_source = parser_raw.decode("utf-8")
    for token in ('quality_grade") or "") != "research"', 'maximum_positional_accuracy_m', 'obscured', 'geoprivacy'):
        require(token in parser_source, f"historical local parser guard missing: {token}")
    producer_source = producer_raw.decode("utf-8")
    require("stable_candidate_query(" in producer_source and "parse_h9_observation(" in producer_source, "historical producer is not wired to frozen query/parser")
    require('"candidate_image_pixels_opened": False' in producer_source and '"flower_colour_used": False' in producer_source, "historical producer firewall fields missing")

    reserve_source = RESERVE_MODULE.read_text(encoding="utf-8")
    reserve_assign = literal_assignments(reserve_source)
    require(reserve_assign.get("SOURCE_COMMIT") == SOURCE_COMMIT, "reserve source commit drift")
    require(reserve_assign.get("CANDIDATE") == "data/frozen/global_monte_carlo_candidate_photos_v1.csv", "reserve candidate path drift")
    require(EXPECTED_CANDIDATE_SHA256 in reserve_source, "reserve candidate hash guard drift")
    require("~candidate.inat_taxon_id.isin(selected.inat_taxon_id)" in reserve_source, "reserve complement rule missing")
    require("target_photos_per_species=target, maximum_species=budget, seed=20260918" in reserve_source, "reserve discovery-budget replay rule drift")
    require("No colour outcome is consulted in selection" in reserve_source, "reserve outcome-blind declaration missing")

    absent_explicit_filters = [name for name in ("native_range", "captive_false", "wild_true")]
    result = {
        "protocol": "polymorphism-acquisition-provenance-audit-v1",
        "status": "pass",
        "source_commit": SOURCE_COMMIT,
        "historical_inputs": {
            "candidate_contract": CONTRACT_PATH,
            "candidate_contract_sha256": sha256(contract_raw),
            "candidate_manifest": MANIFEST_PATH,
            "candidate_manifest_sha256": sha256(manifest_raw),
            "candidate_pool": manifest["files"]["candidate_photo_pool"],
            "candidate_pool_sha256": EXPECTED_CANDIDATE_SHA256,
            "candidate_rows": manifest["candidate_rows"],
            "candidate_species": manifest["full_target_species"],
            "photos_per_species": manifest["capacity_selected_raw_photo_target"],
        },
        "proven_acquisition_conditions": {
            "source": "iNaturalist v1 observations via InaturalistObservationClient",
            "quality_grade": q["quality_grade"],
            "photos_required": q["photos"],
            "georeference_required": q["geo"],
            "taxon_rank": q["rank"],
            "flowering_annotation_term_id": q["flowering_term_id"],
            "flowering_annotation_value_id": q["flowering_term_value_id"],
            "maximum_positional_accuracy_m": q["maximum_positional_accuracy_m"],
            "obscuration": q["obscuration"],
            "allowed_photo_licenses": q["allowed_photo_licenses"],
            "observer_cap_per_species": contract["selection"]["observer_cap_per_species"],
            "final_photo_selection": contract["selection"]["final_photo_selection"],
            "flower_colour_used_for_page_selection": False,
            "flower_colour_used_for_species_selection": False,
        },
        "reserve_lineage": {
            "construction": "entire species-disjoint complement of the original fixed hash-ranked 500-taxon measurement budget from the same frozen 1000-species candidate pool",
            "candidate_pool_sha256_guarded": True,
            "raw_photos_per_species": 100,
            "colour_outcome_used_for_reserve_selection": False,
        },
        "explicit_filters_not_applied_by_acquisition_code": absent_explicit_filters,
        "claim_boundary": "The geographic estimand describes organization among the observed community-photograph records. The acquisition code did not impose a native-range filter or explicit captive=false/wild=true parameter, so native-range-only or exclusively wild-population interpretation is not warranted.",
        "candidate_pixels_opened_by_audit": False,
        "flower_colour_outcomes_opened_by_audit": False,
    }
    return result


def render_md(result: dict[str, Any]) -> str:
    p = result["proven_acquisition_conditions"]
    return "\n".join([
        "# Polymorphism acquisition provenance audit",
        "",
        f"Status: **{result['status'].upper()}**",
        "",
        "## Proven lineage",
        "",
        f"- Historical source commit: `{result['source_commit']}`.",
        f"- Frozen candidate pool: `{result['historical_inputs']['candidate_pool']}`; SHA-256 `{result['historical_inputs']['candidate_pool_sha256']}`.",
        f"- Candidate denominator: {result['historical_inputs']['candidate_species']} species × {result['historical_inputs']['photos_per_species']} photos = {result['historical_inputs']['candidate_rows']:,} rows.",
        "- The reserve is the species-disjoint complement of the original fixed hash-ranked 500-taxon discovery budget within that same frozen candidate pool.",
        "",
        "## Proven acquisition conditions",
        "",
        f"- iNaturalist v1 observation source; `quality_grade={p['quality_grade']}`; photographs and georeferences required; taxon `rank={p['taxon_rank']}`.",
        f"- Flowering annotation gate: term {p['flowering_annotation_term_id']}, value {p['flowering_annotation_value_id']}.",
        f"- Positional accuracy ≤ {p['maximum_positional_accuracy_m']} m; `obscuration={p['obscuration']}`.",
        f"- Allowed photo licences: {', '.join(p['allowed_photo_licenses'])}.",
        f"- Observer cap: {p['observer_cap_per_species']} photos per species per observer; final selection: {p['final_photo_selection']}.",
        "- Page/species selection was outcome-blind to flower colour and candidate pixels remained unopened at acquisition.",
        "",
        "## Explicitly not established",
        "",
        "- No native-range filter was applied by the frozen acquisition query.",
        "- No explicit `captive=false` query parameter was applied.",
        "- No explicit `wild=true` query parameter was applied.",
        "",
        "Therefore Research Grade must not be rewritten as native-range-only or guaranteed wild-population sampling. The paper should interpret the spatial estimand as geographic organization in the observed community-photograph sample.",
        "",
    ])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    result = build_result()
    json_text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    md_text = render_md(result)
    expected = {OUT_DIR / "result.json": json_text, OUT_DIR / "RESULT.md": md_text}
    if args.check:
        missing_or_drifted = [str(path.relative_to(ROOT)) for path, text in expected.items() if not path.exists() or path.read_text(encoding="utf-8") != text]
        if missing_or_drifted:
            raise RuntimeError("provenance outputs missing or drifted: " + ", ".join(missing_or_drifted))
    else:
        OUT_DIR.mkdir(parents=True, exist_ok=True)
        for path, text in expected.items():
            path.write_text(text, encoding="utf-8", newline="\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
