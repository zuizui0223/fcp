#!/usr/bin/env python3
"""Acquire direct GIFT v3.2 Life_form_1 records for the frozen FCP species frame.

This is an outcome-blind source acquisition step: it reads species names only and
never reads D.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "polymorphism_species_attributes_step4_lifeform_source_20260910"
OUT.mkdir(parents=True, exist_ok=True)
SPECIES_METRICS = ROOT / "results" / "polymorphism_directionality_step1_20260909" / "species_metrics.csv"

GIFT_VERSION = "3.2"
GIFT_TRAIT_ID = "2.3.1"
GIFT_TRAIT_NAME = "Life_form_1"
GIFT_ENDPOINT = "https://gift.uni-goettingen.de/api/extended/index3.2.php"
USER_AGENT = "fcp-polymorphism-step4-lifeform/1.0 (direct source-backed acquisition)"
ALLOWED = {"phanerophyte", "chamaephyte", "hemicryptophyte", "cryptophyte", "therophyte"}
EXPECTED_UNITS = "phanerophyte, chamaephyte, hemicryptophyte, cryptophyte, therophyte"

SOURCE_RUN_ID = 29338693288
SOURCE_ARTIFACT_ID = 8313108327
SOURCE_ARTIFACT_NAME = "gift-direct-traits-29338693288"
SOURCE_ARTIFACT_DIGEST = "sha256:ae660ca791bab115f21d8b2ce7c9fdf65446c66365b274b6e33a7af1ada9ccca"


def _text(value: Any) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    return " ".join(str(value).split())


def _flag(value: Any) -> bool:
    return _text(value).casefold() in {"1", "true", "yes"}


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _locate_one(root: Path, name: str) -> Path:
    matches = sorted(root.rglob(name))
    if len(matches) != 1:
        raise RuntimeError(f"expected exactly one {name}, found {len(matches)}: {matches[:5]}")
    return matches[0]


def _load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _fetch_json(url: str, retries: int = 5, timeout: float = 180.0) -> Any:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    last: Exception | None = None
    for attempt in range(retries + 1):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as response:
                return json.loads(response.read().decode("utf-8"))
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            last = exc
            if attempt >= retries:
                raise
            time.sleep(min(30.0, 2.0 ** attempt))
    raise RuntimeError(str(last))


def _query_url(ref_id: str) -> str:
    return GIFT_ENDPOINT + "?" + urllib.parse.urlencode(
        {
            "query": "traits_raw",
            "traitid": GIFT_TRAIT_ID,
            "deriv": 0,
            "biasderiv": 0,
            "refid": ref_id,
        }
    )


def _split_tokens(raw: str) -> set[str]:
    return {x.strip().casefold() for x in raw.split("/") if x.strip()}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gift-artifact-root", type=Path, required=True)
    parser.add_argument("--max-workers", type=int, default=6)
    args = parser.parse_args()

    # Outcome firewall: only species names are read from the frozen D table.
    fcp = pd.read_csv(SPECIES_METRICS, usecols=["species"])
    fcp["species"] = fcp["species"].astype(str).str.strip()
    if len(fcp) != 369 or fcp["species"].nunique() != 369:
        raise RuntimeError("frozen FCP species fingerprint mismatch")
    fcp_names = set(fcp["species"])

    root = args.gift_artifact_root
    versions_path = _locate_one(root, "versions.json")
    traits_meta_path = _locate_one(root, "traits_meta.json")
    reference_traits_path = _locate_one(root, "reference_traits.json")
    references_path = _locate_one(root, "references.json")
    species_path = _locate_one(root, "species.csv.gz")

    versions = _load_json(versions_path)
    traits_meta = _load_json(traits_meta_path)
    reference_traits = _load_json(reference_traits_path)
    references = _load_json(references_path)
    published = {_text(row.get("version")) for row in versions}
    if GIFT_VERSION not in published:
        raise RuntimeError(f"GIFT {GIFT_VERSION} absent from pinned versions metadata")

    meta_hits = [row for row in traits_meta if _text(row.get("Lvl3")) == GIFT_TRAIT_ID]
    if len(meta_hits) != 1:
        raise RuntimeError(f"expected one metadata row for {GIFT_TRAIT_ID}, found {len(meta_hits)}")
    meta = meta_hits[0]
    if _text(meta.get("Trait2")) != GIFT_TRAIT_NAME:
        raise RuntimeError(f"trait name mismatch: {meta.get('Trait2')}")
    if _text(meta.get("type")).casefold() != "categorical":
        raise RuntimeError(f"trait type mismatch: {meta.get('type')}")
    if _text(meta.get("Units")).casefold() != EXPECTED_UNITS.casefold():
        raise RuntimeError(f"trait ontology mismatch: {meta.get('Units')}")

    refs_by_id = {_text(row.get("ref_ID")): row for row in references if _text(row.get("ref_ID"))}
    public_refs = {
        rid
        for rid, row in refs_by_id.items()
        if _text(row.get("restricted")) == "0" and _text(row.get("ref_long"))
    }
    pairs: set[str] = set()
    all_pairs: set[str] = set()
    for row in reference_traits:
        if _text(row.get("bias")) not in {"", "0"}:
            continue
        rid = _text(row.get("ref_ID"))
        if not rid:
            continue
        trait_values = {
            _text(v) for k, v in row.items() if str(k).startswith("trait") and _text(v)
        }
        if GIFT_TRAIT_ID in trait_values:
            all_pairs.add(rid)
            if rid in public_refs:
                pairs.add(rid)
    if not pairs:
        raise RuntimeError("no public unbiased GIFT references for Life_form_1")

    gift_species = pd.read_csv(species_path, usecols=["work_ID", "work_species"])
    work_to_species = {
        _text(row.work_ID): _text(row.work_species)
        for row in gift_species.itertuples(index=False)
        if _text(row.work_ID) and len(_text(row.work_species).split()) == 2
    }

    pair_list = sorted(pairs, key=lambda x: int(x))
    urls = {rid: _query_url(rid) for rid in pair_list}

    def fetch_pair(rid: str) -> tuple[str, list[dict[str, Any]], str | None]:
        try:
            payload = _fetch_json(urls[rid])
            if not isinstance(payload, list):
                raise RuntimeError("response is not a list")
            return rid, payload, None
        except Exception as exc:
            return rid, [], f"{type(exc).__name__}:{exc}"

    fetched: list[tuple[str, list[dict[str, Any]], str | None]] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, args.max_workers)) as ex:
        for item in ex.map(fetch_pair, pair_list):
            fetched.append(item)

    failures = {rid: err for rid, _rows, err in fetched if err}
    if failures:
        # Source acquisition is all-or-nothing so missing references cannot create silent coverage bias.
        _write_json(OUT / "reference_fetch_failures.json", failures)
        raise RuntimeError(f"failed to fetch {len(failures)} / {len(pair_list)} public Life_form_1 references")

    admitted: list[dict[str, Any]] = []
    audit = Counter()
    response_counts: dict[str, int] = {}
    for requested_ref, rows, _err in fetched:
        response_counts[requested_ref] = len(rows)
        reference = refs_by_id[requested_ref]
        citation = _text(reference.get("ref_long"))
        for row in rows:
            if _text(row.get("ref_ID")) != requested_ref or _text(row.get("trait_ID")) != GIFT_TRAIT_ID:
                audit["query_scope_mismatch"] += 1
                continue
            if _flag(row.get("derived")) or _flag(row.get("bias_deriv")):
                audit["derived_or_bias_derived"] += 1
                continue
            if not (_flag(row.get("matched")) and _flag(row.get("resolved"))):
                audit["unmatched_or_unresolved"] += 1
                continue
            if any(_flag(row.get(c)) for c in ("cf_genus", "cf_species", "aff_species")):
                audit["uncertain_source_name"] += 1
                continue
            if _text(row.get("subtaxon")):
                audit["infraspecific_scope"] += 1
                continue
            work_id = _text(row.get("work_ID"))
            species = work_to_species.get(work_id, "")
            if not species:
                audit["missing_binomial_work_species"] += 1
                continue
            raw = _text(row.get("trait_value"))
            tokens = _split_tokens(raw)
            if not tokens or not tokens <= ALLOWED:
                audit["unsupported_trait_value"] += 1
                continue
            if species not in fcp_names:
                continue
            admitted.append(
                {
                    "species": species,
                    "work_ID": work_id,
                    "trait_derived_ID": _text(row.get("trait_derived_ID")),
                    "ref_ID": requested_ref,
                    "trait_ID": GIFT_TRAIT_ID,
                    "trait_value_raw": raw,
                    "tokens": "|".join(sorted(tokens)),
                    "source_citation": citation,
                    "source_url": urls[requested_ref],
                }
            )

    evidence = pd.DataFrame(
        admitted,
        columns=[
            "species", "work_ID", "trait_derived_ID", "ref_ID", "trait_ID",
            "trait_value_raw", "tokens", "source_citation", "source_url",
        ],
    )
    if not evidence.empty:
        evidence = evidence.drop_duplicates().sort_values(["species", "ref_ID", "trait_derived_ID"])
    evidence.to_csv(OUT / "life_form_direct_evidence_fcp.csv.gz", index=False, compression="gzip")

    species_rows: list[dict[str, Any]] = []
    grouped = evidence.groupby("species") if not evidence.empty else None
    for sp in sorted(fcp_names):
        union: set[str] = set()
        n_records = 0
        n_refs = 0
        if grouped is not None and sp in grouped.groups:
            sub = grouped.get_group(sp)
            n_records = len(sub)
            n_refs = sub["ref_ID"].nunique()
            for token_blob in sub["tokens"]:
                union.update(str(token_blob).split("|"))
        if len(union) == 1:
            life = next(iter(union))
        elif len(union) > 1:
            life = "mixed"
        else:
            life = pd.NA
        species_rows.append(
            {
                "species": sp,
                "life_form": life,
                "direct_category_union": "|".join(sorted(union)),
                "n_direct_records": n_records,
                "n_public_references": n_refs,
            }
        )
    life = pd.DataFrame(species_rows)
    life.to_csv(OUT / "life_form_preoutcome.csv", index=False)

    counts = {str(k): int(v) for k, v in life["life_form"].value_counts().items()}
    n_available = int(life["life_form"].notna().sum())
    eligible_levels = sorted(str(k) for k, v in life["life_form"].value_counts().items() if int(v) >= 10)
    gate_total = n_available >= 50
    gate_levels = len(eligible_levels) >= 2
    result = {
        "analysis": "polymorphism_species_attributes_step4_lifeform_source",
        "protocol": "docs/POLYMORPHISM_SPECIES_ATTRIBUTES_STEP4_LIFEFORM_SOURCE_PROTOCOL_20260910.md",
        "date_jst": "2026-09-10",
        "D_association_computed": False,
        "new_image_acquisition": False,
        "fcp_species": len(fcp_names),
        "gift_version": GIFT_VERSION,
        "trait_id": GIFT_TRAIT_ID,
        "trait_name": GIFT_TRAIT_NAME,
        "trait_metadata": meta,
        "source_run_id": SOURCE_RUN_ID,
        "source_artifact_id": SOURCE_ARTIFACT_ID,
        "source_artifact_name": SOURCE_ARTIFACT_NAME,
        "source_artifact_digest": SOURCE_ARTIFACT_DIGEST,
        "metadata_sha256": {
            "versions": _sha256(versions_path),
            "traits_meta": _sha256(traits_meta_path),
            "reference_traits": _sha256(reference_traits_path),
            "references": _sha256(references_path),
            "species": _sha256(species_path),
        },
        "n_unbiased_references_with_trait": len(all_pairs),
        "n_public_unbiased_references_fetched": len(pair_list),
        "raw_response_rows": int(sum(response_counts.values())),
        "raw_response_counts_by_ref": response_counts,
        "admission_audit": dict(sorted(audit.items())),
        "fcp_admitted_evidence_rows": int(len(evidence)),
        "fcp_species_with_life_form": n_available,
        "life_form_counts": counts,
        "eligible_levels_n_ge_10": eligible_levels,
        "gate_total_ge_50": gate_total,
        "gate_at_least_two_levels_ge_10": gate_levels,
        "life_form_gate_pass": bool(gate_total and gate_levels),
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    _write_json(OUT / "result.json", result)

    md = [
        "# Polymorphism species attributes — Step 4 life-form source result",
        "",
        "**Outcome firewall: no D/covariate association was computed.**",
        "",
        f"- FCP species: **{len(fcp_names)}**",
        f"- GIFT trait: **{GIFT_TRAIT_ID} / {GIFT_TRAIT_NAME}**",
        f"- public unbiased GIFT references fetched: **{len(pair_list)}**",
        f"- admitted direct evidence rows for FCP species: **{len(evidence)}**",
        f"- FCP species with direct life form: **{n_available} / {len(fcp_names)}**",
        f"- standardized counts: `{counts}`",
        f"- categories eligible at n>=10: `{eligible_levels}`",
        f"- prospective life-form gate pass: **{bool(gate_total and gate_levels)}**",
        "",
        "No alternate life-form source or finer GIFT life-form trait will be substituted after D is opened.",
        "",
    ]
    (OUT / "RESULT.md").write_text("\n".join(md), encoding="utf-8")


if __name__ == "__main__":
    main()
