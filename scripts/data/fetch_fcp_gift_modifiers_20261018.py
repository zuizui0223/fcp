#!/usr/bin/env python3
"""Freeze exact-name GIFT 3.2 trait matches for the full measured FCP roster.

Public API pages are cached outside the repository; matched source records and
page hashes are retained with the analytical outputs. No synonym guessing or
genus/family-level trait imputation is performed.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pandas as pd

API = "https://gift.uni-goettingen.de/api/extended/index3.2.php"
TRAITS = {"self_fertilization": "3.1.1", "lifecycle": "2.1.1", "woodiness": "1.1.1"}
PAGE_SIZE = 10000


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def get_page(cache: Path, name: str, params: dict) -> tuple[list, dict]:
    url = API + "?" + urlencode(params)
    path = cache / (name + ".json")
    if not path.exists():
        request = Request(url, headers={"User-Agent": "FCP-within-species-traits/1.0", "Accept": "application/json"})
        with urlopen(request, timeout=45) as response:
            raw = response.read()
        rows = json.loads(raw)
        if not isinstance(rows, list):
            raise ValueError(f"GIFT returned a non-list: {url}")
        path.write_bytes(raw)
    rows = json.loads(path.read_bytes())
    if not isinstance(rows, list):
        raise ValueError(f"Invalid cached GIFT page: {path}")
    return rows, {"url": url, "cache_file": path.name, "rows": len(rows), "bytes": path.stat().st_size, "sha256": digest(path)}


def get_series(cache: Path, prefix: str, params: dict) -> tuple[pd.DataFrame, list[dict]]:
    combined, receipts = [], []
    for start in range(0, 1000000, PAGE_SIZE):
        rows, receipt = get_page(cache, f"{prefix}_{start:06d}", {**params, "startat": start, "limit": PAGE_SIZE})
        combined.extend(rows)
        receipts.append(receipt)
        if len(rows) < PAGE_SIZE:
            break
    else:
        raise ValueError(f"Pagination did not terminate for {prefix}")
    frame = pd.DataFrame(combined)
    if frame.empty or "work_ID" not in frame or frame.work_ID.duplicated().any():
        raise ValueError(f"Empty or repeated species identifiers in {prefix}")
    frame["work_ID"] = pd.to_numeric(frame.work_ID, errors="raise").astype("int64")
    return frame, receipts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    args.cache.mkdir(parents=True, exist_ok=True)
    args.outdir.mkdir(parents=True, exist_ok=True)
    roster_path = root / "results/fcp_real_447500_white_axis_decomposition_20261018/species_photo_composition.csv"
    env_path = root / "results/fcp_real_447500_environment_description_20261018/primary_ge5_species_environment.csv"
    roster = pd.read_csv(roster_path)
    env = pd.read_csv(env_path)
    assert len(roster) == 4475 and roster.inat_taxon_id.nunique() == 4475
    assert env.inat_taxon_id.nunique() == 1574
    metadata, meta_receipt = get_page(args.cache, "traits_meta", {"query": "traits_meta"})
    selected_metadata = [r for r in metadata if r["Lvl3"] in TRAITS.values()]
    if len(selected_metadata) != len(TRAITS):
        raise ValueError("Required trait metadata missing")
    jobs = [("species", {"query": "species"})] + [("traits_" + trait_id.replace(".", "_"), {"query": "traits", "traitid": trait_id, "biasref": 0, "biasderiv": 0}) for trait_id in TRAITS.values()]
    with ThreadPoolExecutor(max_workers=4) as executor:
        fetched = list(executor.map(lambda job: get_series(args.cache, *job), jobs))
    species, species_receipts = fetched[0]
    duplicates = species.work_species.duplicated(keep=False)
    unique_names = species.loc[~duplicates, ["work_species", "work_ID", "work_genus", "work_author"]]
    audit = roster[["inat_taxon_id", "species", "source_cohort", "n_classifiable", "n_white", "n_nonwhite"]].merge(unique_names, how="left", left_on="species", right_on="work_species", validate="many_to_one")
    audit["environment_panel"] = audit.inat_taxon_id.isin(set(env.inat_taxon_id))
    ambiguous = set(species.loc[duplicates, "work_species"])
    audit["name_match_status"] = "unmatched_exact_name"
    audit.loc[audit.work_ID.notna(), "name_match_status"] = "unique_exact_name"
    audit.loc[audit.species.isin(ambiguous), "name_match_status"] = "ambiguous_exact_name"
    assert len(audit) == 4475 and audit.inat_taxon_id.nunique() == 4475
    source_records, receipts, coverage = [], [meta_receipt, *species_receipts], []
    for (trait, trait_id), (values, trait_receipts) in zip(TRAITS.items(), fetched[1:]):
        receipts.extend(trait_receipts)
        q = values.loc[values.work_ID.isin(set(audit.work_ID.dropna().astype(int)))].copy()
        q["trait"] = trait
        q["trait_ID"] = trait_id
        q["agreement"] = pd.to_numeric(q.agreement, errors="coerce")
        q["admitted_agreement"] = q.agreement.ge(0.66) & q.trait_value.notna() & q.trait_value.astype(str).str.len().gt(0)
        source_records.append(q)
        admitted = q.loc[q.admitted_agreement].set_index("work_ID")
        audit[trait] = audit.work_ID.map(admitted.trait_value)
        for group, panel in [("all4475", audit), ("environment1574", audit.loc[audit.environment_panel])]:
            for cohort, sub in [("pooled", panel), *list(panel.groupby("source_cohort"))]:
                coverage.append({"panel": group, "cohort": cohort, "trait": trait, "trait_ID": trait_id, "n_species": len(sub), "n_known": int(sub[trait].notna().sum()), "n_missing": int(sub[trait].isna().sum()), "value_counts": json.dumps(sub[trait].value_counts().sort_index().to_dict(), sort_keys=True)})
    audit.to_csv(args.outdir / "species_trait_join.csv", index=False)
    pd.concat(source_records, ignore_index=True).to_csv(args.outdir / "gift_matched_source_records.csv", index=False)
    pd.DataFrame(coverage).to_csv(args.outdir / "trait_coverage.csv", index=False)
    profile = {"completed_utc": datetime.now(timezone.utc).isoformat(), "gift_version": "3.2", "api": API, "roster_species": 4475, "environment_species": 1574, "gift_species_rows": len(species), "match_status_counts": audit.name_match_status.value_counts().to_dict(), "traits": selected_metadata, "agreement_threshold": 0.66, "biasref": 0, "biasderiv": 0, "missing_semantics": "unknown documentation, never biological absence", "matching": "unique exact species string only; no synonym/genus/family imputation", "self_fertilization_interpretation": "GIFT recorded present/absent category; not realized selfing rate or a pigment mechanism", "source_sha256": {str(roster_path.relative_to(root)): digest(roster_path), str(env_path.relative_to(root)): digest(env_path)}, "public_api_pages": receipts, "code_sha256": digest(Path(__file__)), "output_sha256": {p.name: digest(p) for p in sorted(args.outdir.glob("*.csv"))}}
    (args.outdir / "profile.json").write_text(json.dumps(profile, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({"name_matches": profile["match_status_counts"], "pooled_coverage": [r for r in coverage if r["cohort"] == "pooled"]}, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
