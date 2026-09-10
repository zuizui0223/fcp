#!/usr/bin/env python3
"""Build the Step 4 species-attribute covariate panel without opening D associations."""
from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import math
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "polymorphism_species_attributes_step4_preflight_20260910"
OUT.mkdir(parents=True, exist_ok=True)

SPECIES_METRICS = ROOT / "results" / "polymorphism_directionality_step1_20260909" / "species_metrics.csv"
RANGE_PRIMARY = ROOT / "data" / "frozen" / "global_monte_carlo_candidate_species_audit_v1.csv"
RANGE_SENSITIVITY = ROOT / "data" / "frozen" / "global_monte_carlo_capacity_scan_selected_species_v3.csv"
CANDIDATE_PHOTOS = ROOT / "data" / "frozen" / "global_monte_carlo_candidate_photos_v1.csv"

ISLAND_REPO = "zuizui0223/island"
ISLAND_RUN_ID = 29409415292
ISLAND_ARTIFACT_ID = 8340477381
ISLAND_ARTIFACT_NAME = "all-master-trait-ledger-29409415292"
ISLAND_ARTIFACT_DIGEST = "sha256:0cbba21f384d2b595272f13ed7be5a4b89d8ab0b77c6a0dd7a355fb3b3be1665"

INAT_TAXON_URL = "https://api.inaturalist.org/v1/taxa/{taxon_id}"
GIFT_VERSIONS_URL = "https://gift.uni-goettingen.de/api/index.php?query=versions"
GIFT_VERSION = "3.2"
GIFT_API_ROOT = "https://gift.uni-goettingen.de/api/extended"
USER_AGENT = "fcp-polymorphism-step4/1.0 (pre-outcome covariate freeze)"

POLLINATION_MAP = {
    "bees": "bee",
    "bumblebees": "bee",
    "flies": "other_animal",
    "birds": "other_animal",
    "moths": "other_animal",
    "butterflies": "other_animal",
    "wind": "wind",
    "mixed": "mixed",
    "self": "self",
}
LIFEFORM_SEARCH_TERMS = ("life", "growth", "habit", "wood")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")


def boolish(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return False
    return str(value).strip().lower() in {"true", "1", "yes"}


def fetch_json(url: str, *, retries: int = 5, timeout: float = 90.0) -> Any:
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
            time.sleep(min(20.0, 1.5 * (2 ** attempt)))
    raise RuntimeError(str(last))


def locate_one(root: Path, filename: str) -> Path:
    matches = sorted(root.rglob(filename))
    if len(matches) != 1:
        raise RuntimeError(f"expected exactly one {filename} below {root}, found {len(matches)}: {matches[:5]}")
    return matches[0]


def load_species_frame() -> pd.DataFrame:
    # Deliberately never read the D column in preflight.
    sp = pd.read_csv(SPECIES_METRICS, usecols=["species", "n"])
    sp["species"] = sp["species"].astype(str).str.strip()
    if len(sp) != 369 or sp["species"].nunique() != 369:
        raise RuntimeError(f"frozen species fingerprint mismatch: rows={len(sp)}, unique={sp['species'].nunique()}")
    if sp["species"].str.split().str.len().lt(2).any():
        bad = sp.loc[sp["species"].str.split().str.len().lt(2), "species"].tolist()
        raise RuntimeError(f"non-binomial species names in frozen frame: {bad[:10]}")
    sp = sp.rename(columns={"n": "n_classifiable"})
    sp["genus"] = sp["species"].str.split().str[0]
    return sp.sort_values("species").reset_index(drop=True)


def build_geo_covariates(species: pd.DataFrame) -> pd.DataFrame:
    rp = pd.read_csv(
        RANGE_PRIMARY,
        usecols=["species", "inat_taxon_id", "after_observer_cap", "maximum_span_km_after_observer_cap"],
    )
    rs = pd.read_csv(
        RANGE_SENSITIVITY,
        usecols=["species", "inat_taxon_id", "after_observer_cap", "maximum_span_km"],
    ).rename(
        columns={
            "inat_taxon_id": "inat_taxon_id_sensitivity",
            "after_observer_cap": "after_observer_cap_sensitivity",
        }
    )
    photos = pd.read_csv(CANDIDATE_PHOTOS, usecols=["species", "inat_taxon_id", "latitude"])
    photos["latitude"] = pd.to_numeric(photos["latitude"], errors="coerce")
    lat = (
        photos.dropna(subset=["latitude"])
        .groupby("species", as_index=False)
        .agg(
            latitude_n=("latitude", "size"),
            mean_latitude=("latitude", "mean"),
            min_latitude=("latitude", "min"),
            max_latitude=("latitude", "max"),
            inat_taxon_id_photos=("inat_taxon_id", "first"),
        )
    )
    lat["abs_mean_latitude"] = lat["mean_latitude"].abs()

    panel = species.merge(rp, on="species", how="left", validate="one_to_one")
    panel = panel.merge(rs, on="species", how="left", validate="one_to_one")
    panel = panel.merge(lat, on="species", how="left", validate="one_to_one")
    panel["log1p_span_primary"] = np.log1p(pd.to_numeric(panel["maximum_span_km_after_observer_cap"], errors="coerce"))
    panel["log1p_span_sensitivity"] = np.log1p(pd.to_numeric(panel["maximum_span_km"], errors="coerce"))

    ids = panel[["inat_taxon_id", "inat_taxon_id_sensitivity", "inat_taxon_id_photos"]].apply(
        pd.to_numeric, errors="coerce"
    )
    disagree = ids.dropna().nunique(axis=1).gt(1)
    panel["taxon_id_source_disagreement"] = disagree
    if disagree.any():
        rows = panel.loc[disagree, ["species", "inat_taxon_id", "inat_taxon_id_sensitivity", "inat_taxon_id_photos"]]
        rows.to_csv(OUT / "taxon_id_source_disagreements.csv", index=False)
        raise RuntimeError(f"taxon ID disagreement in {int(disagree.sum())} frozen species")
    return panel


def fetch_one_family(row: tuple[str, int]) -> dict[str, Any]:
    species, taxon_id = row
    url = INAT_TAXON_URL.format(taxon_id=int(taxon_id))
    stamp = datetime.now(timezone.utc).isoformat()
    try:
        payload = fetch_json(url, retries=4, timeout=60.0)
        results = payload.get("results", []) if isinstance(payload, dict) else []
        if len(results) != 1:
            return {
                "species": species,
                "inat_taxon_id": int(taxon_id),
                "family": None,
                "family_taxon_id": None,
                "retrieved_taxon_name": None,
                "retrieved_taxon_rank": None,
                "retrieval_status": f"unexpected_result_count_{len(results)}",
                "source_url": url,
                "retrieved_at_utc": stamp,
            }
        taxon = results[0]
        if int(taxon.get("id")) != int(taxon_id):
            raise RuntimeError(f"taxon id mismatch {taxon.get('id')} != {taxon_id}")
        lineage = list(taxon.get("ancestors") or []) + [taxon]
        fam = next((x for x in lineage if str(x.get("rank", "")).lower() == "family"), None)
        return {
            "species": species,
            "inat_taxon_id": int(taxon_id),
            "family": (fam or {}).get("name"),
            "family_taxon_id": (fam or {}).get("id"),
            "retrieved_taxon_name": taxon.get("name"),
            "retrieved_taxon_rank": taxon.get("rank"),
            "retrieval_status": "ok" if fam else "family_missing_in_lineage",
            "source_url": url,
            "retrieved_at_utc": stamp,
        }
    except Exception as exc:  # retain missingness rather than guessing
        return {
            "species": species,
            "inat_taxon_id": int(taxon_id),
            "family": None,
            "family_taxon_id": None,
            "retrieved_taxon_name": None,
            "retrieved_taxon_rank": None,
            "retrieval_status": f"error:{type(exc).__name__}:{exc}",
            "source_url": url,
            "retrieved_at_utc": stamp,
        }


def build_family_mapping(panel: pd.DataFrame) -> pd.DataFrame:
    rows = panel[["species", "inat_taxon_id"]].dropna().copy()
    rows["inat_taxon_id"] = rows["inat_taxon_id"].astype(int)
    tasks = list(rows.itertuples(index=False, name=None))
    out: list[dict[str, Any]] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
        for item in ex.map(fetch_one_family, tasks):
            out.append(item)
    fam = pd.DataFrame(out).sort_values("species").reset_index(drop=True)
    fam.to_csv(OUT / "family_taxonomy_receipt.csv", index=False)
    return fam


def build_pollination(species: pd.DataFrame, artifact_root: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    traits_path = locate_one(artifact_root, "all_species_traits.csv")
    evidence_path = locate_one(artifact_root, "all_species_trait_evidence.csv.gz")
    traits = pd.read_csv(
        traits_path,
        usecols=["species", "pollination_guild", "pollination_notes", "evidence_type", "confidence"],
    )
    ev = pd.read_csv(evidence_path)
    required = {"species", "field", "value", "source_kind", "source_url", "source_record_id", "evidence_type", "confidence", "source_backed"}
    missing = required - set(ev.columns)
    if missing:
        raise RuntimeError(f"island evidence ledger missing columns: {sorted(missing)}")
    ev = ev.loc[ev["field"].astype(str).eq("pollination_guild") & ev["source_backed"].map(boolish)].copy()
    if ev.empty:
        raise RuntimeError("no source-backed pollination evidence in island artifact")

    frame = species[["species"]].merge(traits, on="species", how="left", validate="one_to_one")
    exact = ev[["species", "value"]].drop_duplicates()
    exact["admissible_source_backed_exact"] = True
    frame = frame.merge(
        exact,
        left_on=["species", "pollination_guild"],
        right_on=["species", "value"],
        how="left",
        validate="one_to_one",
    )
    frame["admissible_source_backed_exact"] = frame["admissible_source_backed_exact"].fillna(False).astype(bool)
    frame["pollination_guild_raw"] = np.where(
        frame["admissible_source_backed_exact"] & frame["pollination_guild"].notna() & frame["pollination_guild"].ne("unknown"),
        frame["pollination_guild"],
        pd.NA,
    )
    frame["pollination_guild_mapped"] = frame["pollination_guild_raw"].map(POLLINATION_MAP)

    admissible_species = set(frame.loc[frame["pollination_guild_raw"].notna(), "species"])
    ev_fcp = ev.loc[ev["species"].isin(admissible_species)].copy()
    ev_fcp = ev_fcp.sort_values(["species", "value", "source_kind", "source_record_id"])
    ev_fcp.to_csv(OUT / "pollination_source_backed_evidence_fcp.csv.gz", index=False, compression="gzip")

    coverage = {
        "island_repo": ISLAND_REPO,
        "island_run_id": ISLAND_RUN_ID,
        "island_artifact_id": ISLAND_ARTIFACT_ID,
        "island_artifact_name": ISLAND_ARTIFACT_NAME,
        "island_recorded_artifact_digest": ISLAND_ARTIFACT_DIGEST,
        "traits_file_sha256": sha256(traits_path),
        "evidence_file_sha256": sha256(evidence_path),
        "fcp_species": int(len(species)),
        "species_with_nonunknown_final_pollination": int(frame["pollination_guild"].fillna("unknown").ne("unknown").sum()),
        "species_with_admissible_exact_source_backed_pollination": int(frame["pollination_guild_raw"].notna().sum()),
        "raw_counts": {str(k): int(v) for k, v in frame["pollination_guild_raw"].value_counts(dropna=False).items() if pd.notna(k)},
        "mapped_counts": {str(k): int(v) for k, v in frame["pollination_guild_mapped"].value_counts(dropna=False).items() if pd.notna(k)},
    }
    mapped_counts = frame["pollination_guild_mapped"].value_counts()
    eligible_levels = sorted(str(k) for k, v in mapped_counts.items() if int(v) >= 10)
    coverage["eligible_levels_n_ge_10"] = eligible_levels
    coverage["gate_total_ge_50"] = bool(frame["pollination_guild_raw"].notna().sum() >= 50)
    coverage["gate_at_least_two_levels_ge_10"] = bool(len(eligible_levels) >= 2)
    coverage["pollination_gate_pass"] = bool(coverage["gate_total_ge_50"] and coverage["gate_at_least_two_levels_ge_10"])

    keep = [
        "species", "pollination_guild", "pollination_guild_raw", "pollination_guild_mapped",
        "admissible_source_backed_exact", "pollination_notes", "evidence_type", "confidence",
    ]
    return frame[keep], ev_fcp, coverage


def scan_gift_lifeform_metadata() -> dict[str, Any]:
    versions = fetch_json(GIFT_VERSIONS_URL, retries=4, timeout=90.0)
    if not isinstance(versions, list):
        raise RuntimeError("GIFT versions response is not a list")
    published = {str(row.get("version", "")) for row in versions if isinstance(row, dict)}
    if GIFT_VERSION not in published:
        raise RuntimeError(f"GIFT {GIFT_VERSION} not present in published versions: {sorted(published)}")
    endpoint = f"{GIFT_API_ROOT}/index{GIFT_VERSION}.php"
    meta_url = endpoint + "?" + urllib.parse.urlencode({"query": "traits_meta"})
    reftraits_url = endpoint + "?" + urllib.parse.urlencode({"query": "reference_traits"})
    references_url = endpoint + "?" + urllib.parse.urlencode({"query": "references"})
    traits_meta = fetch_json(meta_url, retries=4, timeout=90.0)
    reference_traits = fetch_json(reftraits_url, retries=4, timeout=90.0)
    references = fetch_json(references_url, retries=4, timeout=90.0)
    if not all(isinstance(x, list) for x in [traits_meta, reference_traits, references]):
        raise RuntimeError("one or more GIFT metadata responses are not lists")

    candidates = []
    for row in traits_meta:
        blob = " ".join(str(v) for v in row.values()).casefold()
        matched = sorted({term for term in LIFEFORM_SEARCH_TERMS if term in blob})
        if matched:
            candidates.append({**row, "_matched_terms": "|".join(matched)})
    cand = pd.DataFrame(candidates)
    cand.to_csv(OUT / "gift_v3_2_lifeform_trait_metadata_candidates.csv", index=False)
    write_json(OUT / "gift_v3_2_versions.json", versions)

    public_refs = {
        str(row.get("ref_ID"))
        for row in references
        if isinstance(row, dict) and str(row.get("restricted", "")) == "0"
    }
    trait_ids: set[str] = set()
    for row in candidates:
        for key, value in row.items():
            if "trait" in str(key).casefold() and value not in (None, ""):
                text = str(value).strip()
                if any(ch.isdigit() for ch in text) and "." in text:
                    trait_ids.add(text)
    # Also retain explicit trait_ID-like values when present.
    for row in candidates:
        for key in ("trait_ID", "trait_id", "ID", "id"):
            if key in row and row[key] not in (None, ""):
                trait_ids.add(str(row[key]).strip())

    availability: list[dict[str, Any]] = []
    for trait_id in sorted(trait_ids):
        refs = set()
        public = set()
        for row in reference_traits:
            ref_id = str(row.get("ref_ID", ""))
            hit = any(str(v).strip() == trait_id for k, v in row.items() if str(k).startswith("trait"))
            if hit and ref_id:
                refs.add(ref_id)
                if ref_id in public_refs:
                    public.add(ref_id)
        availability.append(
            {
                "candidate_trait_id": trait_id,
                "n_reference_rows": len(refs),
                "n_public_reference_rows": len(public),
            }
        )
    pd.DataFrame(availability).to_csv(OUT / "gift_v3_2_lifeform_public_reference_availability.csv", index=False)

    return {
        "gift_version": GIFT_VERSION,
        "versions_url": GIFT_VERSIONS_URL,
        "endpoint": endpoint,
        "traits_meta_url": meta_url,
        "reference_traits_url": reftraits_url,
        "references_url": references_url,
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "n_traits_meta_rows": len(traits_meta),
        "n_candidate_metadata_rows": len(candidates),
        "candidate_trait_ids_detected": sorted(trait_ids),
        "availability": availability,
        "search_terms_frozen": list(LIFEFORM_SEARCH_TERMS),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--island-artifact-root", type=Path, required=True)
    args = parser.parse_args()

    species = load_species_frame()
    panel = build_geo_covariates(species)
    family = build_family_mapping(panel)
    panel = panel.merge(
        family[["species", "family", "family_taxon_id", "retrieval_status"]],
        on="species",
        how="left",
        validate="one_to_one",
    )

    poll, _ev, poll_cov = build_pollination(species, args.island_artifact_root)
    panel = panel.merge(poll, on="species", how="left", validate="one_to_one")
    life_meta = scan_gift_lifeform_metadata()

    # Outcome firewall: the pre-outcome panel must contain no D values.
    forbidden = {c for c in panel.columns if c == "D" or c.startswith("D_")}
    if forbidden:
        raise RuntimeError(f"outcome leaked into preflight panel: {sorted(forbidden)}")

    panel = panel.sort_values("species").reset_index(drop=True)
    panel.to_csv(OUT / "covariate_panel_preoutcome.csv", index=False)

    coverage = {
        "analysis": "polymorphism_species_attributes_step4_preflight",
        "protocol": "docs/POLYMORPHISM_SPECIES_ATTRIBUTES_STEP4_PREFLIGHT_PROTOCOL_20260910.md",
        "date_jst": "2026-09-10",
        "new_image_acquisition": False,
        "D_association_computed": False,
        "species_rows": int(len(panel)),
        "unique_species": int(panel["species"].nunique()),
        "genus_nonmissing": int(panel["genus"].notna().sum()),
        "family_nonmissing": int(panel["family"].notna().sum()),
        "family_retrieval_status": {str(k): int(v) for k, v in panel["retrieval_status"].fillna("missing").value_counts().items()},
        "span_primary_nonmissing": int(panel["log1p_span_primary"].notna().sum()),
        "span_sensitivity_nonmissing": int(panel["log1p_span_sensitivity"].notna().sum()),
        "latitude_nonmissing": int(panel["abs_mean_latitude"].notna().sum()),
        "latitude_photo_count": {
            "min": int(panel["latitude_n"].min()) if panel["latitude_n"].notna().any() else None,
            "median": float(panel["latitude_n"].median()) if panel["latitude_n"].notna().any() else None,
            "max": int(panel["latitude_n"].max()) if panel["latitude_n"].notna().any() else None,
        },
        "pollination": poll_cov,
        "life_form_source_preflight": life_meta,
        "input_sha256": {
            "species_metrics_source": sha256(SPECIES_METRICS),
            "range_primary": sha256(RANGE_PRIMARY),
            "range_sensitivity": sha256(RANGE_SENSITIVITY),
            "candidate_photos": sha256(CANDIDATE_PHOTOS),
        },
    }
    write_json(OUT / "result.json", coverage)

    md = [
        "# Polymorphism species attributes — Step 4 preflight result",
        "",
        "**Outcome firewall: no D/covariate association was computed.**",
        "",
        f"- frozen species: **{coverage['species_rows']}**",
        f"- family resolved: **{coverage['family_nonmissing']} / {coverage['species_rows']}**",
        f"- primary span available: **{coverage['span_primary_nonmissing']} / {coverage['species_rows']}**",
        f"- absolute latitude centroid available: **{coverage['latitude_nonmissing']} / {coverage['species_rows']}**",
        f"- admissible source-backed pollination: **{poll_cov['species_with_admissible_exact_source_backed_pollination']} / {coverage['species_rows']}**",
        f"- pollination gate pass: **{poll_cov['pollination_gate_pass']}**",
        f"- pollination mapped counts: `{poll_cov['mapped_counts']}`",
        f"- GIFT v3.2 life/growth/habit/wood metadata candidates: **{life_meta['n_candidate_metadata_rows']}**",
        f"- detected candidate trait IDs: `{life_meta['candidate_trait_ids_detected']}`",
        "",
        "The next step may freeze the exact life-form trait mapping and the exact taxonomic-clustering statistic, then open the predeclared D associations. No covariate may be added or backfilled based on those outcomes.",
        "",
    ]
    (OUT / "RESULT.md").write_text("\n".join(md), encoding="utf-8")


if __name__ == "__main__":
    main()
