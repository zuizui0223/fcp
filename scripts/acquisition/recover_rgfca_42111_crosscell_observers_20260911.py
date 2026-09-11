#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import time
from itertools import combinations
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
FRAME = ROOT / "results/rgfca_42111_taxon_cell_anchor_step8e2_20260911/taxon_cell_anchor_frame.csv.gz"
FRAME_RESULT = ROOT / "results/rgfca_42111_taxon_cell_anchor_step8e2_20260911/result.json"
PROTOCOL = ROOT / "docs/RGFCA_42111_CROSSCELL_DISCORDANCE_OBSERVER_PROTOCOL_20260911.md"
OUT = ROOT / "results/rgfca_42111_crosscell_observer_preflight_20260911"
API = "https://api.inaturalist.org/v1/observations"
USER_AGENT = "fcp-rgfca-42111-crosscell-observer/1.0 (github.com/zuizui0223/fcp)"
BATCH_SIZE = 200
REQUEST_INTERVAL_SECONDS = 1.05
EXPECTED_ROWS = 85337
SEED = 20260911


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch_batch(ids: list[int]) -> list[dict]:
    url = API + "/" + ",".join(map(str, ids)) + f"?per_page={len(ids)}"
    req = Request(url, headers={"Accept": "application/json", "User-Agent": USER_AGENT})
    try:
        with urlopen(req, timeout=90) as response:
            payload = json.loads(response.read().decode("utf-8"))
    except HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"iNaturalist HTTP {exc.code}: {body[:1000]}") from exc
    results = payload.get("results") or []
    if not isinstance(results, list):
        raise RuntimeError("iNaturalist returned non-list results")
    return results


def pair_hash(taxon_id: int, a: pd.Series, b: pd.Series) -> str:
    c1, c2 = sorted((int(a.cell_id), int(b.cell_id)))
    o1, o2 = sorted((int(a.observation_id), int(b.observation_id)))
    p1, p2 = sorted((int(a.photo_id), int(b.photo_id)))
    text = f"{SEED}|crosscell-pair|{taxon_id}|{c1}|{c2}|{o1}|{o2}|{p1}|{p2}"
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def cell_centroid(cell_id: int, n_lon: int = 18, n_sinlat: int = 9) -> tuple[float, float]:
    row = int(cell_id) // n_lon
    col = int(cell_id) % n_lon
    lon_width = 360.0 / n_lon
    lon = -180.0 + (col + 0.5) * lon_width
    s0 = -1.0 + 2.0 * row / n_sinlat
    s1 = -1.0 + 2.0 * (row + 1) / n_sinlat
    lat = math.degrees(math.asin((s0 + s1) / 2.0))
    return lat, lon


def centroid_distance_km(cell_a: int, cell_b: int) -> float:
    lat1, lon1 = cell_centroid(cell_a)
    lat2, lon2 = cell_centroid(cell_b)
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dl = math.radians(lon2 - lon1)
    dot = math.sin(p1) * math.sin(p2) + math.cos(p1) * math.cos(p2) * math.cos(dl)
    return math.acos(max(-1.0, min(1.0, dot))) * 6371.0088


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    if not PROTOCOL.exists():
        raise RuntimeError("frozen cross-cell protocol missing")
    fr = json.loads(FRAME_RESULT.read_text(encoding="utf-8"))
    if fr.get("status") != "complete_metadata_only_unique_taxon_cell_anchor_freeze":
        raise RuntimeError("taxon-cell frame not formally frozen")
    if fr.get("image_pixels_opened") is not False or fr.get("flower_colour_used") is not False:
        raise RuntimeError("taxon-cell frame already opened forbidden outcome")
    x = pd.read_csv(FRAME)
    if len(x) != EXPECTED_ROWS or x["observation_id"].nunique() != EXPECTED_ROWS or x["photo_id"].nunique() != EXPECTED_ROWS:
        raise RuntimeError("taxon-cell identity denominator drift")
    if x[["inat_taxon_id", "cell_id"]].drop_duplicates().shape[0] != EXPECTED_ROWS:
        raise RuntimeError("taxon-cell pair denominator drift")

    ids = x["observation_id"].astype(int).tolist()
    observer_by_obs: dict[int, str] = {}
    status_by_obs: dict[int, str] = {}
    batch_receipts: list[dict[str, object]] = []
    last_started = 0.0
    for batch_index, start in enumerate(range(0, len(ids), BATCH_SIZE), start=1):
        batch = ids[start:start+BATCH_SIZE]
        elapsed = time.monotonic() - last_started
        if last_started and elapsed < REQUEST_INTERVAL_SECONDS:
            time.sleep(REQUEST_INTERVAL_SECONDS - elapsed)
        last_started = time.monotonic()
        t0 = time.time()
        results = fetch_batch(batch)
        by_id = {int(r.get("id")): r for r in results if isinstance(r, dict) and r.get("id") is not None}
        recovered = 0
        for oid in batch:
            r = by_id.get(int(oid))
            if r is None:
                observer_by_obs[int(oid)] = ""
                status_by_obs[int(oid)] = "observation_not_returned"
                continue
            user = r.get("user") or {}
            uid = user.get("id") if isinstance(user, dict) else None
            if uid in (None, ""):
                observer_by_obs[int(oid)] = ""
                status_by_obs[int(oid)] = "observer_id_missing"
                continue
            observer_by_obs[int(oid)] = str(uid)
            status_by_obs[int(oid)] = "observer_resolved"
            recovered += 1
        batch_receipts.append({
            "batch_index": batch_index,
            "requested": len(batch),
            "api_results_returned": len(results),
            "observer_ids_resolved": recovered,
            "elapsed_seconds": float(time.time() - t0),
        })
        print(json.dumps(batch_receipts[-1], sort_keys=True), flush=True)

    x["observer_id"] = x["observation_id"].astype(int).map(observer_by_obs).fillna("")
    x["observer_resolution_status"] = x["observation_id"].astype(int).map(status_by_obs).fillna("unmapped")
    x["observer_resolved"] = x["observer_id"].astype(str).ne("")

    pair_rows: list[dict[str, object]] = []
    eligible_species = 0
    for (species, taxon_id), g in x.groupby(["species", "inat_taxon_id"], sort=False, observed=True):
        g = g.loc[g["observer_resolved"]].copy().reset_index(drop=True)
        choices: list[tuple[str, pd.Series, pd.Series]] = []
        for i, j in combinations(range(len(g)), 2):
            a = g.iloc[i]
            b = g.iloc[j]
            if int(a.cell_id) == int(b.cell_id):
                continue
            if str(a.observer_id) == str(b.observer_id):
                continue
            choices.append((pair_hash(int(taxon_id), a, b), a, b))
        if not choices:
            continue
        eligible_species += 1
        h, a, b = min(choices, key=lambda z: z[0])
        pair_rows.append({
            "species": str(species),
            "inat_taxon_id": int(taxon_id),
            "pair_hash": h,
            "cell_id_1": int(a.cell_id),
            "cell_id_2": int(b.cell_id),
            "observation_id_1": int(a.observation_id),
            "observation_id_2": int(b.observation_id),
            "photo_id_1": int(a.photo_id),
            "photo_id_2": int(b.photo_id),
            "observer_id_1": str(a.observer_id),
            "observer_id_2": str(b.observer_id),
            "cell_centroid_distance_km": centroid_distance_km(int(a.cell_id), int(b.cell_id)),
            "candidate_eligible_pairs": int(len(choices)),
        })

    pairs = pd.DataFrame(pair_rows)
    if len(pairs):
        pairs = pairs.sort_values(["inat_taxon_id", "pair_hash"], kind="mergesort").reset_index(drop=True)
        if pairs["inat_taxon_id"].nunique() != len(pairs):
            raise RuntimeError("selected pair frame is not one row per species")
        if (pairs["cell_id_1"] == pairs["cell_id_2"]).any() or (pairs["observer_id_1"] == pairs["observer_id_2"]).any():
            raise RuntimeError("selected pair violates distinct-cell or distinct-observer rule")

    enriched_path = OUT / "taxon_cell_observer_metadata_85337.csv.gz"
    pair_path = OUT / "crosscell_observer_disjoint_pairs.csv.gz"
    receipt_path = OUT / "observer_batch_receipts.csv"
    x.to_csv(enriched_path, index=False, compression="gzip", lineterminator="\n")
    pairs.to_csv(pair_path, index=False, compression="gzip", lineterminator="\n")
    pd.DataFrame(batch_receipts).to_csv(receipt_path, index=False, lineterminator="\n")

    n_resolved = int(x["observer_resolved"].sum())
    species_ge2cells = int((x.groupby("inat_taxon_id", observed=True)["cell_id"].nunique() >= 2).sum())
    result = {
        "analysis": "rgfca_42111_crosscell_observer_preflight",
        "status": "complete_metadata_only_observer_disjoint_pair_freeze",
        "taxon_cell_rows": EXPECTED_ROWS,
        "species_universe": int(x["inat_taxon_id"].nunique()),
        "species_with_ge2_frozen_cells": species_ge2cells,
        "observer_ids_resolved": n_resolved,
        "observer_ids_unresolved": int(EXPECTED_ROWS - n_resolved),
        "observer_resolution_fraction": float(n_resolved / EXPECTED_ROWS),
        "observer_api_requests": int(len(batch_receipts)),
        "observer_disjoint_crosscell_pair_species": int(len(pairs)),
        "pair_eligibility_fraction_of_ge2cell_species": float(len(pairs) / species_ge2cells) if species_ge2cells else None,
        "image_pixels_opened": False,
        "flower_colour_used": False,
        "replacement_after_observer_failure": False,
        "lineage": {
            "protocol_sha256": sha256_file(PROTOCOL),
            "taxon_cell_frame_sha256": sha256_file(FRAME),
            "observer_metadata_sha256": sha256_file(enriched_path),
            "pair_frame_sha256": sha256_file(pair_path),
            "batch_receipts_sha256": sha256_file(receipt_path),
        },
        "claim_boundary": "Pair frame identifies observer-disjoint cross-cell measurement opportunity only; no flower-colour discordance has been opened."
    }
    (OUT / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (OUT / "RESULT.md").write_text(
        "# RGFCA 42,111 — cross-cell observer-disjoint pair preflight\n\n"
        f"- taxon×cell rows: **{EXPECTED_ROWS:,}**\n"
        f"- species with >=2 frozen discovery cells: **{species_ge2cells:,}**\n"
        f"- observer IDs resolved: **{n_resolved:,} / {EXPECTED_ROWS:,} ({n_resolved/EXPECTED_ROWS:.3%})**\n"
        f"- frozen distinct-cell, distinct-observer pairs: **{len(pairs):,} species**\n"
        f"- pair eligibility among >=2-cell species: **{len(pairs)/species_ge2cells:.3%}**\n"
        "- image pixels opened: **false**\n"
        "- flower colour used: **false**\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
