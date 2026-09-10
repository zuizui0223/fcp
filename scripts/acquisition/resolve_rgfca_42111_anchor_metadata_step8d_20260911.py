#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request, urlopen

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
ALLOC = ROOT / "results/rgfca_42111_tiered_measurement_step8c_20260911/species_anchor_allocation.csv.gz"
SMOKE = ROOT / "results/rgfca_42111_anchor_resolution_step8d_20260911/result.json"
OUT = ROOT / "results/rgfca_42111_anchor_resolution_full_step8d_20260911"
API = "https://api.inaturalist.org/v1/observations"
USER_AGENT = "fcp-rgfca-42111-anchor-resolver/1.0 (github.com/zuizui0223/fcp)"
BATCH_SIZE = 200
REQUEST_INTERVAL_SECONDS = 1.05
N_SPECIES = 42111


def large_url(url: str) -> str:
    for size in ("square", "small", "medium"):
        token = f"/{size}."
        if token in url:
            return url.replace(token, "/large.")
    return url


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
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


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    smoke = json.loads(SMOKE.read_text(encoding="utf-8"))
    if smoke.get("status") != "complete_metadata_only_batch_resolution_smoke":
        raise RuntimeError("Step 8D smoke result is missing or incomplete")
    if smoke.get("requested_observations") != 200 or smoke.get("resolution_counts", {}).get("resolved") != 200:
        raise RuntimeError("Step 8D smoke did not resolve 200/200 exact anchors")
    if smoke.get("image_pixels_opened") is not False or smoke.get("flower_colour_used") is not False:
        raise RuntimeError("Step 8D smoke violated the metadata-only firewall")

    a = pd.read_csv(ALLOC)
    if len(a) != N_SPECIES or a["inat_taxon_id"].nunique() != N_SPECIES:
        raise RuntimeError("anchor allocation is not the exact 42,111-species universe")
    q = a.sort_values("breadth_rank", kind="mergesort").copy()
    if q["anchor_any_observation_id"].isna().any() or q["anchor_any_photo_id"].isna().any():
        raise RuntimeError("one or more species lack their frozen discovery anchor")

    expected = {
        int(row.anchor_any_observation_id): {
            "species": str(row.species),
            "inat_taxon_id": int(row.inat_taxon_id),
            "breadth_rank": int(row.breadth_rank),
            "observation_id": int(row.anchor_any_observation_id),
            "photo_id": int(row.anchor_any_photo_id),
            "anchor_source": str(row.anchor_any_source),
            "opened_existing_species": bool(row.opened_existing_species),
        }
        for row in q.itertuples(index=False)
    }
    if len(expected) != N_SPECIES:
        raise RuntimeError("breadth anchor observation IDs are not unique across species")

    rows: list[dict[str, object]] = []
    batch_receipts: list[dict[str, object]] = []
    ids = list(expected)
    last_started = 0.0
    started_all = time.time()
    for batch_index, start in enumerate(range(0, len(ids), BATCH_SIZE), start=1):
        batch = ids[start : start + BATCH_SIZE]
        elapsed = time.monotonic() - last_started
        if last_started and elapsed < REQUEST_INTERVAL_SECONDS:
            time.sleep(REQUEST_INTERVAL_SECONDS - elapsed)
        last_started = time.monotonic()
        t0 = time.time()
        results = fetch_batch(batch)
        by_id = {int(r.get("id")): r for r in results if isinstance(r, dict) and r.get("id") is not None}
        resolved = 0
        for oid in batch:
            base = expected[oid]
            r = by_id.get(oid)
            if r is None:
                rows.append({**base, "status": "observation_not_returned", "photo_url_large": "", "photo_license": "", "attribution": ""})
                continue
            photos = r.get("photos") or []
            match = next((p for p in photos if isinstance(p, dict) and p.get("id") is not None and int(p["id"]) == base["photo_id"]), None)
            if match is None:
                rows.append({**base, "status": "expected_photo_not_attached", "photo_url_large": "", "photo_license": "", "attribution": ""})
                continue
            purl = str(match.get("url") or "")
            if not purl:
                rows.append({**base, "status": "photo_url_missing", "photo_url_large": "", "photo_license": str(match.get("license_code") or ""), "attribution": str(match.get("attribution") or "")})
                continue
            rows.append({**base, "status": "resolved", "photo_url_large": large_url(purl), "photo_license": str(match.get("license_code") or ""), "attribution": str(match.get("attribution") or "")})
            resolved += 1
        batch_receipts.append({
            "batch_index": batch_index,
            "requested": len(batch),
            "api_results_returned": len(results),
            "resolved_exact_photo": resolved,
            "elapsed_seconds": float(time.time() - t0),
            "first_breadth_rank": int(expected[batch[0]]["breadth_rank"]),
            "last_breadth_rank": int(expected[batch[-1]]["breadth_rank"]),
        })
        print(json.dumps(batch_receipts[-1], sort_keys=True), flush=True)

    out = pd.DataFrame(rows).sort_values("breadth_rank", kind="mergesort").reset_index(drop=True)
    if len(out) != N_SPECIES or out["inat_taxon_id"].nunique() != N_SPECIES:
        raise RuntimeError("full resolver did not retain exactly 42,111 species")
    if out["breadth_rank"].tolist() != list(range(1, N_SPECIES + 1)):
        raise RuntimeError("breadth rank changed during resolution")

    csv_path = OUT / "anchor_metadata_42111.csv.gz"
    receipt_path = OUT / "batch_receipts.csv"
    out.to_csv(csv_path, index=False, compression="gzip")
    pd.DataFrame(batch_receipts).to_csv(receipt_path, index=False)
    counts = {str(k): int(v) for k, v in out["status"].value_counts().to_dict().items()}
    result = {
        "analysis": "rgfca_42111_anchor_resolution_full_step8d",
        "status": "complete_metadata_only_full_anchor_resolution",
        "species_universe": N_SPECIES,
        "batch_size": BATCH_SIZE,
        "requests": len(batch_receipts),
        "resolution_counts": counts,
        "resolved_species": int((out["status"] == "resolved").sum()),
        "unresolved_species": int((out["status"] != "resolved").sum()),
        "all_species_retained": True,
        "image_pixels_opened": False,
        "flower_colour_used": False,
        "elapsed_seconds": float(time.time() - started_all),
        "lineage": {
            "allocation_sha256": sha256_file(ALLOC),
            "smoke_result_sha256": sha256_file(SMOKE),
            "resolved_csv_sha256": sha256_file(csv_path),
            "batch_receipts_sha256": sha256_file(receipt_path),
        },
        "boundary": "Current API metadata resolution only. Unresolved frozen anchors remain in the 42,111-species denominator and are not substituted here.",
    }
    (OUT / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (OUT / "RESULT.md").write_text(
        "# RGFCA Step 8D — full 42,111-species anchor metadata resolution\n\n"
        f"- species retained: **{N_SPECIES:,} / {N_SPECIES:,}**\n"
        f"- API requests: **{len(batch_receipts)}** (<=200 exact IDs/request)\n"
        f"- exact anchors resolved: **{result['resolved_species']:,}**\n"
        f"- unresolved frozen anchors: **{result['unresolved_species']:,}**\n"
        f"- status counts: **{counts}**\n"
        f"- image pixels opened: **false**\n"
        f"- flower colour used: **false**\n\n"
        "Unresolved records are retained in the denominator and were not replaced.\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
