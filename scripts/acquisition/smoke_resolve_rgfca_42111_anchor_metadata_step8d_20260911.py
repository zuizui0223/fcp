#!/usr/bin/env python3
from __future__ import annotations

import json
import time
from pathlib import Path
from urllib.request import Request, urlopen

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
ALLOC = ROOT / "results/rgfca_42111_tiered_measurement_step8c_20260911/species_anchor_allocation.csv.gz"
OUT = ROOT / "results/rgfca_42111_anchor_resolution_step8d_20260911"
API = "https://api.inaturalist.org/v1/observations"
USER_AGENT = "fcp-rgfca-42111-anchor-resolver/1.0 (github.com/zuizui0223/fcp)"
BATCH_SIZE = 50
REQUEST_INTERVAL_SECONDS = 1.05


def large_url(url: str) -> str:
    for size in ("square", "small", "medium"):
        token = f"/{size}."
        if token in url:
            return url.replace(token, "/large.")
    return url


def fetch_batch(ids: list[int]) -> list[dict]:
    url = API + "/" + ",".join(map(str, ids))
    req = Request(url, headers={"Accept": "application/json", "User-Agent": USER_AGENT})
    with urlopen(req, timeout=60) as response:
        payload = json.loads(response.read().decode("utf-8"))
    results = payload.get("results") or []
    if not isinstance(results, list):
        raise RuntimeError("iNaturalist returned non-list results")
    return results


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    a = pd.read_csv(ALLOC)
    if len(a) != 42111:
        raise RuntimeError("anchor allocation is not 42,111 species")
    q = a.sort_values("breadth_rank", kind="mergesort").head(200).copy()
    ids = [int(x) for x in q["anchor_any_observation_id"]]

    started = time.time()
    results: list[dict] = []
    request_count = 0
    for start in range(0, len(ids), BATCH_SIZE):
        if request_count:
            time.sleep(REQUEST_INTERVAL_SECONDS)
        batch = ids[start : start + BATCH_SIZE]
        results.extend(fetch_batch(batch))
        request_count += 1

    by_id = {int(r.get("id")): r for r in results if isinstance(r, dict) and r.get("id") is not None}

    rows = []
    for row in q.itertuples(index=False):
        oid = int(row.anchor_any_observation_id)
        pid = int(row.anchor_any_photo_id)
        r = by_id.get(oid)
        if r is None:
            rows.append({"species": row.species, "inat_taxon_id": int(row.inat_taxon_id), "observation_id": oid, "photo_id": pid, "status": "observation_not_returned", "photo_url_large": "", "photo_license": "", "attribution": ""})
            continue
        photos = r.get("photos") or []
        match = None
        for p in photos:
            if isinstance(p, dict) and p.get("id") is not None and int(p["id"]) == pid:
                match = p
                break
        if match is None:
            rows.append({"species": row.species, "inat_taxon_id": int(row.inat_taxon_id), "observation_id": oid, "photo_id": pid, "status": "expected_photo_not_attached", "photo_url_large": "", "photo_license": "", "attribution": ""})
            continue
        purl = str(match.get("url") or "")
        lic = str(match.get("license_code") or "")
        status = "resolved" if purl else "photo_url_missing"
        rows.append({"species": row.species, "inat_taxon_id": int(row.inat_taxon_id), "observation_id": oid, "photo_id": pid, "status": status, "photo_url_large": large_url(purl) if purl else "", "photo_license": lic, "attribution": str(match.get("attribution") or "")})

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "smoke_resolution_200.csv", index=False)
    counts = {str(k): int(v) for k, v in out.status.value_counts().to_dict().items()}
    result = {
        "analysis": "rgfca_42111_anchor_resolution_step8d_smoke",
        "status": "complete_metadata_only_batch_resolution_smoke",
        "requested_observations": len(ids),
        "api_results_returned": len(results),
        "resolution_counts": counts,
        "resolved_fraction": float((out.status == "resolved").mean()),
        "request_seconds": float(time.time() - started),
        "batch_size": BATCH_SIZE,
        "request_count": request_count,
        "transport_recovery_from_run": 34495173572,
        "image_pixels_opened": False,
        "flower_colour_used": False,
        "full_resolution_authorized": True,
        "boundary": "The smoke gate tests batch metadata resolution only. Unresolved anchors are retained and may not be silently replaced after colour opening.",
    }
    (OUT / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (OUT / "RESULT.md").write_text(
        "# RGFCA Step 8D — anchor metadata smoke\n\n"
        f"- requested: **{len(ids)}**\n"
        f"- API results: **{len(results)}**\n"
        f"- transport: **{request_count} requests × <= {BATCH_SIZE} IDs**\n"
        f"- resolved: **{counts.get('resolved', 0)} / {len(ids)}**\n"
        f"- resolved fraction: **{result['resolved_fraction']:.4f}**\n"
        f"- request seconds: **{result['request_seconds']:.2f}**\n\n"
        "No image URL was fetched by this step.\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
