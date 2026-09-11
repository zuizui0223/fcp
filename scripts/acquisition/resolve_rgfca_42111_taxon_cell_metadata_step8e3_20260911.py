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
FRAME = ROOT / "results/rgfca_42111_taxon_cell_anchor_step8e2_20260911/taxon_cell_anchor_frame.csv.gz"
FRAME_RESULT = ROOT / "results/rgfca_42111_taxon_cell_anchor_step8e2_20260911/result.json"
SPECIES_META = ROOT / "results/rgfca_42111_anchor_resolution_full_step8d_20260911/anchor_metadata_42111.csv.gz"
OUT = ROOT / "results/rgfca_42111_taxon_cell_resolution_step8e3_20260911"
API = "https://api.inaturalist.org/v1/observations"
USER_AGENT = "fcp-rgfca-42111-taxon-cell-resolver/1.0 (github.com/zuizui0223/fcp)"
BATCH_SIZE = 200
REQUEST_INTERVAL_SECONDS = 1.05
EXPECTED_PAIRS = 85337


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def large_url(url: str) -> str:
    for size in ("square", "small", "medium"):
        token = f"/{size}."
        if token in url:
            return url.replace(token, "/large.")
    return url


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
    fr = json.loads(FRAME_RESULT.read_text(encoding="utf-8"))
    if fr.get("status") != "complete_metadata_only_unique_taxon_cell_anchor_freeze":
        raise RuntimeError("formal unique taxon-cell frame is not frozen")
    if int(fr.get("taxon_cell_anchors", -1)) != EXPECTED_PAIRS or int(fr.get("unresolved_pairs", -1)) != 0:
        raise RuntimeError("formal taxon-cell frame denominator is incomplete")
    if fr.get("image_pixels_opened") is not False or fr.get("flower_colour_used") is not False:
        raise RuntimeError("formal frame opened forbidden outcomes")

    frame = pd.read_csv(FRAME)
    if len(frame) != EXPECTED_PAIRS:
        raise RuntimeError("taxon-cell frame row count drift")
    if frame["observation_id"].nunique() != EXPECTED_PAIRS or frame["photo_id"].nunique() != EXPECTED_PAIRS:
        raise RuntimeError("taxon-cell frame identities are not globally unique")

    species_meta = pd.read_csv(SPECIES_META)
    species_meta = species_meta.loc[species_meta["status"].astype(str).eq("resolved")].copy()
    reuse_cols = ["observation_id", "photo_id", "photo_url_large", "photo_license", "attribution"]
    reusable = species_meta[reuse_cols].drop_duplicates(["observation_id", "photo_id"])
    out = frame.merge(reusable, on=["observation_id", "photo_id"], how="left", validate="one_to_one")
    reused_mask = out["photo_url_large"].fillna("").astype(str).ne("")
    out["resolution_status"] = "pending_api_resolution"
    out.loc[reused_mask, "resolution_status"] = "reused_species_anchor_resolution"

    pending = out.loc[~reused_mask, ["observation_id", "photo_id"]].copy()
    expected_photo = dict(zip(pending["observation_id"].astype(int), pending["photo_id"].astype(int)))
    pending_ids = list(expected_photo)
    if len(pending_ids) != pending["observation_id"].nunique():
        raise RuntimeError("pending observation IDs are not unique")

    resolved_rows: dict[int, dict[str, str]] = {}
    batch_receipts: list[dict[str, object]] = []
    last_started = 0.0
    started_all = time.time()
    for batch_index, start in enumerate(range(0, len(pending_ids), BATCH_SIZE), start=1):
        batch = pending_ids[start:start + BATCH_SIZE]
        elapsed = time.monotonic() - last_started
        if last_started and elapsed < REQUEST_INTERVAL_SECONDS:
            time.sleep(REQUEST_INTERVAL_SECONDS - elapsed)
        last_started = time.monotonic()
        t0 = time.time()
        results = fetch_batch(batch)
        by_id = {int(r.get("id")): r for r in results if isinstance(r, dict) and r.get("id") is not None}
        resolved = 0
        for oid in batch:
            r = by_id.get(int(oid))
            if r is None:
                resolved_rows[int(oid)] = {"resolution_status": "observation_not_returned", "photo_url_large": "", "photo_license": "", "attribution": ""}
                continue
            pid = int(expected_photo[int(oid)])
            photos = r.get("photos") or []
            match = next((p for p in photos if isinstance(p, dict) and p.get("id") is not None and int(p["id"]) == pid), None)
            if match is None:
                resolved_rows[int(oid)] = {"resolution_status": "expected_photo_not_attached", "photo_url_large": "", "photo_license": "", "attribution": ""}
                continue
            purl = str(match.get("url") or "")
            if not purl:
                resolved_rows[int(oid)] = {"resolution_status": "photo_url_missing", "photo_url_large": "", "photo_license": str(match.get("license_code") or ""), "attribution": str(match.get("attribution") or "")}
                continue
            resolved_rows[int(oid)] = {"resolution_status": "resolved_by_api", "photo_url_large": large_url(purl), "photo_license": str(match.get("license_code") or ""), "attribution": str(match.get("attribution") or "")}
            resolved += 1
        batch_receipts.append({
            "batch_index": batch_index,
            "requested": len(batch),
            "api_results_returned": len(results),
            "resolved_exact_photo": resolved,
            "elapsed_seconds": float(time.time() - t0),
        })
        print(json.dumps(batch_receipts[-1], sort_keys=True), flush=True)

    for idx in out.index[~reused_mask]:
        oid = int(out.at[idx, "observation_id"])
        r = resolved_rows[oid]
        out.at[idx, "resolution_status"] = r["resolution_status"]
        out.at[idx, "photo_url_large"] = r["photo_url_large"]
        out.at[idx, "photo_license"] = r["photo_license"]
        out.at[idx, "attribution"] = r["attribution"]

    out["resolved_exact_photo"] = out["resolution_status"].isin(["reused_species_anchor_resolution", "resolved_by_api"])
    if len(out) != EXPECTED_PAIRS or out[["inat_taxon_id", "cell_id"]].drop_duplicates().shape[0] != EXPECTED_PAIRS:
        raise RuntimeError("resolution changed the taxon-cell denominator")

    csv_path = OUT / "taxon_cell_anchor_metadata_85337.csv.gz"
    receipt_path = OUT / "batch_receipts.csv"
    out.to_csv(csv_path, index=False, compression="gzip", lineterminator="\n")
    pd.DataFrame(batch_receipts).to_csv(receipt_path, index=False, lineterminator="\n")
    counts = {str(k): int(v) for k, v in out["resolution_status"].value_counts().to_dict().items()}
    resolved_n = int(out["resolved_exact_photo"].sum())
    result = {
        "analysis": "rgfca_42111_taxon_cell_resolution_step8e3",
        "status": "complete_metadata_only_taxon_cell_resolution",
        "taxon_cell_pairs": EXPECTED_PAIRS,
        "species": int(out["inat_taxon_id"].nunique()),
        "occupied_cells": int(out["cell_id"].nunique()),
        "reused_species_anchor_rows": int(reused_mask.sum()),
        "api_resolution_rows_requested": int(len(pending_ids)),
        "api_requests": int(len(batch_receipts)),
        "resolved_exact_photo": resolved_n,
        "unresolved_rows": int(EXPECTED_PAIRS - resolved_n),
        "resolved_fraction": float(resolved_n / EXPECTED_PAIRS),
        "resolution_status_counts": counts,
        "image_pixels_opened": False,
        "flower_colour_used": False,
        "replacement_after_resolution_failure": False,
        "lineage": {
            "frame_sha256": sha256_file(FRAME),
            "species_metadata_sha256": sha256_file(SPECIES_META),
            "resolved_csv_sha256": sha256_file(csv_path),
            "batch_receipts_sha256": sha256_file(receipt_path),
        },
        "boundary": "Metadata transport only. Unresolved taxon-cell anchors remain unresolved and are not replaced."
    }
    (OUT / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (OUT / "RESULT.md").write_text(
        "# RGFCA Step 8E3 — taxon×cell anchor metadata resolution\n\n"
        f"- taxon×cell denominator: **{EXPECTED_PAIRS:,}**\n"
        f"- reused from species-anchor resolver: **{int(reused_mask.sum()):,}**\n"
        f"- additional exact IDs queried: **{len(pending_ids):,}**\n"
        f"- API requests: **{len(batch_receipts):,}**\n"
        f"- exact photo URLs resolved: **{resolved_n:,} / {EXPECTED_PAIRS:,}**\n"
        f"- unresolved: **{EXPECTED_PAIRS-resolved_n:,}**\n"
        f"- resolved fraction: **{resolved_n / EXPECTED_PAIRS:.6f}**\n"
        f"- status counts: **{counts}**\n"
        "- image pixels opened: **false**\n"
        "- flower colour used: **false**\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
