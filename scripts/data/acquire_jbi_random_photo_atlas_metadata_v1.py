#!/usr/bin/env python3
"""Acquire metadata only for the prospective species-unfixed random photo atlas v1."""
from __future__ import annotations
import argparse, csv, json, time
from pathlib import Path
from typing import Any
from urllib.parse import urlencode
from urllib.request import Request, urlopen

USER_AGENT = "zuizui0223-fcp-random-photo-atlas-v1/1.0 (research reproducibility)"
DEFAULT_BASE = "https://api.inaturalist.org/v1"
SPECIES_RANKS = {"species", "subspecies", "variety", "form", "hybrid"}


def get_json(url: str, pause: float, retries: int = 5) -> dict[str, Any]:
    last = None
    for attempt in range(retries):
        try:
            req = Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
            with urlopen(req, timeout=90) as r:
                out = json.load(r)
            if pause: time.sleep(pause)
            return out
        except Exception as exc:
            last = exc
            if attempt + 1 < retries: time.sleep(min(2 ** attempt, 16))
    raise RuntimeError(f"iNaturalist request failed: {last}")


def coords(obs: dict[str, Any]) -> tuple[float | None, float | None]:
    gj = obs.get("geojson") or {}
    xy = gj.get("coordinates") or []
    if isinstance(xy, list) and len(xy) >= 2:
        try: return float(xy[1]), float(xy[0])
        except Exception: pass
    return None, None


def large_url(url: str | None) -> str | None:
    if not url: return None
    for token in ("square", "small", "medium"):
        if f"/{token}." in url: return url.replace(f"/{token}.", "/large.")
    return url


def species_key(taxon: dict[str, Any]) -> str | None:
    rank = str(taxon.get("rank") or "").lower()
    name = str(taxon.get("name") or "").strip()
    parts = name.split()
    if rank not in SPECIES_RANKS or len(parts) < 2: return None
    return " ".join(parts[:2])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--qc-out", type=Path, required=True)
    ap.add_argument("--created-d1", default="2026-08-28")
    ap.add_argument("--created-d2", default="2026-09-02")
    ap.add_argument("--base-url", default=DEFAULT_BASE)
    ap.add_argument("--pause", type=float, default=0.25)
    ap.add_argument("--per-page", type=int, default=200)
    ap.add_argument("--max-eligible", type=int, default=100000)
    args = ap.parse_args()

    fixed = {
        "taxon_id": 47125,
        "quality_grade": "research",
        "photos": "true",
        "geo": "true",
        "captive": "false",
        "term_id": 12,
        "term_value_id": 13,
        "created_d1": args.created_d1,
        "created_d2": args.created_d2,
        "per_page": args.per_page,
        "order_by": "id",
        "order": "asc",
    }
    rows: list[dict[str, Any]] = []
    seen_obs: set[int] = set(); seen_photo: set[int] = set()
    id_above = 0; requests = 0; raw_results = 0
    rejected = {"obscured":0,"bad_coords":0,"bad_taxon":0,"no_photo":0,"duplicate":0}
    complete = False
    while True:
        params = dict(fixed)
        if id_above: params["id_above"] = id_above
        url = args.base_url.rstrip("/") + "/observations?" + urlencode(params)
        payload = get_json(url, args.pause); requests += 1
        results = payload.get("results", []) or []; raw_results += len(results)
        if not results:
            complete = True; break
        max_id = id_above
        for obs in results:
            try: oid = int(obs["id"])
            except Exception: continue
            max_id = max(max_id, oid)
            if oid in seen_obs:
                rejected["duplicate"] += 1; continue
            seen_obs.add(oid)
            if bool(obs.get("obscured")) or str(obs.get("geoprivacy") or "").lower() not in {"", "open"}:
                rejected["obscured"] += 1; continue
            lat, lon = coords(obs)
            if lat is None or lon is None or not (-90 <= lat <= 90 and -180 <= lon <= 180):
                rejected["bad_coords"] += 1; continue
            taxon = obs.get("taxon") or {}; skey = species_key(taxon)
            if not skey:
                rejected["bad_taxon"] += 1; continue
            photos = sorted((obs.get("photos") or []), key=lambda p: int(p.get("id") or 10**30))
            if not photos:
                rejected["no_photo"] += 1; continue
            photo = photos[0]
            try: pid = int(photo["id"])
            except Exception:
                rejected["no_photo"] += 1; continue
            if pid in seen_photo:
                rejected["duplicate"] += 1; continue
            seen_photo.add(pid)
            created = str(obs.get("created_at") or "")
            if not created.startswith((args.created_d1, args.created_d2)):
                # Dates between endpoints are accepted below; this guard only protects malformed fields.
                try:
                    day = created[:10]
                    if not (args.created_d1 <= day <= args.created_d2): continue
                except Exception: continue
            user = obs.get("user") or {}
            rows.append({
                "observation_id": oid, "photo_id": pid, "species_key": skey,
                "inat_taxon_id": taxon.get("id"), "inat_taxon_name": taxon.get("name"), "inat_rank": taxon.get("rank"),
                "latitude": f"{lat:.8f}", "longitude": f"{lon:.8f}",
                "positional_accuracy_m": obs.get("positional_accuracy"), "created_at": created,
                "observed_on": obs.get("observed_on"), "quality_grade": obs.get("quality_grade"),
                "photo_url": large_url(photo.get("url")), "photo_url_api": photo.get("url"),
                "photo_license": photo.get("license_code"), "attribution": photo.get("attribution"),
                "observer_id": user.get("id"), "observer": user.get("login"),
            })
            if len(rows) >= args.max_eligible:
                break
        if len(rows) >= args.max_eligible:
            break
        if max_id <= id_above:
            raise RuntimeError("iNaturalist id_above pagination failed to advance")
        id_above = max_id
        if len(results) < args.per_page:
            complete = True; break
        if requests % 25 == 0: print({"requests":requests,"eligible":len(rows),"last_id":id_above}, flush=True)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else ["observation_id","photo_id","species_key","latitude","longitude","photo_url"]
    with args.out.open("w", encoding="utf-8", newline="") as h:
        w=csv.DictWriter(h, fieldnames=fields, lineterminator="\n"); w.writeheader(); w.writerows(rows)
    qc = {
        "status": "complete" if complete else "not_evaluable_metadata_ceiling",
        "protocol": "jbi-random-photo-atlas-metadata-v1",
        "created_d1": args.created_d1, "created_d2": args.created_d2,
        "terminal_snapshot_exclusion": "source window begins after 2026-08-27 terminal Open Data snapshot",
        "filters": fixed, "requests": requests, "raw_results": raw_results,
        "eligible_observations": len(rows), "unique_species_keys": len({r['species_key'] for r in rows}),
        "unique_photos": len(seen_photo), "rejected": rejected,
        "metadata_complete": complete, "hard_ceiling": args.max_eligible,
        "pixel_accessed": False,
    }
    args.qc_out.write_text(json.dumps(qc, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    print(json.dumps(qc, indent=2, sort_keys=True)); return 0

if __name__ == "__main__": raise SystemExit(main())
