#!/usr/bin/env python3
from __future__ import annotations

import argparse
import io
import json
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests
from requests.adapters import HTTPAdapter

DEFAULT_VARS = ["tmin", "tmax", "ppt", "def", "vpd"]
DEFAULT_BASE = "http://thredds.northwestknowledge.net:8080/thredds/ncss/agg_terraclimate_{var}_1950_CurrentYear_GLOBE.nc"

_TLS = threading.local()


def get_session() -> requests.Session:
    session = getattr(_TLS, "session", None)
    if session is None:
        session = requests.Session()
        adapter = HTTPAdapter(pool_connections=2, pool_maxsize=2, max_retries=0)
        session.mount("http://", adapter)
        session.mount("https://", adapter)
        session.headers.update({"Accept-Encoding": "gzip", "User-Agent": "fcp-terraclimate-partition-v2"})
        _TLS.session = session
    return session


def request_one(var: str, lat: float, lon: float, start_year: int, end_year: int, base: str, timeout: int, retries: int) -> pd.DataFrame:
    params = {
        "var": var,
        "latitude": f"{lat:.8f}",
        "longitude": f"{lon:.8f}",
        "time_start": f"{start_year}-01-01T00:00:00Z",
        "time_end": f"{end_year}-12-31T23:59:59Z",
        "accept": "csv",
    }
    last = None
    for attempt in range(retries):
        try:
            r = get_session().get(base.format(var=var), params=params, timeout=timeout)
            r.raise_for_status()
            d = pd.read_csv(io.StringIO(r.text))
            if d.empty:
                raise RuntimeError("empty NCSS response")
            time_col = next(
                (c for c in d.columns if c.lower() in {"date", "time"} or c.lower().startswith("date")),
                None,
            )
            value_col = next(
                (c for c in d.columns if c == var or c.lower().startswith(var.lower() + "[")),
                None,
            )
            if time_col is None or value_col is None:
                raise RuntimeError(f"unexpected columns: {list(d.columns)}")
            z = pd.DataFrame(
                {
                    "time": pd.to_datetime(d[time_col], utc=True, errors="coerce"),
                    "value": pd.to_numeric(d[value_col], errors="coerce"),
                }
            ).dropna(subset=["time"])
            z = (
                z.loc[z.time.dt.year.between(start_year, end_year)]
                .drop_duplicates("time")
                .sort_values("time")
                .reset_index(drop=True)
            )
            if z.empty:
                raise RuntimeError("no requested years in NCSS response")
            return z
        except Exception as exc:  # transport failures are retried; scientific request is unchanged
            last = exc
            if attempt + 1 < retries:
                time.sleep(min(2**attempt, 30))
    raise RuntimeError(f"{var}@({lat},{lon}) failed after {retries} attempts: {last}")


def canonical_jobs(points: pd.DataFrame, variables: list[str]) -> list[dict]:
    unique = (
        points[["tc_lat_idx", "tc_lon_idx", "terraclimate_lat", "terraclimate_lon"]]
        .drop_duplicates()
        .sort_values(["tc_lat_idx", "tc_lon_idx"], kind="stable")
        .reset_index(drop=True)
    )
    jobs: list[dict] = []
    rank = 0
    for r in unique.itertuples(index=False):
        for var in variables:
            jobs.append(
                {
                    "request_rank": rank,
                    "var": var,
                    "lat": float(r.terraclimate_lat),
                    "lon": float(r.terraclimate_lon),
                    "tc_lat_idx": int(r.tc_lat_idx),
                    "tc_lon_idx": int(r.tc_lon_idx),
                }
            )
            rank += 1
    return jobs


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--points", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--qc-out", required=True)
    p.add_argument("--partition-count", type=int, required=True)
    p.add_argument("--partition-index", type=int, required=True)
    p.add_argument("--start-year", type=int, default=1958)
    p.add_argument("--end-year", type=int, default=2025)
    p.add_argument("--variables", nargs="+", default=DEFAULT_VARS)
    p.add_argument("--base-template", default=DEFAULT_BASE)
    p.add_argument("--workers", type=int, default=6)
    p.add_argument("--timeout", type=int, default=120)
    p.add_argument("--retries", type=int, default=8)
    a = p.parse_args()

    if a.partition_count < 1:
        raise SystemExit("partition-count must be >=1")
    if not 0 <= a.partition_index < a.partition_count:
        raise SystemExit("partition-index outside [0, partition-count)")

    pts = pd.read_csv(a.points)
    required = {
        "canonical_name",
        "family",
        "selection_rank",
        "terraclimate_lat",
        "terraclimate_lon",
        "tc_lat_idx",
        "tc_lon_idx",
    }
    if not required <= set(pts):
        raise SystemExit(f"missing point columns: {sorted(required - set(pts))}")

    jobs = canonical_jobs(pts, list(a.variables))
    part = [j for j in jobs if j["request_rank"] % a.partition_count == a.partition_index]
    if not part:
        raise SystemExit("empty compute partition")

    rows: list[pd.DataFrame] = []
    failures: list[dict] = []

    def run(job: dict) -> pd.DataFrame:
        z = request_one(
            job["var"],
            job["lat"],
            job["lon"],
            a.start_year,
            a.end_year,
            a.base_template,
            a.timeout,
            a.retries,
        )
        z.insert(0, "request_rank", int(job["request_rank"]))
        z.insert(1, "var", job["var"])
        z.insert(2, "tc_lat_idx", int(job["tc_lat_idx"]))
        z.insert(3, "tc_lon_idx", int(job["tc_lon_idx"]))
        z.insert(4, "terraclimate_lat", float(job["lat"]))
        z.insert(5, "terraclimate_lon", float(job["lon"]))
        return z

    with ThreadPoolExecutor(max_workers=max(1, a.workers)) as ex:
        futs = {ex.submit(run, job): job for job in part}
        for i, fut in enumerate(as_completed(futs), 1):
            job = futs[fut]
            try:
                rows.append(fut.result())
            except Exception as exc:
                failures.append(
                    {
                        "request_rank": int(job["request_rank"]),
                        "var": job["var"],
                        "tc_lat_idx": int(job["tc_lat_idx"]),
                        "tc_lon_idx": int(job["tc_lon_idx"]),
                        "error": str(exc)[:1000],
                    }
                )
            if i % 20 == 0 or i == len(futs):
                print(
                    {
                        "partition": a.partition_index,
                        "requests_done": i,
                        "requests_total": len(futs),
                        "failures": len(failures),
                    },
                    flush=True,
                )

    qc_path = Path(a.qc_out)
    qc_path.parent.mkdir(parents=True, exist_ok=True)
    if failures:
        qc = {
            "status": "failed",
            "protocol": "terraclimate-fixed-points-partition-v2",
            "partition_count": a.partition_count,
            "partition_index": a.partition_index,
            "global_request_count": len(jobs),
            "assigned_request_count": len(part),
            "assigned_request_rank_min": min(j["request_rank"] for j in part),
            "assigned_request_rank_max": max(j["request_rank"] for j in part),
            "failures": failures,
            "source": a.base_template,
        }
        qc_path.write_text(json.dumps(qc, indent=2) + "\n", encoding="utf-8")
        raise SystemExit(f"partition {a.partition_index}: {len(failures)} requests failed")

    out = pd.concat(rows, ignore_index=True)
    out["time"] = pd.to_datetime(out["time"], utc=True)
    out = out.sort_values(["request_rank", "time"], kind="stable").reset_index(drop=True)

    ranks = sorted(out.request_rank.unique().tolist())
    expected_ranks = sorted(j["request_rank"] for j in part)
    if ranks != expected_ranks:
        raise SystemExit("partition request-rank coverage mismatch")
    if out.duplicated(["request_rank", "time"]).any():
        raise SystemExit("duplicate request-rank/time rows")

    per_request = out.groupby("request_rank").time.nunique()
    op = Path(a.out)
    op.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(op, index=False, compression="gzip")

    qc = {
        "status": "complete",
        "protocol": "terraclimate-fixed-points-partition-v2",
        "partition_count": a.partition_count,
        "partition_index": a.partition_index,
        "global_request_count": len(jobs),
        "assigned_request_count": len(part),
        "assigned_request_rank_min": min(expected_ranks),
        "assigned_request_rank_max": max(expected_ranks),
        "rows": len(out),
        "min_months_per_request": int(per_request.min()),
        "max_months_per_request": int(per_request.max()),
        "start_year": a.start_year,
        "end_year": a.end_year,
        "variables": list(a.variables),
        "source": a.base_template,
        "scientific_guard": "compute partition only; coordinates, variables, source, time interval and values are unchanged from frozen dynamic-climate v1",
    }
    qc_path.write_text(json.dumps(qc, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(qc, indent=2))


if __name__ == "__main__":
    main()
