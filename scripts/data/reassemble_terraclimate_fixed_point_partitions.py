#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

DEFAULT_VARS = ["tmin", "tmax", "ppt", "def", "vpd"]


def canonical_jobs(points: pd.DataFrame, variables: list[str]) -> pd.DataFrame:
    unique = (
        points[["tc_lat_idx", "tc_lon_idx", "terraclimate_lat", "terraclimate_lon"]]
        .drop_duplicates()
        .sort_values(["tc_lat_idx", "tc_lon_idx"], kind="stable")
        .reset_index(drop=True)
    )
    rows = []
    rank = 0
    for r in unique.itertuples(index=False):
        for var in variables:
            rows.append(
                {
                    "request_rank": rank,
                    "var": var,
                    "tc_lat_idx": int(r.tc_lat_idx),
                    "tc_lon_idx": int(r.tc_lon_idx),
                    "terraclimate_lat": float(r.terraclimate_lat),
                    "terraclimate_lon": float(r.terraclimate_lon),
                }
            )
            rank += 1
    return pd.DataFrame(rows)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--points", required=True)
    p.add_argument("--shard-dir", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--qc-out", required=True)
    p.add_argument("--partition-count", type=int, required=True)
    p.add_argument("--start-year", type=int, default=1958)
    p.add_argument("--end-year", type=int, default=2025)
    p.add_argument("--variables", nargs="+", default=DEFAULT_VARS)
    a = p.parse_args()

    pts = pd.read_csv(a.points)
    expected = canonical_jobs(pts, list(a.variables))
    shard_dir = Path(a.shard_dir)
    csvs = sorted(shard_dir.glob("partition_*.csv.gz"))
    qcs = sorted(shard_dir.glob("partition_*_qc.json"))

    if len(csvs) != a.partition_count:
        raise SystemExit(f"expected {a.partition_count} partition CSVs, found {len(csvs)}")
    if len(qcs) != a.partition_count:
        raise SystemExit(f"expected {a.partition_count} partition QCs, found {len(qcs)}")

    qc_records = [json.loads(p.read_text(encoding="utf-8")) for p in qcs]
    if {q.get("status") for q in qc_records} != {"complete"}:
        raise SystemExit("one or more partition QCs are not complete")
    indices = sorted(int(q["partition_index"]) for q in qc_records)
    if indices != list(range(a.partition_count)):
        raise SystemExit("partition-index coverage is not exact")
    if {int(q["partition_count"]) for q in qc_records} != {a.partition_count}:
        raise SystemExit("partition-count mismatch among shards")
    if {int(q["global_request_count"]) for q in qc_records} != {len(expected)}:
        raise SystemExit("global request-count mismatch among shards")

    pieces = []
    for path in csvs:
        d = pd.read_csv(
            path,
            usecols=[
                "request_rank",
                "var",
                "tc_lat_idx",
                "tc_lon_idx",
                "terraclimate_lat",
                "terraclimate_lon",
                "time",
                "value",
            ],
        )
        d["time"] = pd.to_datetime(d["time"], utc=True, errors="coerce")
        pieces.append(d)
    long = pd.concat(pieces, ignore_index=True)
    if long.time.isna().any():
        raise SystemExit("unparseable time in partition output")
    if long.duplicated(["request_rank", "time"]).any():
        raise SystemExit("duplicate request-rank/time rows across partitions")

    actual_keys = (
        long[["request_rank", "var", "tc_lat_idx", "tc_lon_idx", "terraclimate_lat", "terraclimate_lon"]]
        .drop_duplicates()
        .sort_values("request_rank", kind="stable")
        .reset_index(drop=True)
    )
    expected_cmp = expected.sort_values("request_rank", kind="stable").reset_index(drop=True)
    if len(actual_keys) != len(expected_cmp):
        raise SystemExit(f"request-key count mismatch: {len(actual_keys)} != {len(expected_cmp)}")
    for col in ["request_rank", "var", "tc_lat_idx", "tc_lon_idx"]:
        if not actual_keys[col].equals(expected_cmp[col]):
            raise SystemExit(f"request-key mismatch in {col}")
    for col in ["terraclimate_lat", "terraclimate_lon"]:
        if not (actual_keys[col].round(10).to_numpy() == expected_cmp[col].round(10).to_numpy()).all():
            raise SystemExit(f"request coordinate mismatch in {col}")

    request_months = long.groupby("request_rank").time.nunique()
    if request_months.index.nunique() != len(expected):
        raise SystemExit("not every frozen request has a time series")

    keycols = ["tc_lat_idx", "tc_lon_idx", "terraclimate_lat", "terraclimate_lon", "time"]
    wide = (
        long.pivot(index=keycols, columns="var", values="value")
        .reset_index()
        .rename_axis(columns=None)
    )
    missing_vars = [v for v in a.variables if v not in wide.columns]
    if missing_vars:
        raise SystemExit(f"missing variables after pivot: {missing_vars}")

    out = pts.merge(
        wide,
        on=["tc_lat_idx", "tc_lon_idx", "terraclimate_lat", "terraclimate_lon"],
        how="left",
        validate="many_to_many",
    )
    out["time"] = pd.to_datetime(out["time"], utc=True, errors="coerce")
    out["year"] = out.time.dt.year
    out["month"] = out.time.dt.month
    per_point = out.groupby(["canonical_name", "selection_rank"]).time.nunique()

    # Preserve the exact v1 workflow coverage boundary: at least 800 distinct
    # monthly timestamps per frozen point, with no value imputation.
    if int(per_point.min()) < 800:
        raise SystemExit(f"v1 monthly coverage gate failed: min={int(per_point.min())}")
    if len(per_point) != len(pts):
        raise SystemExit("not every frozen species-location row is represented")

    op = Path(a.out)
    op.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(op, index=False, compression="gzip")

    qc = {
        "status": "complete",
        "protocol": "terraclimate-fixed-points-reassembly-v2",
        "partition_count": a.partition_count,
        "partitions_present": len(csvs),
        "fixed_points": len(pts),
        "unique_grid_cells": int(expected[["tc_lat_idx", "tc_lon_idx"]].drop_duplicates().shape[0]),
        "requests_expected": len(expected),
        "requests_reassembled": int(actual_keys.request_rank.nunique()),
        "start_year": a.start_year,
        "end_year": a.end_year,
        "variables": list(a.variables),
        "min_months_per_request": int(request_months.min()),
        "max_months_per_request": int(request_months.max()),
        "min_months_per_point": int(per_point.min()),
        "max_months_per_point": int(per_point.max()),
        "rows": len(out),
        "value_imputation": False,
        "scientific_guard": "exact union of all frozen requests; partition topology cannot select or replace values",
    }
    Path(a.qc_out).write_text(json.dumps(qc, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(qc, indent=2))


if __name__ == "__main__":
    main()
