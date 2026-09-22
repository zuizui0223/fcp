#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

import numpy as np
import pandas as pd

EXPECTED_ROWS = 49_900
EXPECTED_PARTITIONS = 256

def main() -> None:
    p=argparse.ArgumentParser()
    p.add_argument("--results-dir", type=Path, required=True)
    p.add_argument("--firewall-dir", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    args=p.parse_args()

    tech_files=sorted(args.results_dir.rglob("technical_b*_s*_p*.csv"))
    acq_files=sorted(args.results_dir.rglob("acquisition_b*_s*_p*.csv"))
    if len(tech_files) != EXPECTED_PARTITIONS or len(acq_files) != EXPECTED_PARTITIONS:
        raise RuntimeError(
            f"partition census incomplete: technical={len(tech_files)} acquisition={len(acq_files)}"
        )

    acq=pd.concat([pd.read_csv(x,dtype=str).fillna("") for x in acq_files], ignore_index=True)
    if len(acq) != EXPECTED_ROWS or acq["measurement_id"].nunique() != EXPECTED_ROWS:
        raise RuntimeError("acquisition receipt census is incomplete or duplicated")

    tech=pd.concat([pd.read_csv(x,dtype={"measurement_id":str}) for x in tech_files], ignore_index=True)
    if tech["measurement_id"].nunique() != len(tech):
        raise RuntimeError("technical rows contain duplicate IDs")
    if not set(tech["measurement_id"]).issubset(set(acq["measurement_id"])):
        raise RuntimeError("technical rows are outside frozen acquisition census")

    merged=acq.merge(
        tech,
        on=["measurement_id","reacquired_image_sha256"],
        how="left",
        validate="one_to_one",
        suffixes=("","_technical"),
    )
    if len(merged) != EXPECTED_ROWS:
        raise RuntimeError("technical merge lost frozen rows")
    failed=merged["acquisition_status"].ne("acquired_and_decode_verified")
    merged.loc[failed,"technical_status"]="image_acquisition_failed"
    merged.loc[failed,"failure_reason"]=merged.loc[failed,"failure_reason"].where(
        merged.loc[failed,"failure_reason"].notna() & merged.loc[failed,"failure_reason"].ne(""),
        merged.loc[failed,"failure_reason_technical"] if "failure_reason_technical" in merged else ""
    )

    near=pd.to_numeric(merged["near_clip_fraction"], errors="coerce")
    available=near.notna() & np.isfinite(near)
    if int(available.sum()) == 0:
        raise RuntimeError("no technical highlight diagnostics available")
    q95=float(np.quantile(near[available].to_numpy(float),0.95))
    threshold=float(max(0.01,q95))
    high=available & (near > threshold)

    out=args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out/"technical_table.csv.gz",index=False,compression="gzip",lineterminator="\n")
    pd.DataFrame({"measurement_id":merged.loc[high,"measurement_id"].astype(str)}).to_csv(
        out/"high_clip_ids.csv",index=False,lineterminator="\n"
    )
    shutil.copy2(args.firewall_dir/"sealed_join_key.csv", out/"sealed_join_key.csv")
    summary={
        "schema":"third_cohort_highlight_technical_seal_v1",
        "status":"response_blind_high_clip_set_frozen",
        "rows":EXPECTED_ROWS,
        "partition_receipts":EXPECTED_PARTITIONS,
        "acquired_rows":int((merged["acquisition_status"]=="acquired_and_decode_verified").sum()),
        "acquisition_failed_rows":int(failed.sum()),
        "highlight_metrics_available_rows":int(available.sum()),
        "highlight_metrics_unavailable_rows":int((~available).sum()),
        "near_clip_q95":q95,
        "high_clip_threshold":threshold,
        "high_clip_rule":"near_clip_fraction > max(0.01, q95)",
        "high_clip_rows":int(high.sum()),
        "biological_outcomes_opened":False,
        "species_opened":False,
        "morph_opened":False,
    }
    (out/"technical_summary.json").write_text(
        json.dumps(summary,indent=2,sort_keys=True)+"\n",encoding="utf-8"
    )
    print(json.dumps(summary,indent=2,sort_keys=True))

if __name__=="__main__":
    main()
