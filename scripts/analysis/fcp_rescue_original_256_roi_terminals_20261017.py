#!/usr/bin/env python3
"""Fail-closed technical salvage of *original* pixel terminals, never scientific inference.

Both real 256-worker ROI production Actions may lose their assembled photo CSV
if the dependent long-running D/50-km/W statistical job fails. If and ONLY if
all 256 historically identical ROI workers succeeded but the final inference
job failed/cancelled, authorize reassembling those original sealed packets
without downloading images again or changing the colour/biological null.
"""
from __future__ import annotations
import argparse
import json
import re
from pathlib import Path

PINNED={
    38997169400:{
        "label":"new_U100",
        "head_sha":"b79d57686394c1bcbb32fcc2142f0dd19db4ced9",
        "event":"push",
        "head_branch":"main",
        "name":"FCP full 4730 original flower ROI and H1 D 50km W replay",
        "source_photos":248100,
        "preflight_name":"original-import-preflight-no-image-pixels",
        "firewall_artifact":"fcp-all4730-unchanged-original-blind-firewall",
        "partition_artifact_pattern":"fcp-all4730-original-b*",
    },
    38992543115:{
        "label":"wave4",
        "head_sha":"ef448c9514e0c702e7887e00b288ee966728b0d5",
        "event":"pull_request",
        "head_branch":"analysis/fcp-wave4-roi-import-repair-20261017",
        "name":"FCP wave4 identical original photo ROI and D W spatial replay",
        "source_photos":49800,
        "preflight_name":None,
        "firewall_artifact":"fcp-wave4-unchanged-original-blind-firewall",
        "partition_artifact_pattern":"fcp-wave4-original-b*",
    },
}
PIXEL_JOB=re.compile(r"unmodified-original-roi-workers \(([01]), (\d{1,2}), ([0-3])\)")
FINAL="same-original-H1-50km-W"


def decide(run:dict,job_pages:list[dict])->dict:
    rid=run.get("id")
    if rid not in PINNED:
        raise RuntimeError("not one of the two exact frozen original ROI runs")
    x=PINNED[rid]
    for key in ("head_sha","event","head_branch","name"):
        if run.get(key)!=x[key]:
            raise RuntimeError(f"{x['label']}: immutable original {key} drift")
    if (run.get("run_attempt")!=1
        or run.get("head_repository",{}).get("full_name")!="zuizui0223/fcp"):
        raise RuntimeError("original-source first-run or repository integrity drift")
    jobs=[j for page in job_pages for j in page.get("jobs",[])]
    by_name={}
    for j in jobs:
        name=j.get("name","")
        if name in by_name:
            raise RuntimeError("duplicate original photo ROI job records")
        by_name[name]=j
    if len({j.get("id") for j in jobs})!=len(jobs):
        raise RuntimeError("duplicate original source worker job IDs")
    actual=set()
    for j in jobs:
        m=PIXEL_JOB.fullmatch(j.get("name",""))
        if not m:
            continue
        b,s,p=(int(v) for v in m.groups())
        if not(0<=b<2 and 0<=s<32 and 0<=p<4):
            raise RuntimeError("original pixel worker matrix index drift")
        if (b,s,p) in actual:
            raise RuntimeError("duplicate original pixel ROI matrix assignment")
        actual.add((b,s,p))
    expected={(b,s,p) for b in range(2) for s in range(32) for p in range(4)}
    if actual!=expected:
        raise RuntimeError("original 256-ROI worker matrix identities changed")
    status={j["name"]:j.get("conclusion") for j in jobs}
    workers_all_success=all(
        status[f"unmodified-original-roi-workers ({b}, {s}, {p})"]=="success"
        for b,s,p in expected
    )
    frozen_success=(status.get("freeze-new-original-pixel-IDs")=="success")
    preflight_ok=(x["preflight_name"] is None or status.get(x["preflight_name"])=="success")
    final_status=status.get(FINAL)
    complete=(run.get("status")=="completed")
    if (complete and run.get("conclusion")=="success"
        and workers_all_success and frozen_success and preflight_ok
        and final_status=="success"):
        reason="FULL_ACTUAL_ORIGINAL_PIXEL_AND_SCIENCE_ALREADY_COMPLETE"
        ready=False
    elif (complete and run.get("conclusion") in ("failure","cancelled","timed_out")
          and workers_all_success and frozen_success and preflight_ok
          and final_status in ("failure","cancelled","timed_out")):
        reason="ALL_256_ORIGINAL_PIXEL_TERMINALS_SUCCESSFUL_ONLY_STATISTICS_FAILED"
        ready=True
    else:
        reason="NO_COMPLETE_AUTHENTIC_256_ROI_PIXEL_TERMINAL_RESCUE_AUTHORIZED"
        ready=False
    return {
      "schema":"fcp_restore_authentic_original_roi_256_terminal_if_statistics_failed_v1",
      "ready":ready,"reason":reason,"source_run":rid,"cohort":x["label"],
      "source_original_photo_IDs":x["source_photos"],
      "terminal_worker_jobs":len(actual),
      "original_worker_jobs_successful":sum(
          status.get(f"unmodified-original-roi-workers ({b}, {s}, {p})")=="success"
          for b,s,p in expected
      ),
      "original_final_science_job_conclusion":final_status,
      "original_freeze_job_successful":frozen_success,
      "original_preflight_job_successful":preflight_ok,
      "original_pixel_measurement_method_replaced":False,
      "original_H1_D_50km_H2_W_not_rerun":True,
      "forbidden_early_success_claim":not ready,
      "firewall_artifact":x["firewall_artifact"],
      "partition_pattern":x["partition_artifact_pattern"],
    }


def main()->None:
    ap=argparse.ArgumentParser()
    ap.add_argument("--source",type=Path,required=True)
    ap.add_argument("--page",type=Path,action="append",required=True)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--github-output",type=Path,required=True)
    a=ap.parse_args()
    rec=decide(json.loads(a.source.read_text()),[json.loads(f.read_text()) for f in a.page])
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(rec,indent=2,sort_keys=True)+"\n")
    with a.github_output.open("a") as f:
        for k,v in (
            ("ready",str(rec["ready"]).lower()),
            ("label",rec["cohort"]),
            ("source_run",str(rec["source_run"])),
            ("firewall_artifact",rec["firewall_artifact"]),
            ("partition_pattern",rec["partition_pattern"]),
        ):
            f.write(k+"="+v+"\n")
    print("FCP_AUTHENTIC_ORIGINAL_256_TERMINAL_RESCUE_GATE",json.dumps(rec,sort_keys=True))


if __name__=="__main__":
    main()
