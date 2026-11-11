#!/usr/bin/env python3
"""Produce a fixed, colour-label-free metadata-check queue for FCP site-year gaps.

The source is the measured-photo *opportunity* CSV from an already completed
historical audit. Queue all gap-1 and gap-2 original taxa, not merely the
minimum number that would reach a sample-size threshold, so future checks
cannot choose successes after seeing new biological outcomes.
NO IMAGE API REQUEST OR NEW PHOTO COLOUR IS MADE BY THIS SCRIPT.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import pandas as pd

COHORTS = ("discovery","validation","third")
REQUIRED_ORIGINAL = {"discovery":369,"validation":363,"third":377}
EXPECTED_FROZEN_GAP = {"discovery":{"0":10,"1":37,"2":49},
                       "validation":{"0":10,"1":22,"2":60},
                       "third":{"0":20,"1":33,"2":62}}
GOAL = 30


def build_queue(table: pd.DataFrame, receipt: dict) -> tuple[dict,pd.DataFrame]:
    if receipt.get("schema")!="fcp_same_siteyear_same_month_replenishment_gap_v1":
        raise ValueError("Wrong frozen source receipt")
    if receipt.get("original_all_cohorts_gate_hold") is not True:
        raise ValueError("Historical feasibility HOLD unexpectedly changed")
    required={"cohort","inat_taxon_id","species","two_year_month_already_observed",
              "minimum_additional_observer_photo_slots","candidate_years",
              "candidate_month","anchor_photo_id","existing_observer_counts",
              "original_anchor_eligible"}
    if not required.issubset(table):
        raise ValueError("Missing opportunity columns: "+str(sorted(required-set(table))))
    if table.duplicated(["cohort","inat_taxon_id"]).any():
        raise ValueError("Duplicate historical cohort-taxon identity")
    if len(table)!=sum(REQUIRED_ORIGINAL.values()):
        raise ValueError("Incomplete original 1109 taxon opportunity sample")
    if set(table.cohort)!=set(COHORTS):
        raise ValueError("Historical cohort mismatch")
    if {"morph","white","flower_colour","color","colour"} & set(table):
        raise ValueError("Photo-colour labels must not be available to queue construction")
    requested=[]
    cohorts={}
    for cohort in COHORTS:
        group=table.loc[table.cohort==cohort].copy().sort_values("inat_taxon_id",kind="stable")
        if len(group)!=REQUIRED_ORIGINAL[cohort]:
            raise ValueError(f"{cohort}: frozen sample denominator changed")
        recorded=receipt["cohorts"][cohort]["minimal_additional_observer_photo_slots_histogram"]
        if {str(j):int(recorded.get(str(j),0)) for j in range(3)}!=EXPECTED_FROZEN_GAP[cohort]:
            raise ValueError(f"{cohort}: source CI gap counts changed")
        actual={str(i):0 for i in range(3)}
        for _,s in group.iterrows():
            gap=s.minimum_additional_observer_photo_slots
            if pd.isna(gap):
                continue
            if float(gap) not in (0.,1.,2.):
                raise ValueError("Unexpected missing-observer slot count")
            gap=int(gap)
            actual[str(gap)]+=1
            if gap==0:
                if not bool(s.original_anchor_eligible):
                    raise ValueError("Gap 0 is not qualified")
                continue
            if bool(s.original_anchor_eligible):
                raise ValueError("Nonzero gap erroneously marked qualified")
            years=[int(x) for x in str(s.candidate_years).split(",")]
            obs=[int(x) for x in str(s.existing_observer_counts).split(",")]
            if len(years)!=2 or len(set(years))!=2 or len(obs)!=2 or years!=sorted(years):
                raise ValueError("Malformed year/observer matching")
            if not 1<=int(s.candidate_month)<=12:
                raise ValueError("Bad calendar month")
            needs=[max(0,2-n) for n in obs]
            if sum(needs)!=gap:
                raise ValueError("Opportunity gap and year observer counts disagree")
            for y,need,have in zip(years,needs,obs):
                if need==0:
                    continue
                requested.append({
                    "priority_tier":gap,
                    "cohort":cohort,
                    "inat_taxon_id":str(s.inat_taxon_id),
                    "species":str(s.species),
                    "source_anchor_photo_id":str(s.anchor_photo_id),
                    "calendar_month":int(s.candidate_month),
                    "target_year":int(y),
                    "existing_distinct_observers_in_target_year":int(have),
                    "required_additional_distinct_observer_photos":int(need),
                    "new_photo_or_observation_id": "NOT_RETRIEVED",
                    "live_metadata_status":"NOT_CHECKED",
                })
        if actual!=EXPECTED_FROZEN_GAP[cohort]:
            raise ValueError(f"{cohort}: source row counts disagree with frozen receipt")
        have=actual["0"]
        need=max(0,GOAL-have)
        gapone=actual["1"]
        cohorts[cohort]={
            "original_anchored_source_species":len(group),
            "already_complete":have,
            "one_slot_species_available_to_check":gapone,
            "two_slot_species_available_to_check":actual["2"],
            "minimum_hypothetical_new_observer_photo_slots_to_goal":need,
            "goal_theoretically_reachable_from_gap1_only":gapone>=need,
            "required_fraction_of_gap1_candidates_successful_if_gap2_unused":need/gapone if gapone else None,
        }
    frame=pd.DataFrame(requested).sort_values(
        ["priority_tier","cohort","inat_taxon_id","target_year"],kind="stable"
    ).reset_index(drop=True)
    summary={
        "schema":"fcp_replenishment_metadata_check_queue_v1",
        "status":"HISTORIC_PHOTO_METADATA_CHECK_QUEUE_ONLY_NO_API_REQUESTS",
        "original_hold_unchanged":True,
        "confirmatory_decisions_changed":False,
        "cohorts":cohorts,
        "n_gap1_species":sum(v["one_slot_species_available_to_check"] for v in cohorts.values()),
        "n_gap2_species":sum(v["two_slot_species_available_to_check"] for v in cohorts.values()),
        "total_required_observer_photo_slots_queued":int(frame.required_additional_distinct_observer_photos.sum()),
        "minimum_hypothetical_new_photo_slots_for_all_three_cohorts":sum(v["minimum_hypothetical_new_observer_photo_slots_to_goal"] for v in cohorts.values()),
        "request_selection":"all historical gap1 and gap2 taxa, sorted without outcomes, source-photo anchor/matched source-year/month only",
        "nonclaims":[
            "Queue rows describe necessary photo and independent observer metadata, not images found",
            "Original gap0 10km source eligibility HOLD is unchanged",
            "No new photo ID, live API response, pixel, pigment, climate anomaly or fitness was opened",
            "The sample is from historically outcome-exposed species and is not the frozen 2000+730 independent confirmation",
            "Any eventual photo must be verified for exact species identity, observer separation, date, distance, licence and true flower visibility before analysis",
        ]
    }
    if summary["n_gap1_species"]!=92 or summary["n_gap2_species"]!=171:
        raise ValueError("Historic three-cohort candidate counts changed")
    return summary,frame


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source-csv",required=True,type=Path)
    ap.add_argument("--source-receipt",required=True,type=Path)
    ap.add_argument("--outdir",required=True,type=Path)
    a=ap.parse_args()
    src=pd.read_csv(a.source_csv,dtype={"cohort":"string","inat_taxon_id":"string","anchor_photo_id":"string"},low_memory=False)
    receipt=json.loads(a.source_receipt.read_text())
    report,queue=build_queue(src,receipt)
    a.outdir.mkdir(parents=True,exist_ok=True)
    queue.to_csv(a.outdir/"source_only_observer_photo_query_queue.csv",index=False)
    (a.outdir/"result.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2,sort_keys=True))

if __name__=="__main__":
    main()
