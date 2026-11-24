#!/usr/bin/env python3
"""FCP: partial-identification/MNAR sensitivity for 62 unreadable photo pairs.

Replays, and does not alter, Table S11's 128 selected thermal-tail pair draws.
The visual scores are from ONE AI reader, NOT independent botanical ground truth.

Use from repository root:
  python scripts/analysis/fcp_organ_unreadable_mnar_bounds_20261124.py \
    --out-dir /tmp/fcp-petal-mnar --bootstrap 5000
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import random
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT/"manuscripts/climate_colour/Table_S11_pair_audit.csv"
SEED = 20261122
FULL_BASE = {"white": 2.63, "mismatch": -0.01}
EXPECTED = {
    "white": {"original_draw": 1.4151585511609,
              "retain_original": 0.3697729911886505},
    "mismatch": {"original_draw": -0.013722614152666956,
                 "retain_original": -1.7331135069932402},
}


def require(ok, explanation):
    if not ok:
        raise AssertionError("FCP_UNREADABLE_SENSITIVITY: "+str(explanation))


def load():
    with SOURCE.open(newline="",encoding="utf-8") as f:
        src=list(csv.DictReader(f))
    require(len(src)==128, "not 128 source pair draws")
    ids=[r["id"] for r in src]
    require(ids==[f"P{i:03d}" for i in range(1,129)], "original blind IDs changed")
    stratum=Counter((int(r["tail"]),int(r["original_pair_mismatch"])) for r in src)
    require(stratum=={(-1,0):32,(-1,1):32,(1,0):32,(1,1):32},
            "not four equal original thermal/label strata")
    require(sum(r["valid"]=="True" for r in src)==66,
            "single-reader visual completeness changed")
    require(sum(r["valid"]=="False" for r in src)==62,
            "unreadable visual count changed")
    for r in src:
        require(r["valid"] in ("True","False"), "readability missing")
        require(float(r["stratum_weight"])>0, "nonpositive source design weight")
        require(r["white"] in ("0.0","0.5","1.0"),
                "source pair white fraction changed")
        if r["valid"]=="True":
            require(r["viswhite"] in ("0.0","0.5","1.0") and
                    r["visdiff"] in ("True","False"),"readable pair lacks visual score")
        else:
            require(r["viswhite"]=="" and r["visdiff"]=="",
                    "unreadable pair incorrectly given a visual score")
    return src


def label(r,kind,visual=False):
    if kind=="white":
        return float(r["viswhite"] if visual else r["white"])
    return float((r["visdiff"]=="True") if visual
                 else r["original_pair_mismatch"]=="1")


def tail_difference(rows,values):
    group={}
    for t in (-1,1):
        rr=[(r,x) for r,x in zip(rows,values) if int(r["tail"])==t]
        denom=sum(float(r["stratum_weight"]) for r,x in rr)
        require(denom>0,"empty thermal tail")
        group[t]=sum(float(r["stratum_weight"])*x for r,x in rr)/denom
    return 100*(group[1]-group[-1])


def scoring(rows,kind):
    # All coefficients below are on the sample-standardized pair contrast
    # scale. The correction is applied to the separately measured rounded
    # full-pair-frame mean, NOT the raw audit draw contrast itself.
    original=[label(r,kind) for r in rows]
    original_delta=tail_difference(rows,original)
    visual=[]
    bystratum={}
    for t in (-1,1):
        for orig in (0,1):
            valid=[r for r in rows if int(r["tail"])==t
                   and int(r["original_pair_mismatch"])==orig
                   and r["valid"]=="True"]
            require(valid,"no readable photo in source audit stratum")
            bystratum[(t,orig)]=sum(label(r,kind,True) for r in valid)/len(valid)
    for mode in ("retain_original","missing_as_half","mar_within_four_strata",
                 "all_missing_lower","all_missing_upper"):
        corrected=[]
        for r in rows:
            t=int(r["tail"])
            orig=int(r["original_pair_mismatch"])
            if r["valid"]=="True":
                x=label(r,kind,True)
            elif mode=="retain_original":
                x=label(r,kind)
            elif mode=="missing_as_half":
                x=.5
            elif mode=="mar_within_four_strata":
                x=bystratum[(t,orig)]
            elif mode=="all_missing_lower":
                x=0. if t==1 else 1.
            else:
                x=1. if t==1 else 0.
            corrected.append(x)
        visual.append((mode,FULL_BASE[kind]+tail_difference(rows,corrected)-original_delta))
    result=dict(visual)
    require(result["all_missing_lower"]<=result["retain_original"]<=result["all_missing_upper"],
            "unreadable bounds do not include retain-original estimate")
    require(result["all_missing_lower"]<=result["mar_within_four_strata"]<=result["all_missing_upper"],
            "unreadable bounds do not include stratified MAR estimate")
    require(result["all_missing_lower"]<=result["missing_as_half"]<=result["all_missing_upper"],
            "unreadable bounds do not include neutral estimate")
    # The MAR interpolation is just a sensitivity parameter, NOT an estimated
    # proportion of errorful images or proof the missingness process is random.
    denominator=result["mar_within_four_strata"]-result["retain_original"]
    tipping=(-result["retain_original"]/denominator) if denominator else None
    return {
        "original_draw_contrast_pp":original_delta,
        "source_full_frame_baseline_pp_rounded":FULL_BASE[kind],
        "model_point_estimates_pp":result,
        "stratum_readable_visual_mean":{
            f"tail_{t}_original_mismatch_{m}":v
            for (t,m),v in bystratum.items()
        },
        "formal_marginal_unreadable_bounds_pp":{
            "lower":result["all_missing_lower"],"upper":result["all_missing_upper"]
        },
        "retain_to_stratified_mar_zero_crossing_fraction":tipping,
        "assumptions":{
            "retain_original":"Unreadable pairs retain original photo labels.",
            "missing_as_half":"All unreadable pair scores are set to 0.5; not observed.",
            "mar_within_four_strata":"An unreadable pair has the same mean hypothetical petal score as readable pairs with identical warm/cold tail and original label agreement. Unverifiable assumption.",
            "all_missing_lower":"For lower warm-minus-cold bound, all unreadable warm scores are 0 and cold scores 1.",
            "all_missing_upper":"For upper warm-minus-cold bound, all unreadable warm scores are 1 and cold scores 0.",
            "formal_bounds_scope":"Mathematical marginal completions of 62 undefined single-reader petal scores conditional on 66 fixed AI scores, source draw weights and rounded full-frame baseline. Not population-wide petal-colour prevalence nor biological identified bounds if no homologous petals exist."
        }
    }


def bootstrap(rows, n):
    # Reproduce the prior retention point and bootstrap source, then evaluate
    # alternative hypothetical completions using exactly the same resampled
    # source audit draws, with four fixed strata and 32 draws each.
    groups={}
    for r in rows:
        groups.setdefault((int(r["tail"]),int(r["original_pair_mismatch"])),[]).append(r)
    ordered=[groups[k] for k in sorted(groups)]
    rng=random.Random(SEED)
    samples={kind:{mode:[] for mode in
                   ("retain_original","missing_as_half","mar_within_four_strata")}
             for kind in FULL_BASE}
    invalid=0
    for _ in range(n):
        draw=[rng.choice(group) for group in ordered for _ in range(32)]
        # Bootstrap draws with no visually readable pair in a source stratum
        # are undefined for the MAR estimator; report and drop such replicates
        # for ALL scenarios to preserve paired comparisons.
        if any(not any(r["valid"]=="True" for r in draw if
                       (int(r["tail"]),int(r["original_pair_mismatch"]))==k)
               for k in groups):
            invalid+=1
            continue
        for kind in FULL_BASE:
            point=scoring(draw,kind)["model_point_estimates_pp"]
            for mode in samples[kind]:
                samples[kind][mode].append(point[mode])
    require(invalid<=max(3,n//100),
            "too many undefined bootstrap source strata")
    out={"seed":SEED,"resamples_requested":n,"undefined_stratum_draws":invalid,
         "resamples_accepted":n-invalid,
         "interval_type":"percentile_95",
         "intervals":{},
         "method":"Independent resampling of 32 audited draws within each of four source strata; single AI fixed, rounded full-frame anchor fixed"}
    for kind in FULL_BASE:
        out["intervals"][kind]={}
        for mode,vals in samples[kind].items():
            vals=sorted(vals)
            def q(p):
                k=(len(vals)-1)*p
                a=int(math.floor(k));b=min(a+1,len(vals)-1)
                return vals[a]+(vals[b]-vals[a])*(k-a)
            out["intervals"][kind][mode]=[q(.025),q(.975)]
    # Check the *retention* CI matches the published 2026-11-22 replay to
    # tolerance of the same random draws and source stratum order.
    lo,hi=out["intervals"]["white"]["retain_original"]
    require(abs(lo-(-9.859244398876653))<1e-6 and
            abs(hi-10.314225197556857)<1e-6 if n==5000 else True,
            "archived 5k paired-boot CI changed")
    return out


def main(out,bs):
    rows=load()
    white=scoring(rows,"white")
    mismatch=scoring(rows,"mismatch")
    require(abs(white["original_draw_contrast_pp"] -
                EXPECTED["white"]["original_draw"])<1e-8,
            "original source warm/cold white contrast drift")
    for kind,point in (("white",white),("mismatch",mismatch)):
        require(abs(point["model_point_estimates_pp"]["retain_original"]-
                    EXPECTED[kind]["retain_original"])<1e-8,
                kind+" retention source audit mismatch")
    result={
        "schema":"fcp_128_pair_unreadable_partial_identification_20261124",
        "status":"exploratory_source_table_sensitivity_NOT_human_petal_ground_truth",
        "source_path":str(SOURCE.relative_to(ROOT)),
        "source_sha256":hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "n_selected_draws":128,
        "n_unique_taxa":len({r["inat_taxon_id"] for r in rows}),
        "n_readable_by_single_AI":66,
        "n_unreadable_by_single_AI":62,
        "weight":"stratum_weight; source-pair draw sample, not species-equal global atlas",
        "stratum_counts":{
            f"tail_{t}_original_mismatch_{m}":{
                "n":sum(int(r["tail"])==t and int(r["original_pair_mismatch"])==m for r in rows),
                "n_readable":sum(int(r["tail"])==t and int(r["original_pair_mismatch"])==m and r["valid"]=="True" for r in rows)
            } for t in (-1,1) for m in (0,1)
        },
        "white_pair_share_pp":white,
        "category_mismatch_pp":mismatch,
        "bootstrap":bootstrap(rows,bs),
        "scientific_limits":[
            "All measurements derive from one earlier AI visual reader, not two independent botanists.",
            "62 pairs cannot be visually compared; some may not have any comparable open floral organs, so completions are formal numerical sensitivity not observed biological outcomes.",
            "Sample is selected by original photo-label agreement and warm/cold extremes; do not treat as global confusion matrix.",
            "The +2.63 pp full-frame white contrast is a rounded external anchor; it is not the 2837-species conditional thermal-white coefficient.",
            "Do not use these bounds to invalidate or calibrate the original 3600-species photographic spatial estimator or a pigment fitness mechanism.",
            "The within-stratum MAR scenario is an unsupported sensitivity assumption, not a missingness correction."
        ]
    }
    out.mkdir(parents=True,exist_ok=True)
    (out/"fcp_unreadable_visual_pair_mnar_sensitivity_20261124.json").write_text(
        json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({k:result[k] for k in
                      ["n_selected_draws","n_readable_by_single_AI",
                       "n_unreadable_by_single_AI","stratum_counts",
                       "white_pair_share_pp","category_mismatch_pp",
                       "bootstrap"]},indent=2,ensure_ascii=False))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out-dir",type=Path,default=ROOT/"output/fcp_mnar")
    p.add_argument("--bootstrap",type=int,default=5000)
    args=p.parse_args()
    main(args.out_dir,args.bootstrap)
