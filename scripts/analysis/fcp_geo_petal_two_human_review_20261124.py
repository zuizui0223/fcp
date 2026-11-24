#!/usr/bin/env python3
"""Unblind new original-photo geography source audit ONLY after 2 human reviews.

Inputs are separate from frozen source keys, and all 320 source items are
blank when generated. A complete human open-petal estimate is impossible until
all 192 pair and 128 single-photo reviews are independently scored and any
discordance is adjudicated by a third reader. Unit tests use synthetic sheets
and MUST NOT be reported as biological results.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

COLORS={"white","yellow_orange","red_pink","blue_purple"}
BOOL={"yes","no","uncertain"}
CONF={"low","medium","high","uncertain"}

def require(test,msg):
    if not test:
        raise AssertionError("FCP_GEO_PETAL_TWO_HUMAN_SCORE: "+str(msg))

def read(path):
    with path.open(newline="",encoding="utf8") as f:
        return list(csv.DictReader(f))

def keyed(path,n):
    rows=read(path)
    ids=[r["blind_id"] for r in rows]
    require(len(rows)==n and len(set(ids))==n,"duplicate/missing blind IDs: "+str(path))
    return {r["blind_id"]:r for r in rows}

def present(r,keys):
    return any((r.get(k) or "").strip() for k in keys)

def val(r,name,allowed,mandatory=True):
    v=(r.get(name) or "").strip().lower()
    if mandatory:
        require(v in allowed,"incomplete/invalid "+name+" field")
    elif v:
        require(v in allowed,"unexpected "+name+" token")
    return v

def evaluate(r,kind):
    if kind=="pairs":
        fields=["photo_A_open_petal_colour","photo_B_open_petal_colour",
                "both_photos_comparable_open_corollas","same_developmental_stage",
                "non_petal_dominance","photo_A_visible_open_corolla",
                "photo_B_visible_open_corolla","photo_quality_sufficient"]
        if not present(r,fields):
            return None
        flags={k:val(r,k,BOOL) for k in fields if k not in
               ("photo_A_open_petal_colour","photo_B_open_petal_colour")}
        ca=val(r,"photo_A_open_petal_colour",COLORS|{"unknown"},False)
        cb=val(r,"photo_B_open_petal_colour",COLORS|{"unknown"},False)
        confidence=val(r,"reviewer_confidence",CONF,False)
        comparable=(all(flags[k]=="yes" for k in flags if k!="non_petal_dominance")
                    and flags["non_petal_dominance"]=="no"
                    and ca in COLORS and cb in COLORS)
        if comparable:
            require(confidence in CONF,"comparable pair missing reviewer confidence")
            return {"readable":True,"ca":ca,"cb":cb,"mismatch":int(ca!=cb)}
        return {"readable":False,"ca":None,"cb":None,"mismatch":None}
    fields=["open_corolla_visible","predominant_open_petal_colour",
            "other_organ_dominant","background_or_nonflower_contribution",
            "flower_bud_or_fruit","photo_quality_sufficient"]
    if not present(r,fields):
        return None
    flags={k:val(r,k,BOOL) for k in fields if k!="predominant_open_petal_colour"}
    c=val(r,"predominant_open_petal_colour",COLORS|{"unknown"},False)
    confidence=val(r,"reviewer_confidence",CONF,False)
    good=(flags["open_corolla_visible"]=="yes" and
          flags["photo_quality_sufficient"]=="yes" and
          all(flags[k]=="no" for k in (
              "other_organ_dominant","background_or_nonflower_contribution",
              "flower_bud_or_fruit")) and c in COLORS)
    if good:
        require(confidence in CONF,"readable photo missing reviewer confidence")
        return {"readable":True,"colour":c}
    return {"readable":False,"colour":None}

def validate_metadata(a,b,key,kind):
    ids=("photo_A_page","photo_B_page") if kind=="pairs" else ("photo_page",)
    for k in ids:
        require(a[k]==b[k],"reviewers saw different photo/order: "+k)
    if kind=="pairs":
        require(a["photo_A_page"].endswith("/"+key["photo_A_id"]) and
                a["photo_B_page"].endswith("/"+key["photo_B_id"]),
                "pair photo/key URL mismatch")
    else:
        require(a["photo_page"].endswith("/"+key["photo_id"]),
                "single photo/key URL mismatch")

def review(kind,A,B,K,third=None):
    require(A.keys()==B.keys()==K.keys(),"blind IDs cannot be joined")
    if third is not None:
        require(third.keys()==K.keys(),"third-reader blind IDs changed")
    pending=[];status=Counter();result=[]
    for uid in sorted(K):
        a,b=A[uid],B[uid]
        validate_metadata(a,b,K[uid],kind)
        x,y=evaluate(a,kind),evaluate(b,kind)
        if x is None or y is None:
            state="pending"
            decided=None
        elif x==y:
            state="two_reader_consensus"
            decided=x
        else:
            z=evaluate(third[uid],kind) if third else None
            if z is None:
                state="third_reader_required"
                decided=None
            else:
                state="third_reader_adjudicated"
                decided=z
        status[state]+=1
        if decided is None:
            pending.append({"blind_id":uid,"kind":kind,"state":state,
                            "photo_A":a.get("photo_A_page",a.get("photo_page","")),
                            "photo_B":a.get("photo_B_page","")})
        result.append({"source":K[uid],"status":state,
                       "score":decided,"blind_id":uid})
    return result,pending,dict(status)

def petal_geo_desc(pair_result):
    strata=defaultdict(dict)
    for rec in pair_result:
        k=rec["source"]
        group=(k["source_stratum"],k["taxon_id"])
        band=k["distance_group"]
        require(band in ("near","far"),"unrecognized original geographic band")
        require(band not in strata[group],"duplicate same species-stratum near/far row")
        strata[group][band]=rec
    require(len(strata)==96 and all(set(v)=={"near","far"} for v in strata.values()),
            "expected 96 dual-distance source taxa draws")
    out={}
    for key in sorted({group[0] for group in strata}):
        draws=[v for (k,_),v in strata.items() if k==key]
        require(len(draws)==24,"not 24 species-stratum draws")
        complete=[]
        lower=[];upper=[]
        near_valid=far_valid=0
        for v in draws:
            a,b=v["near"]["score"],v["far"]["score"]
            near_valid+=int(a["readable"])
            far_valid+=int(b["readable"])
            amin=amax=a["mismatch"] if a["readable"] else None
            bmin=bmax=b["mismatch"] if b["readable"] else None
            lower.append((bmin if bmin is not None else 0) -
                         (amax if amax is not None else 1))
            upper.append((bmax if bmax is not None else 1) -
                         (amin if amin is not None else 0))
            if a["readable"] and b["readable"]:
                complete.append(b["mismatch"]-a["mismatch"])
        out[key]={
            "n_species_source_draws":24,"near_readable":near_valid,
            "far_readable":far_valid,"n_both_readable":len(complete),
            "far_minus_near_mismatch_complete_case_mean":
                sum(complete)/len(complete) if complete else None,
            "complete_case_is_selective":True,
            "conditional_formal_pair_mismatch_lower":sum(lower)/24,
            "conditional_formal_pair_mismatch_upper":sum(upper)/24,
        }
    return out

def single_photo_desc(single_result):
    out={}
    for label in ("True","False"):
        v=[r for r in single_result
           if str(r["source"]["original_classifiable"])==label]
        require(len(v)==64, "expected 64 single controls per source status")
        readable=[r for r in v if r["score"]["readable"]]
        agreement=(sum(r["score"]["colour"]==
                       r["source"]["original_coarse_photo_label"] for r in readable)
                   if label=="True" else None)
        out[("classified" if label=="True" else "unclassified")]={
            "n_source_photo_draws":len(v),
            "n_comparable_open_corolla":len(readable),
            "n_agree_original_label_among_readable":agreement,
            "readable_source_status_not_population_error_rate":True
        }
    return out

def main(args):
    keys_dir=args.keys_dir
    forms_dir=args.forms_dir
    pairs_key=keyed(keys_dir/"UNBLIND_AFTER_HUMAN_REVIEWS_pairs_key.csv",192)
    singles_key=keyed(keys_dir/"UNBLIND_AFTER_HUMAN_REVIEWS_singles_key.csv",128)
    review_states={}
    all_pending=[]
    pair_result=None
    single_result=None
    for kind,keys in (("pairs",pairs_key),("singles",singles_key)):
        a=keyed(forms_dir/f"reviewer_A_geographic_{kind}_blank.csv",len(keys))
        b=keyed(forms_dir/f"reviewer_B_geographic_{kind}_blank.csv",len(keys))
        third=(keyed(args.third_dir/f"third_reader_geographic_{kind}.csv",len(keys))
               if args.third_dir and
                  (args.third_dir/f"third_reader_geographic_{kind}.csv").exists()
               else None)
        rows,pending,status=review(kind,a,b,keys,third)
        all_pending.extend(pending)
        review_states[kind]={"n_source_items":len(keys),"status_counts":status,
                             "n_unresolved":len(pending)}
        if kind=="pairs":pair_result=rows
        else:single_result=rows
    complete=not all_pending
    data={
        "schema":"fcp_source_blinded_geographic_pair_and_photo_two_human_adjudication_v1",
        "status":"HUMAN_PETAL_REVIEW_COMPLETE" if complete else "NO_HUMAN_PETAL_RESULT",
        "review_status":review_states,
        "n_missing_or_discordant_items":len(all_pending),
        "n_total_items":320,
        "source_stratum_labels_are_NOT_biological_petals":True,
        "scientific_limits":[
            "Blinded review sheets are not results until independently scored",
            "The selected 96 dual-geographic-opportunity taxon-stratum draws are not the full 3600-species P2 estimand",
            "Source-colour stratification and unreadable images make complete-case results selection-dependent",
            "Double-human review verifies appearance on open homologous petals, not genotype, selection or fitness"
        ],
        "open_petal_pair_geography":petal_geo_desc(pair_result) if complete else None,
        "single_image_classification":single_photo_desc(single_result) if complete else None
    }
    args.out_dir.mkdir(parents=True,exist_ok=True)
    (args.out_dir/"human_geographic_petal_review_status.json").write_text(
        json.dumps(data,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    with (args.out_dir/"REQUIRES_INDEPENDENT_REVIEW_OR_ADJUDICATION.csv").open(
            "w",newline="",encoding="utf8") as f:
        fieldnames=["blind_id","kind","state","photo_A","photo_B"]
        writer=csv.DictWriter(f,fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_pending)
    print(json.dumps({k:data[k] for k in
          ("status","review_status","n_missing_or_discordant_items",
           "open_petal_pair_geography","single_image_classification")},indent=2))
    if args.require_complete and not complete:
        raise SystemExit("HUMAN_PETAL_VALIDATION_MISSING: signed reviewer data required")

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--forms-dir",type=Path,required=True)
    p.add_argument("--keys-dir",type=Path,required=True)
    p.add_argument("--third-dir",type=Path)
    p.add_argument("--out-dir",type=Path,required=True)
    p.add_argument("--require-complete",action="store_true")
    main(p.parse_args())
