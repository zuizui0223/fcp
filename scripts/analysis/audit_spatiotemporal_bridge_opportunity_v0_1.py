#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DESIGN=ROOT/"data"/"spatiotemporal_bridge_opportunity_design_v0_1.json"


def sha256(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv_species(path:Path)->set[tuple[str,str]]:
    out=set()
    with path.open(newline="",encoding="utf-8") as f:
        r=csv.DictReader(f,delimiter="\t")
        req={"inat_taxon_id","species"}
        if not req.issubset(r.fieldnames or []):
            raise RuntimeError(f"selected manifest missing {sorted(req)}; got {r.fieldnames}")
        for row in r:
            out.add((str(row["inat_taxon_id"]).strip(),str(row["species"]).strip()))
    return out


def select_hash(salt:str,taxon_id:str,species:str)->str:
    return hashlib.sha256(f"{salt}|{taxon_id}|{species}".encode()).hexdigest()


def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--candidate",type=Path,required=True)
    ap.add_argument("--third-selected",type=Path,required=True)
    ap.add_argument("--out-json",type=Path,required=True)
    ap.add_argument("--out-selection",type=Path,required=True)
    a=ap.parse_args()

    design=json.loads(DESIGN.read_text())
    fs=design["fcp_source"]
    q=design["qualification"]
    target=set(design["temporal_source"]["clades"])

    if sha256(a.candidate)!=fs["candidate_frame_sha256"]:
        raise RuntimeError("candidate frame SHA256 mismatch")
    if sha256(a.third_selected)!=fs["third_selection_sha256"]:
        raise RuntimeError("third selection SHA256 mismatch")

    selected_before=read_tsv_species(a.third_selected)
    if len(selected_before)!=fs["third_selected_species"]:
        raise RuntimeError(f"unexpected prior selected count {len(selected_before)}")

    rows=[]
    with a.candidate.open(newline="",encoding="utf-8") as f:
        r=csv.DictReader(f)
        fields=set(r.fieldnames or [])
        required={"inat_taxon_id","species","after_observer_cap"}
        if not required.issubset(fields):
            raise RuntimeError(f"candidate frame missing {sorted(required-fields)}; got {sorted(fields)}")
        for row in r:
            tax=str(row["inat_taxon_id"]).strip()
            sp=str(row["species"]).strip()
            cap=int(float(row["after_observer_cap"]))
            rows.append((tax,sp,cap))

    if len(rows)!=fs["candidate_frame_species"]:
        raise RuntimeError(f"candidate frame count {len(rows)} != {fs['candidate_frame_species']}")
    if any(cap<fs["minimum_after_observer_cap"] for _,_,cap in rows):
        raise RuntimeError("candidate frame contains species below frozen U100 gate")

    candidate_keys={(tax,sp) for tax,sp,_ in rows}
    if not selected_before.issubset(candidate_keys):
        raise RuntimeError("third selected manifest is not a subset of candidate frame")

    unused=[x for x in rows if (x[0],x[1]) not in selected_before]
    if len(unused)!=fs["future_unused_pool_expected"]:
        raise RuntimeError(f"unused count {len(unused)} != {fs['future_unused_pool_expected']}")

    by=defaultdict(list)
    unmatched=0
    for tax,sp,cap in unused:
        genus=sp.split()[0] if sp else ""
        if genus in target:
            by[genus].append((tax,sp,cap))
        else:
            unmatched+=1

    counts={clade:len(by.get(clade,[])) for clade in sorted(target)}
    qualifying=[clade for clade,n in counts.items() if n>=q["minimum_unused_species_per_clade"]]

    salt=design["selection"]["salt"]
    selected=[]
    for clade in qualifying:
        ordered=sorted(
            by[clade],
            key=lambda x:(select_hash(salt,x[0],x[1]),int(x[0]))
        )
        for tax,sp,cap in ordered[:q["maximum_selected_species_per_clade"]]:
            selected.append({
                "clade":clade,
                "inat_taxon_id":tax,
                "species":sp,
                "after_observer_cap":cap,
                "selection_hash":select_hash(salt,tax,sp)
            })

    gate=(
        len(qualifying)>=q["minimum_qualifying_clades"]
        and len(selected)>=q["minimum_total_selected_species"]
    )
    status="PASS_PROSPECTIVE_SPATIOTEMPORAL_BRIDGE_OPPORTUNITY" if gate else "HOLD_INSUFFICIENT_FRESH_CLADE_COVERAGE"

    a.out_json.parent.mkdir(parents=True,exist_ok=True)
    a.out_selection.parent.mkdir(parents=True,exist_ok=True)

    with a.out_selection.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=["clade","inat_taxon_id","species","after_observer_cap","selection_hash"],delimiter="\t")
        w.writeheader()
        if gate:
            w.writerows(selected)

    result={
        "version":"v0.1",
        "status":status,
        "design":"data/spatiotemporal_bridge_opportunity_design_v0_1.json",
        "outcome_blind":True,
        "candidate_frame_sha256":fs["candidate_frame_sha256"],
        "prior_third_selection_sha256":fs["third_selection_sha256"],
        "candidate_species":len(rows),
        "prior_third_selected_species":len(selected_before),
        "unused_u100_species":len(unused),
        "target_temporal_clades":len(target),
        "matched_unused_species":sum(counts.values()),
        "unmatched_unused_species":unmatched,
        "unused_species_per_clade":counts,
        "qualifying_clades":qualifying,
        "n_qualifying_clades":len(qualifying),
        "selected_species_if_pass":len(selected) if gate else 0,
        "qualification_rule":q,
        "selection_manifest":str(a.out_selection),
        "selection_manifest_sha256":sha256(a.out_selection),
        "forbidden_outcomes_opened":False,
        "bridge_outcome_computed":False,
        "interpretation":(
            "Fresh high-depth taxonomic coverage is sufficient to preselect a genuinely new cross-scale spatial cohort without using colour or temporal-effect values. The biological bridge outcome remains unopened."
            if gate else
            "The existing untouched U100 pool does not meet the frozen clade-coverage gate. Do not relax the genus mapping or minimum-clade threshold post hoc and do not compute a bridge correlation from the sparse existing cohorts."
        )
    }
    a.out_json.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,indent=2,sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
