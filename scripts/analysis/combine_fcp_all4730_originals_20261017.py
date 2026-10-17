#!/usr/bin/env python3
"""Final exact-original FCP source extension to the entire U100=4730 species frame.

The old 149900 photographic rows and frozen H1, 50-km, H2/W science functions
are retained. Wave4 and final all-U100 new actual pixel-terminal records are
strictly added; nonclassifiable photos remain valid terminal photo records,
while species with <100 authorized H9 photos are explicitly not invented.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from scripts.analysis.append_wave4_to_original_149900_20261017 import (
    OLD_SOURCE, enforce_and_tag, _check_separately_sampled,
)
from scripts.analysis.replay_fcp_149900_plus_wave4_same_spatial_20261017 import (
    ORIGINAL_EXPECTED, original_eligible,
)
from scripts.analysis.run_polymorphism_distributed_polymorphism_posthoc_20261007 import (
    RADII_KM, analyse_species, cohort_summary,
)

COHORTS=("discovery","validation","prospective","wave4","new_U100")
MAX_TAXA=4730


def sha(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda:handle.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def authenticate_terminal(path:Path, measured:Path, *, wave4:bool)->tuple[int,int]:
    r=json.loads(path.read_text())
    if wave4:
        if (r.get("schema")!="fcp_same_149900_method_wave4_real_photo_terminals_v1"
            or r.get("original_results_changed") is not False
            or r.get("measurement_original_fixed_colour_worker") is not True
            or r.get("terminal_partitions")!=256
            or r.get("measured_sha256")!=sha(measured)):
            raise RuntimeError("wave4 original classifier 256-partition SHA gate failed")
        return int(r["wave4_100photo_species"]),int(r["wave4_new_photo_ids"])
    if (r.get("schema")!="fcp_all4730_new_original_method_pixel_terminal_v1"
        or r.get("original_roi_v4_and_efficientsam_unchanged") is not True
        or r.get("original_frozen_149900_results_changed") is not False
        or r.get("measurement_terminal_partitions")!=256
        or r.get("measured_photo_file_sha256")!=sha(measured)):
        raise RuntimeError("all4730 new original ROI pixel terminal SHA gate failed")
    return int(r["new_complete_100_photo_species"]),int(r["measurement_photo_rows"])


def combine(
    old_paths:dict[str,Path],wave4:Path,wave4_terminal:Path,new:Path,
    new_terminal:Path,all_U100:Path,metadata_receipt:Path,
    *,
    original_contract:dict=OLD_SOURCE,
) ->tuple[pd.DataFrame,dict]:
    if set(old_paths)!=set(OLD_SOURCE):
        raise RuntimeError("all 3 authentic original measured photo cohorts mandatory")
    frames=[];fingerprints={};audits=[]
    photos:set[str]=set();taxa:set[int]=set();species:set[str]=set()
    for name in COHORTS[:3]:
        p=old_paths[name]
        if sha(p)!=original_contract[name]["sha256"]:
            raise RuntimeError("original unmodified 149900 input SHA drift: "+name)
        d,a=enforce_and_tag(pd.read_csv(p,low_memory=False),label=name,
            rows=original_contract[name]["rows"],species=original_contract[name]["species"],target=100)
        _check_separately_sampled(d,photos,taxa,species,name)
        frames.append(d);audits.append(a);fingerprints[name]=sha(p)
    if sum(len(f) for f in frames)!=149900 or sum(a["named_species"] for a in audits)!=1499:
        raise RuntimeError("old original 149900 photographic denominator drift")
    for name,path,terminal in (("wave4",wave4,wave4_terminal),("new_U100",new,new_terminal)):
        num_taxa,num_photos=authenticate_terminal(terminal,path,wave4=name=="wave4")
        if num_photos!=num_taxa*100:
            raise RuntimeError("not exactly 100 authentic new image terminals per selected species")
        if name=="wave4" and (num_taxa,num_photos)!=(498,49800):
            raise RuntimeError("frozen wave4 498 original photo taxa missing")
        if name=="new_U100" and not 1<=num_taxa<=2733:
            raise RuntimeError("all4730 new U100 source support outside authentic range")
        d,a=enforce_and_tag(pd.read_csv(path,low_memory=False),label=name,
                             rows=num_photos,species=num_taxa,target=100)
        _check_separately_sampled(d,photos,taxa,species,name)
        frames.append(d);audits.append(a);fingerprints[name]=sha(path)
    full=pd.read_csv(all_U100)
    source_receipt=json.loads(metadata_receipt.read_text())
    if (source_receipt.get("schema")!="fcp_u100_all4730_original_h9_only_v1"
        or source_receipt.get("stage")!="all_4730_H9_metadata_terminal"
        or source_receipt.get("new_attempted_species")!=2733
        or source_receipt.get("U100_identity_species")!=4730
        or source_receipt.get("old_149900_method_and_results_changed") is not False
        or source_receipt.get("outputs_sha256",{}).get(all_U100.name)!=sha(all_U100)):
        raise RuntimeError("ALL4730 original H9 source ledger not authenticated")
    if len(full)!=MAX_TAXA or full.inat_taxon_id.nunique()!=MAX_TAXA:
        raise RuntimeError("4730 U100 original identities not complete")
    union=set(full.inat_taxon_id.astype(int))
    if taxa-union:
        raise RuntimeError("measured original taxon not in authentic U100")
    if len(taxa)!=1499+498+source_receipt["new_complete_100_photo_species"]:
        raise RuntimeError("all H9-complete species not pixel-measured exactly once")
    if len(taxa)+int(source_receipt["new_unmet_100_photo_species"])!=MAX_TAXA:
        raise RuntimeError("full 4730 taxon denominator was inflated or lost")
    for f in frames:
        pair=f[["inat_taxon_id","species"]].drop_duplicates()
        expected=full.set_index("inat_taxon_id").species
        for tid,name in pair.itertuples(index=False):
            if expected.loc[tid]!=name:raise RuntimeError("historical taxon names changed on extension")
    merged=pd.concat(frames,ignore_index=True,sort=False)
    if len(merged)!=149900+49800+source_receipt["new_original_photo_IDs"]:
        raise RuntimeError("all source photographic terminal row count incomplete")
    if merged.photo_id.nunique()!=len(merged) or merged.inat_taxon_id.nunique()!=len(taxa):
        raise RuntimeError("duplicate or missing original image/photo identity")
    result={
      "schema":"fcp_all4730_original_149900_plus_all_new_pixel_measurements_v1",
      "full_original_eligible_U100_taxa":4730,
      "total_original_100photo_measured_species":len(taxa),
      "total_original_pixel_terminal_photo_rows":len(merged),
      "remaining_U100_under_100_original_H9_photos":int(source_receipt["new_unmet_100_photo_species"]),
      "original_149900_photo_rows_unchanged":True,
      "old_50km_and_W_estimators_unchanged":True,
      "new_original_photo_rows":len(merged)-149900,
      "wave4_authentic_photo_rows":49800,
      "final_H9_additional_photo_rows":source_receipt["new_original_photo_IDs"],
      "cohorts":audits,
      "source_sha256":fingerprints,
      "all_U100_species_status_sha256":sha(all_U100),
      "all_U100_metadata_receipt_sha256":sha(metadata_receipt),
      "no_unmeasured_species_imputed":True,
      "same_original_iNaturalist_source_not_independent_provider":True,
    }
    return merged,result


def replay_same_old_50km(photos:pd.DataFrame)->tuple[pd.DataFrame,dict]:
    results=[];nulls={r:[] for r in RADII_KM}
    cohorts={}
    for cohort in COHORTS:
        group=original_eligible(photos.loc[photos.fcp_source_cohort.eq(cohort)])
        vals=[];cohort_null={r:[] for r in RADII_KM}
        for _,g in group.groupby("inat_taxon_id",sort=True):
            g=g.dropna(subset=["latitude","longitude"])
            if len(g)<40:continue
            seed_label="third" if cohort in ("prospective","wave4","new_U100") else cohort
            item,permutations=analyse_species(g,seed_label)
            item["source_cohort"]=cohort
            vals.append(item)
            for radius in RADII_KM:
                nulls[radius].append(permutations[radius])
                cohort_null[radius].append(permutations[radius])
        if not vals:
            raise RuntimeError("no spatially eligible originally measured species in "+cohort)
        d=pd.DataFrame(vals)
        summary=cohort_summary(d,cohort_null)
        if cohort in ORIGINAL_EXPECTED:
            n,effect,p=ORIGINAL_EXPECTED[cohort]
            old=summary["50"]
            if (old["n_species"]!=n or abs(old["mean_depletion"]-effect)>1e-10
                or old["mean_depletion_p_upper"]!=p):
                raise RuntimeError("original 149900 50km frozen statistic changed: "+cohort)
        cohorts[cohort]=summary
        results.extend(vals)
    out=pd.DataFrame(results)
    combined=cohort_summary(out,nulls)
    record={
       "schema":"fcp_all4730_same_original_50km_vertex_null_v1",
       "old_three_original_149900_effects_exactly_reproduced":True,
       "original_50km_radius_changed":False,
       "original_species_equal_weighting_changed":False,
       "original_199_vertex_null_changed":False,
       "all_sources":cohorts,
       "pooled_original_50km":combined["50"],
       "all_evaluable_species":int(out.inat_taxon_id.nunique()),
       "cohort_distribution":{key:int(out.source_cohort.eq(key).sum()) for key in COHORTS},
       "still_posthoc_50km":True,
       "photographic_observation_biological_population_frequency_not_assumed":True,
    }
    return out,record


def main()->None:
    parser=argparse.ArgumentParser()
    sub=parser.add_subparsers(dest="action",required=True)
    a=sub.add_parser("combine")
    for name in COHORTS[:3]:
        a.add_argument("--"+name,type=Path,required=True)
    for opt in ("wave4","wave4-terminal","new","new-terminal","all-U100","metadata-receipt","out"):
        a.add_argument("--"+opt,type=Path,required=True)
    b=sub.add_parser("50km")
    b.add_argument("--combined",type=Path,required=True)
    b.add_argument("--receipt",type=Path,required=True)
    b.add_argument("--out",type=Path,required=True)
    c=sub.add_parser("H1-D-W")
    c.add_argument("--combined",type=Path,required=True)
    c.add_argument("--receipt",type=Path,required=True)
    c.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    if args.action=="combine":
        data,r=combine(
           {key:getattr(args,key) for key in COHORTS[:3]},
           args.wave4,args.wave4_terminal,args.new,args.new_terminal,
           getattr(args,"all_U100"),args.metadata_receipt)
        p=args.out/"fcp_all4730_original_method_measured_photos.csv.gz"
        data.to_csv(p,index=False,lineterminator="\n",
                 compression={"method":"gzip","compresslevel":9,"mtime":0})
        r["merged_photo_table_sha256"]=sha(p)
        (args.out/"fcp_all4730_original_photo_combine_receipt.json").write_text(
            json.dumps(r,indent=2,sort_keys=True,allow_nan=False)+"\n")
    else:
        r=json.loads(args.receipt.read_text())
        if (r.get("schema")!="fcp_all4730_original_149900_plus_all_new_pixel_measurements_v1"
            or r.get("merged_photo_table_sha256")!=sha(args.combined)):
            raise RuntimeError("total 4730 original photo table not authenticated")
        if args.action=="50km":
            species,metric=replay_same_old_50km(pd.read_csv(args.combined,low_memory=False))
            species.to_csv(args.out/"fcp_all4730_unchanged_original_50km_species.csv",index=False)
            (args.out/"fcp_all4730_unchanged_original_50km_results.json").write_text(
               json.dumps(metric,indent=2,sort_keys=True,allow_nan=False)+"\n")
            r=metric
        else:
            # Only the H1/D/W job restores the historical 149900-photo modules.
            # Combine and original 50km jobs do not need legacy colour/ROI imports.
            from scripts.analysis.run_fcp_wave4_original_roi_and_statistics_20261017 import (
                analyze_original_methods,
            )
            # This is exactly the existing high-depth historical H1/D/H2-W code.
            original=analyze_original_methods(args.combined,args.out)
            r={"schema":"fcp_all4730_unchanged_historical_H1_D_H2_W_v1",
               "all_real_photo_rows":r["total_original_pixel_terminal_photo_rows"],
               "existing_method_called_directly":True,
               "old_frozen_149900_inference_not_overwritten":True,
               "H1_original_observer_disjoint":original["h1_original_method"],
               "H2_W_original_fixed_white_axis":original["h2_original_W_method"],
               "same_iNaturalist_source":True}
            (args.out/"fcp_all4730_H1_D_W_results.json").write_text(
               json.dumps(r,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("FCP_ALL4730_ORIGINAL_"+args.action.upper()+"_DONE",json.dumps(r,sort_keys=True,allow_nan=False))


if __name__=="__main__":
    main()
