#!/usr/bin/env python3
"""Full FCP 42,111-species global flower-photo OPPORTUNITY atlas, no pixels.

Never equate the original >=100-photo U100 gate with the global universe.
Every taxon stays in the exported denominator, including 789 with no eligible
photographs. Original 162-cell photo-discovery metadata and capacity audit are
read from fixed, SHA256-verified historical bytes. No biological colours or
fitness traits are inspected or inferred.

Separate broad colour *observation* opportunity from credible intraspecific
colour diversity or 50km within-species spatial inference.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

SPECIES_SHA = "5d871bd1f1190c68fd513bb527e6c93de1606dc35ff8cfff5f6187850dc382cc"
CAPACITY_SHA = "12d94100d6343597aa2b87670555e801ce0681f4b69e3ad634fa6d80409008e8"
V1_CELL_INDEX_SHA = "b60a8b1b98fcf6745344d6acc394bb6ff905b9d5bded134031402892ff2efa2e"
V2_CELL_INDEX_SHA = "d7f08601ccc4d7dc9ff1943ea3e58ed903b9759277f7b536bc7423bc08189b3a"
FRAME_SIZE = 42111
N_CELLS = 162
CAPACITY_TIERS = (1,2,5,10,20,30,40,50,60,80,100)
CAPACITY_EXPECTED = {
    1:41322, 2:33944, 5:24612, 10:18301, 20:12985,
    30:10330, 40:8753, 50:7632, 60:6770, 80:5558, 100:4730
}
PHOTO_QUOTAS = (1,5,20,40,100)
NATURAL_POLYMORPHISM_CONFIRMED_FROM_METADATA = False


def sha256(path:Path)->str:
    dig=hashlib.sha256()
    with path.open("rb") as h:
        for c in iter(lambda:h.read(1<<20), b""):
            dig.update(c)
    return dig.hexdigest()


def checked(path:Path, fingerprint:str)->None:
    if sha256(path)!=fingerprint:
        raise RuntimeError(f"Immutable historical metadata SHA256 mismatch: {path.name}")


def availability_band(n: int)->str:
    if n==0: return "0_no_metadata_eligible_photos"
    if n==1: return "1_photo"
    if n<5: return "2_to_4_photos"
    if n<10: return "5_to_9_photos"
    if n<20: return "10_to_19_photos"
    if n<40: return "20_to_39_photos"
    if n<100: return "40_to_99_photos"
    return "100_plus_photos"


def load_universe(species_path:Path, capacity_path:Path,
                  links_v1:Path, links_v2:Path)->tuple[pd.DataFrame,dict]:
    for p,hashval in [
        (species_path,SPECIES_SHA),(capacity_path,CAPACITY_SHA),
        (links_v1,V1_CELL_INDEX_SHA),(links_v2,V2_CELL_INDEX_SHA)
    ]:
        checked(p,hashval)

    species=pd.read_csv(species_path,low_memory=False)
    cap=pd.read_csv(capacity_path,low_memory=False)
    if len(species)!=FRAME_SIZE or len(cap)!=FRAME_SIZE:
        raise RuntimeError("Original full 42,111 frame not restored")
    if not {"inat_taxon_id","species"}.issubset(species):
        raise RuntimeError("Source species taxon identities unavailable")
    required={"inat_taxon_id","after_observer_cap","maximum_span_km","request_error"}
    if not required.issubset(cap):
        raise RuntimeError("Original whole-species photographic opportunity fields missing")
    for name,table in [("species",species),("capacity",cap)]:
        if table.inat_taxon_id.isna().any() or table.inat_taxon_id.duplicated().any():
            raise RuntimeError(f"Repeated or missing original taxon ID: {name}")
    if cap.request_error.fillna("").astype(str).str.strip().ne("").any():
        raise RuntimeError("The capacity frame still has unresolved source query failures")

    # Never use any flower-colour outcome to determine who belongs to atlas.
    a=species[["inat_taxon_id","species"]].merge(
        cap[["inat_taxon_id","after_observer_cap","maximum_span_km"]],
        on="inat_taxon_id",how="left",validate="one_to_one")
    if len(a)!=FRAME_SIZE or a.after_observer_cap.isna().any():
        raise RuntimeError("Unresolved capacity record in full world species frame")
    a["after_observer_cap"]=pd.to_numeric(a.after_observer_cap,errors="raise").astype(int)
    if (a.after_observer_cap<0).any():
        raise RuntimeError("Negative observer-capped photo availability not meaningful")
    a["maximum_span_km"]=pd.to_numeric(a.maximum_span_km,errors="coerce")
    if a.species.astype(str).str.strip().eq("").any():
        raise RuntimeError("Blank source taxon name")

    links=[]
    for f in (links_v1,links_v2):
        z=pd.read_csv(f,low_memory=False,usecols=["inat_taxon_id","cell_id"])
        z=z.dropna(subset=["inat_taxon_id","cell_id"])
        z["inat_taxon_id"]=pd.to_numeric(z.inat_taxon_id,errors="raise").astype(int)
        z["cell_id"]=pd.to_numeric(z.cell_id,errors="raise").astype(int)
        if len(z) and (z.cell_id.lt(0)|z.cell_id.ge(N_CELLS)).any():
            raise RuntimeError("Photo-discovery cells outside frozen equal-area world geometry")
        links.append(z)
    all_links=pd.concat(links,ignore_index=True).drop_duplicates(
        ["inat_taxon_id","cell_id"],keep="first")
    full_ids=set(a.inat_taxon_id.astype(int))
    if not set(all_links.inat_taxon_id.astype(int)).issubset(full_ids):
        raise RuntimeError("Cell-geography metadata contain undiscovered taxon IDs")
    if all_links.inat_taxon_id.nunique()!=FRAME_SIZE:
        raise RuntimeError("Some discovered global taxa have no historically recorded geographic cell")
    num=all_links.groupby("inat_taxon_id",sort=False).cell_id.nunique()
    a["discovery_equal_area_cells_162"]=a.inat_taxon_id.map(num).astype(int)
    a["photo_opportunity_band"]=a.after_observer_cap.map(availability_band)
    # Nested planned photo depth: all taxa are retained, even 0-photo species.
    for q in PHOTO_QUOTAS:
        a[f"max_photo_rows_at_cap_{q}"]=np.minimum(a.after_observer_cap,q)
    a["atlas_label_status"]="NOT_MEASURED_FROM_THIS_METADATA"
    a["white_achromatic_state_validated"]=False
    a["genetic_FCP_polymorphism_validated"]=False
    a["observed_no_eligible_photo_metadata"]=a.after_observer_cap.eq(0)
    a["source_only_descriptive_photo_possible"]=a.after_observer_cap.ge(1)
    a["potential_visible_colour_variation_screen"]=a.after_observer_cap.ge(5)
    a["potential_exploratory_colour_diversity"]=a.after_observer_cap.ge(20)
    a["potential_40plus_ITV_measurement"]=a.after_observer_cap.ge(40)
    a["potential_100plus_high_depth_geography"]=a.after_observer_cap.ge(100)
    a=a.sort_values("inat_taxon_id",kind="stable").reset_index(drop=True)

    th={str(i):int(a.after_observer_cap.ge(i).sum()) for i in CAPACITY_TIERS}
    if th!={str(i):n for i,n in CAPACITY_EXPECTED.items()}:
        raise RuntimeError(f"Frozen abundance thresholds no longer reproduce published U100 frame: {th}")
    if len(a.loc[a.after_observer_cap.eq(0)])!=789:
        raise RuntimeError("Previously frozen no-eligible-photo group is not 789 taxa")
    if a.inat_taxon_id.nunique()!=FRAME_SIZE or a.species.nunique()!=FRAME_SIZE:
        raise RuntimeError("FCP global atlas taxa or species-name identity not unique")

    metadata={
        "photo_depth_threshold_counts":th,
        "source_taxa_with_ge1_eligible_photo":int(a.after_observer_cap.ge(1).sum()),
        "source_taxa_with_no_current_eligible_photo":int(a.after_observer_cap.eq(0).sum()),
        "photographic_preprocessing_plan_by_max_per_species":{
            str(q):{
                "number_of_species_with_at_least_one_photo":int(a.after_observer_cap.gt(0).sum()),
                "sum_of_possible_nested_photo_rows":int(a[f"max_photo_rows_at_cap_{q}"].sum()),
                "species_reaching_all_requested_rows":int(a.after_observer_cap.ge(q).sum())
            } for q in PHOTO_QUOTAS
        },
        "species_geographic_cell_counts":{
            "min":int(a.discovery_equal_area_cells_162.min()),
            "median":float(a.discovery_equal_area_cells_162.median()),
            "max":int(a.discovery_equal_area_cells_162.max()),
            "n_2plus_cells":int(a.discovery_equal_area_cells_162.ge(2).sum())
        },
        "n_world_equal_area_cells":N_CELLS
    }
    return a,metadata


def main()->None:
    p=argparse.ArgumentParser()
    p.add_argument("--species",type=Path,required=True)
    p.add_argument("--capacity",type=Path,required=True)
    p.add_argument("--links-v1",type=Path,required=True)
    p.add_argument("--links-v2",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    args=p.parse_args()
    a,stats=load_universe(args.species,args.capacity,args.links_v1,args.links_v2)
    args.outdir.mkdir(parents=True,exist_ok=True)
    fname=args.outdir/"fcp_all42111_species_opportunity_ledger.csv.gz"
    a.to_csv(fname,index=False,lineterminator="\n",
             compression={"method":"gzip","compresslevel":9,"mtime":0})
    totals=a.groupby("photo_opportunity_band",sort=True,observed=True).agg(
        n_species=("inat_taxon_id","size"),
        n_species_with_2plus_geocells=("discovery_equal_area_cells_162",
                                       lambda x:int((x>=2).sum())),
        med_geographic_cells=("discovery_equal_area_cells_162","median"),
        historical_photo_capacity_sum=("after_observer_cap","sum")
    ).reset_index()
    totals.to_csv(args.outdir/"whole_42111_photo_opportunity_bands.csv",index=False)
    report={
        "schema":"fcp_complete_42111_photo_opportunity_nested_inference_v1",
        "date_jst":"2026-10-08",
        "status":"COMPLETE_WORLD_SPECIES_FRAME_METADATA_ONLY_NO_IMAGE_PIXELS",
        "full_species_discovery_frame":FRAME_SIZE,
        "source_sha256":{
            "species":SPECIES_SHA,"capacity":CAPACITY_SHA,
            "photo_discovery_index_v1":V1_CELL_INDEX_SHA,
            "photo_discovery_index_v2":V2_CELL_INDEX_SHA
        },
        "full_atlas_species_coverage":"ALL 42111, INCLUDING NO-ELIGIBLE-PHOTO",
        "photo_measurement_started_by_this_workflow":False,
        "colour_and_fitness_outcome_opened":False,
        "nested_photo_depth_opportunity":stats,
        "recorded_full_species_ledger_SHA256":sha256(fname),
        "observation_scope_interpretation":{
            "n_ge1":"species with any historical photo opportunity; image colour not necessarily classifiable",
            "n_ge5":"screened visible colour occurrence only; failure to see rare morph not natural monomorphism",
            "n_ge20":"comparative visible-colour diversity has broad but imperfect photographic support",
            "n_ge40":"candidate species for minimum measured within-species ITV",
            "n_ge100":"candidate species for higher-depth geographic colour association, not necessarily enough local photo pairs"
        },
        "claim_limits":[
            "metadata-only observer-capped photo availability is not actual measured flower-colour photos",
            "zero-photo species remain in the denominator as NOT OBSERVED, not white/monomorphic",
            "single-photo species yield no within-species colour variation estimates",
            "two or five photo labels cannot reliably prove the absence of rare flower colour variants",
            "presence in >1 equal-area cell is not proof of photographed reproductive populations or ecological adaptation",
            "the 42111 species are an outcome-blind iNaturalist discovered opportunity frame, not all world angiosperms",
            "primary benefit-cost/selection mechanisms, genotypes and white anthocyanin absence are not measured"
        ],
        "confirmed_global_photo_colour_prevalence_estimable_now":False,
        "confirmatory_decisions_changed":False
    }
    (args.outdir/"result.json").write_text(
        json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+"\n",
        encoding="utf-8")
    print(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False))


if __name__=="__main__":
    main()
