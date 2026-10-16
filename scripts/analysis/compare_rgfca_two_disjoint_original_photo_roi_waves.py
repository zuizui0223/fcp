#!/usr/bin/env python3
"""Compare exactly two DISJOINT six-taxon genuine ROI-v4 technical waves.

Requires independent source/measurement SHA-bound outputs from
audit_rgfca_real_photo_observer_geocell_mixture.py. Both source wave photo IDs
originated from a LIVE iNaturalist API with a fixed created-before date,
not confirmed membership in the official 2026-09-27 monthly archive.

This is quantitative measurement-bias calibration: no heritable flower
polymorphism, shared geography, W, H2 or causal environmental effect inference.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from collections import Counter,defaultdict
from pathlib import Path

# Original RGFCA 1000 selected-species CSV SHA from independent metadata-only CI.
ROSTER_SHA = "88aba565e1bdeef48ad3dc5eec0121c52448fda434aca849e542bb615f9d547c"
DEPTHS=(20,32,40,44,48)
COLOURS=("white","yellow_orange","red_pink","blue_purple")


def sha256(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):
            h.update(b)
    return h.hexdigest()


def read_csv(path:Path)->list[dict]:
    with path.open("r",newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))


def read_wave(root:Path, wave:int, *, synthetic:bool=False)->list[dict]:
    report=root/f"wave{wave}_original_ROI_observer_geocell_mixed_colour_QC.json"
    tab=root/f"wave{wave}_original_ROI_observer_geocell_mixed_colour_QC.csv"
    meta=json.loads(report.read_text(encoding="utf-8"))
    if (meta.get("schema")!="fcp_real_rgfca_original_photo_technical_observer_geo_organ_QC_v1"
        or meta.get("status")!="TECHNICAL_REAL_ORIGINAL_ROI_ONLY_NOT_POLYMORPHISM_ESTIMATE"
        or meta.get("pilot_wave")!=wave
        or meta.get("n_taxa")!=6
        or meta.get("source_pilot_random_photo_draws_per_taxon")!=200
        or meta.get("pilot_resampling_fold_depths")!=list(DEPTHS)
        or meta.get("CSV_sha256")!=sha256(tab)
        or meta.get("original_source_date_bulk_archive_membership_verified") is not False
        or meta.get("main_1000_RGFCA_P_poly_estimated") is not False
        or meta.get("no_genetic_FCP_species_asserted") is not True
        or meta.get("no_old_colour_measurements_reused") is not True):
        raise ValueError(f"Original RGFCA wave {wave} observer/site ROI label receipt not admitted")
    rows=read_csv(tab)
    if len(rows)!=6:
        raise ValueError("Each technical-only RGFCA wave requires six actual source taxa")
    total=valid=0
    seen=set()
    for row in rows:
        tid=row["inat_taxon_id"]
        if not tid.isdecimal() or tid in seen or not row["species"].strip():
            raise ValueError("Bad original technical source taxon identity")
        seen.add(tid)
        n,v=int(row["source_photo_ids_selected"]),int(row["classifiable_ROI"])
        if n<0 or not 0<=v<=n or row["pilot_wave"]!=str(wave):
            raise ValueError("Invalid or missing original ROI source terminal counts")
        counts=json.loads(row["photo_label_counts"])
        if (set(counts)-set(COLOURS)
            or sum(counts.values())!=v
            or int(row["n_classifiable_original_observers"])>v
            or int(row["n_classifiable_5degree_cells"])>v
            or row["genetic_or_population_colour_polymorphism_confirmed"]!="False"):
            raise ValueError("Original ROI floral colours/cell/observer audit contradicted")
        for dep in DEPTHS:
            prefix=f"depth_{dep}_"
            eligible=int(row[prefix+"evaluable_repetitions"])
            pos=int(row[prefix+"positive_photo_label_repetitions"])
            ratio=row[prefix+"conditional_photo_label_detection"]
            if not 0<=pos<=eligible<=200:
                raise ValueError("Original technical photo Monte Carlo sample counts invalid")
            if (eligible==0 and ratio!=""
                or eligible and abs(float(ratio)-pos/eligible)>1e-8):
                raise ValueError("Repeated observer-disjoint photo-label denominator drifted")
        total+=n;valid+=v
    if total!=meta["n_original_ROI_photos"] or valid!=meta["n_classifiable_flower_ROIs"]:
        raise ValueError("Original ROI per-species terminal counts do not reconcile to wave receipt")
    return rows


def compare(wave1:Path,wave2:Path,roster:Path,outdir:Path,*,synthetic=False)->dict:
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError("Cannot overwrite original two-wave photo QC")
    if not synthetic and sha256(roster)!=ROSTER_SHA:
        raise ValueError("Original 1000 RGFCA source-only study enrolment roster drifted")
    r=read_csv(roster)
    registry={x["inat_taxon_id"]:x for x in r}
    if len(registry)!=len(r):
        raise ValueError("Duplicate taxa in original balanced RGFCA study frame")
    one=read_wave(wave1,1,synthetic=synthetic)
    two=read_wave(wave2,2,synthetic=synthetic)
    ids=set(x["inat_taxon_id"] for x in one)
    if len(ids&{x["inat_taxon_id"] for x in two}):
        raise ValueError("Wave2 biological species are not independent of first six")
    joined=[]
    tiers=defaultdict(lambda:Counter())
    for row in one+two:
        key=row["inat_taxon_id"]
        if key not in registry or registry[key]["species"]!=row["species"]:
            raise ValueError("Photo pilot taxon not in immutable random 1000 enrolment")
        tier=registry[key]["opportunity_stratum"]
        counts=json.loads(row["photo_label_counts"])
        n,v=int(row["source_photo_ids_selected"]),int(row["classifiable_ROI"])
        structural=bool(
            v>=40 and counts.get("white",0)>0 and counts.get("yellow_orange",0)>0
            and float(row["within_one_ROI_white_and_yellow_fraction"])>=.50)
        capable=int(row["depth_48_evaluable_repetitions"])>=30
        photo_signal=int(row["depth_48_positive_photo_label_repetitions"])>0
        record={
            "pilot_wave":int(row["pilot_wave"]),
            "inat_taxon_id":int(key),"species":row["species"],
            "original_source_photo_count":int(registry[key]["source_photo_opportunities"]),
            "pre_colour_photo_opportunity_tier":tier,
            "original_pixel_terminals":n,"classifiable_original_ROI":v,
            "classifiable_rate":round(v/n,8) if n else 0.0,
            "classified_distinct_observers":int(row["n_classifiable_original_observers"]),
            "classified_geo_5degree_cells":int(row["n_classifiable_5degree_cells"]),
            "ROI_white_yellow_dual_fraction":float(row["within_one_ROI_white_and_yellow_fraction"]),
            "RGFCA_48_per_fold_evaluable_repetitions":int(row["depth_48_evaluable_repetitions"]),
            "RGFCA_48_per_fold_photo_LABEL_detection":(
                row["depth_48_conditional_photo_label_detection"]),
            "photolabel_repetition_support_not_biological_FCP":capable and photo_signal,
            "within_one_ROI_structural_review_hold":structural,
            "phenotypic_polymorphism_confirmed":False,
            "genetic_polymorphism_confirmed":False,
        }
        joined.append(record)
        tierstats=tiers[tier]
        tierstats["n_taxa"]+=1
        tierstats["n_original_photos"]+=n
        tierstats["n_classifiable_ROI"]+=v
        tierstats["n_taxa_ge40_valid"]+=int(v>=40)
        tierstats["n_taxa_structural_review_hold"]+=int(structural)
    if not synthetic and (
        len(joined)!=12
        or set(tiers)!={"shallow_128_255","medium_256_1023","deep_ge1024"}
        or any(x["n_taxa"]!=4 for x in tiers.values())):
        raise ValueError("Original 4-per-depth-tier 12-species technical cohort changed")
    outdir.mkdir(parents=True,exist_ok=True)
    outcsv=outdir/"rgfca_12_real_technical_taxa_photo_roi_observer_bias.csv"
    with outcsv.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(joined[0]),lineterminator="\n")
        w.writeheader()
        w.writerows(sorted(joined,key=lambda v:(v["pilot_wave"],v["inat_taxon_id"])))
    nphoto=sum(x["original_pixel_terminals"] for x in joined)
    nvalid=sum(x["classifiable_original_ROI"] for x in joined)
    doc={
        "schema":"fcp_rgfca_two_six_taxa_original_roi_photo_bias_calibration_v1",
        "status":"REAL_PHOTO_PIXELS_TWO_TECHNICAL_WAVES_NOT_ORIGINAL_DATED_BIOLOGY",
        "n_original_source_taxa_pilot":len(joined),
        "n_distinct_technical_waves":2,
        "source_taxon_overlaps_between_waves":0,
        "all_taxa_preselected_from_original_random_1000":True,
        "total_original_new_image_pixels_ROI_measured":nphoto,
        "total_classifiable_original_flower_ROIs":nvalid,
        "image_classifiable_fraction_descriptive_only":nvalid/nphoto if nphoto else 0,
        "n_species_with_at_least_40_classifiable_ROIs":sum(
            x["classifiable_original_ROI"]>=40 for x in joined),
        "n_species_flagged_for_flower_organ_structure_review":sum(
            x["within_one_ROI_structural_review_hold"] for x in joined),
        "available_photo_opportunity_strata":{
            tier:{**dict(count), "photo_ROI_classifiable_rate": (
                count["n_classifiable_ROI"]/count["n_original_photos"]
                if count["n_original_photos"] else 0)}
            for tier,count in sorted(tiers.items())
        },
        "real_photo_QC_table_sha256":sha256(outcsv),
        "same_platform_and_live_API_source_NOT_independent_dated_snapshot":True,
        "actual_20260927_original_source_photo_ID_membership_verified":False,
        "actual_1000_taxon_200x250_formal_Ppoly_estimated":False,
        "any_natural_population_flower_colour_polymorphism_proven":False,
        "biological_monomorphism_assumed_for_non_detected_taxa":False,
        "repeated_photo_resampling_does_not_expand_biological_sample":True,
    }
    (outdir/"rgfca_12_real_technical_taxa_photo_roi_observer_bias_receipt.json").write_text(
        json.dumps(doc,indent=2,ensure_ascii=False,sort_keys=True)+"\n",encoding="utf-8")
    return doc


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--wave1-qc",required=True,type=Path)
    p.add_argument("--wave2-qc",required=True,type=Path)
    p.add_argument("--original-1000-roster",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    x=p.parse_args(argv)
    print("FCP_RGFCA_12_REAL_TECHNICAL_PIXELS",
          json.dumps(compare(x.wave1_qc,x.wave2_qc,x.original_1000_roster,x.outdir),sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
