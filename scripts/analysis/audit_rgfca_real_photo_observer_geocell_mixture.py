#!/usr/bin/env python3
"""Independent observer/geocell and within-ROI colour QC on LIVE-API FCP pilots.

Strictly technical, never production 2026-09-27 archive FCP:
  * only authentic newly measured original ROI-v4/EfficientSAM photo IDs;
  * matches every photo ID/source observer/species and all receipt hashes;
  * reports classification/observer/geography/flower-organ colour mixing;
  * 2+ independent observers or 2+ geocells do NOT establish local FCP.

Works on the first and SECOND DISJOINT six-taxon outcome-blind technical waves.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

from scripts.analysis.audit_rgfca_real_six_taxa_photo_repeatability_20261016 import (
    DEPTHS, repeat_descriptive,
)

MORPHS=("white","yellow_orange","red_pink","blue_purple")
PALETTE=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
ORIGINAL_ROSTER_SHA="88aba565e1bdeef48ad3dc5eec0121c52448fda434aca849e542bb615f9d547c"
ROI_SHA={
    "detector_sha256":"f1aaeec4664fe2c178e5cf2bc1f508977bef3e4aa7b40613026cb8ae3de789d5",
    "efficient_sam_encoder_sha256":"84ed466ffcc5c1f8d08409bc34a23bb364ab2c15e402cb12d4335a42be0e0951",
    "efficient_sam_decoder_sha256":"a62f8fa5ea080447c0689418d69e58f1e83e0b7adf9c142e2bd9bcc8045c0b11",
}
API_KIND="iNaturalist_API_created_before_20260927_live_snapshot_NOT_dated_monthly_archive"
VALID_SOURCE_SCHEMAS={
    "fcp_rgfca_1000_enrolled_real_original_API_six_taxon_pilot_v1":1,
    "fcp_rgfca_1000_enrolled_real_original_API_disjoint_second_six_taxon_pilot_v1":2,
}
VALID_ROI_SCHEMAS={
    "fcp_rgfca_new_original_ROI_v4_live_API_technical_pilot_v1":1,
    "fcp_rgfca_new_original_ROI_v4_live_API_technical_second_disjoint_pilot_v1":2,
}
ROI_STATUSES={
    1:"REAL_NEW_ORIGINAL_PIXELS_TECHNICAL_PILOT_NOT_RGFCA_PPOLY",
    2:"REAL_NEW_ORIGINAL_PIXELS_DISJOINT_WAVE2_TECHNICAL_NOT_RGFCA_PPOLY",
}


def sha256(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def csv_rows(path:Path)->list[dict]:
    with path.open("r",newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))


def summarize(pilot_root:Path, roster:Path, outdir:Path, *, synthetic=False)->dict:
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError("Cannot overwrite real source-to-ROI observer/geocell QC")
    def expect_one(glob):
        paths=list(pilot_root.rglob(glob))
        if len(paths)!=1:
            raise ValueError(f"Expected one authentic technical source file: {glob}")
        return paths[0]
    source_path=expect_one("rgfca_six*_source_photo_candidates_before_colour.csv")
    source_doc_path=expect_one("rgfca_six*_source_precolour_receipt.json")
    measured_path=expect_one("rgfca_original_new_photo_six_taxa*live_API_technical.csv")
    roi_doc_path=expect_one("rgfca_original_new_photo_six_taxa*live_API_technical_receipt.json")
    source_doc=json.loads(source_doc_path.read_text(encoding="utf-8"))
    roi_doc=json.loads(roi_doc_path.read_text(encoding="utf-8"))
    wave=VALID_SOURCE_SCHEMAS.get(source_doc.get("schema"))
    if wave is None or VALID_ROI_SCHEMAS.get(roi_doc.get("schema"))!=wave:
        raise ValueError("Technical source/ROI waves mixed or biology receipt overwritten")
    if (source_doc.get("status")!="SELECTED_ORIGINAL_PROVIDER_API_PHOTO_IDS_NO_IMAGE_PIXELS_YET"
        or source_doc.get("NOT_THE_IMMUTABLE_20260927_METADATA_CENSUS") is not True
        or source_doc.get("EXACT_20260927_ARCHIVE_PHOTO_ID_MEMBERSHIP_VERIFIED") is not False
        or source_doc.get("photo_selection_used_literature_phenotype") is not False
        or source_doc.get("photo_selection_used_climate_soil_or_prior_photo_colours") is not False
        or source_doc.get("image_bytes_opened") is not False
        or source_doc.get("RGFCA_200_realization_photo_polymorphism_Ppoly_computed") is not False
        or source_doc.get("source_ROSTER_sha256")!=sha256(roster)
        or source_doc.get("selected_metadata_sha256")!=sha256(source_path)
        or roi_doc.get("status")!=ROI_STATUSES[wave]
        or roi_doc.get("source_kind")!=API_KIND
        or roi_doc.get("exact_dated_20260927_source_membership_admitted") is not False
        or roi_doc.get("full_1000_species_ROI_study_complete") is not False
        or roi_doc.get("prior_FCP_measured_colours_loaded") is not False
        or roi_doc.get("new_terminal_photo_csv_sha256")!=sha256(measured_path)
        or roi_doc.get("source_image_ID_selection_sha256")!=sha256(source_path)
        or roi_doc.get("source_API_selection_receipt_sha256")!=sha256(source_doc_path)
        or source_doc.get("pilot_taxa")!=6
        or roi_doc.get("number_taxa_selected")!=6
        or (wave==2 and (source_doc.get("technical_only_disjoint_previous_six_taxa") is not True
                        or roi_doc.get("technical_only_disjoint_previous_six_taxa") is not True))):
        raise ValueError("Technical source/ROI byte identity or negative FCP claims missing")
    method=roi_doc.get("model_attestation") or {}
    if (not synthetic and
        (sha256(roster)!=ORIGINAL_ROSTER_SHA
         or any(method.get(k)!=v for k,v in ROI_SHA.items())
         or method.get("location_blind") is not True
         or method.get("no_old_measurement_csv_inputs") is not True)):
        raise ValueError("Original RGFCA roster or original ROI-v4/SAM model pin is not valid")
    roster_data={r["inat_taxon_id"]:r for r in csv_rows(roster)}
    source=csv_rows(source_path)
    real=csv_rows(measured_path)
    lookup={(r["inat_taxon_id"],r["photo_id"]):r for r in source}
    if (len(lookup)!=len(source)
        or len(real)!=len(source)
        or roi_doc.get("selected_real_original_photo_IDs")!=len(real)
        or roi_doc.get("original_pixels_ROI_v4_evaluated")!=len(real)):
        raise ValueError("Original selected photo ID ledger is incomplete or duplicate")
    done=set()
    groups=defaultdict(list)
    for row in real:
        key=(row["inat_taxon_id"],row["photo_id"])
        if key in done or key not in lookup:
            raise ValueError("Repeated or replaced original photo ID after source blind selection")
        done.add(key)
        original=lookup[key]
        for col in ("species","observation_id","observer_id","screen_fold",
                    "latitude","longitude","photo_url_large","photo_license"):
            if row[col]!=original[col]:
                raise ValueError("Original photo attribution changed across blind ROI")
        tid,sp=key[0],row["species"]
        if tid not in roster_data or roster_data[tid]["species"]!=sp:
            raise ValueError("Taxon outside fixed outcome-blind 1000 RGFCA enrolment")
        good=row["global_classifiable"]=="true"
        if (row["global_classifiable"] not in ("true","false")
            or good!=(row["measurement_status"]=="classified_four_state_morph"
                      and row["morph"] in MORPHS)
            or (not good and row["morph"]!="mixed_uncertain")):
            raise ValueError("Actual ROI-v4 four-state QC malformed")
        try:
            lat,lon=float(row["latitude"]),float(row["longitude"])
        except (ValueError,TypeError) as e:
            raise ValueError("Original photo coordinates invalid") from e
        if not math.isfinite(lat) or not math.isfinite(lon):
            raise ValueError("Original coordinates not finite")
        palette=None
        if good:
            palette={c:float(row[f"flower_fraction_{c}"]) for c in PALETTE}
            if any(not math.isfinite(x) or x<0 or x>1 for x in palette.values()):
                raise ValueError("Measured ROI-v4 colour palette fractions corrupt")
            if abs(sum(palette.values())-1)>1e-7:
                raise ValueError("Original nine-palette ROI fractions not conserved")
        groups[(int(tid),sp)].append({
            "photo_id":row["photo_id"], "good":good,"morph":row["morph"],
            "observer":row["observer_id"],
            "observation":row["observation_id"],"fold":str(row["screen_fold"]),
            "cell":(math.floor(lat/5),math.floor(lon/5)),"palette":palette,
        })
    if len(groups)!=6 or len(done)!=len(source):
        raise ValueError("Pilot six taxa or selected original photo IDs missing")
    data=[]
    for (tid,sp),rows in sorted(groups.items()):
        usable=[x for x in rows if x["good"]]
        colour_counts=Counter(x["morph"] for x in usable)
        obs_by={m:set() for m in MORPHS}
        watchers={m:set() for m in MORPHS}
        cells={m:set() for m in MORPHS}
        mixed={m:0 for m in MORPHS}
        folds={m:set() for m in MORPHS}
        for x in usable:
            k=x["morph"]
            obs_by[k].add(x["observation"])
            watchers[k].add(x["observer"])
            cells[k].add(x["cell"])
            folds[k].add(x["fold"])
            if x["palette"]["white"]>=0.1 and x["palette"]["yellow"]>=0.1:
                mixed[k]+=1
        comp_total=sum(mixed.values())
        # The same photo URLs and ROI pixel labels are NOT downloaded again.
        # Repeated draws are technical reweighting of the fixed source photos,
        # never 200 independent natural populations or a biological FCP test.
        repeated=repeat_descriptive(tid,sp,[{
            "photo_id":r["photo_id"],"fold":int(r["fold"]),
            "observation_id":r["observation"],"observer_id":r["observer"],
            "classifiable":r["good"],
            "morph":r["morph"] if r["good"] else None,
            "fractions":r["palette"] if r["good"] else {},
        } for r in rows])
        repeats={
            f"depth_{n}_{field}":repeated[f"depth_{n}_{field}"]
            for n in DEPTHS
            for field in (
                "evaluable_repetitions",
                "positive_photo_label_repetitions",
                "conditional_photo_label_detection",
                "realized_photo_ID_overlap_fraction",
                "distinct_source_photo_ID_draw_sets",
            )
        }
        data.append({
            "inat_taxon_id":tid,"species":sp,"pilot_wave":wave,
            "source_photo_ids_selected":len(rows),"classifiable_ROI":len(usable),
            "classifiable_fraction":round(len(usable)/len(rows),8),
            "n_classifiable_original_observations":len({x["observation"] for x in usable}),
            "n_classifiable_original_observers":len({x["observer"] for x in usable}),
            "n_classifiable_5degree_cells":len({x["cell"] for x in usable}),
            "photo_label_counts":json.dumps(dict(colour_counts),sort_keys=True),
            "per_label_independent_observers":json.dumps({k:len(v) for k,v in watchers.items()},sort_keys=True),
            "per_label_independent_observations":json.dumps({k:len(v) for k,v in obs_by.items()},sort_keys=True),
            "per_label_5degree_cells":json.dumps({k:len(v) for k,v in cells.items()},sort_keys=True),
            "per_label_observer_folds":json.dumps({k:sorted(v) for k,v in folds.items()},sort_keys=True),
            "within_one_ROI_white_and_yellow_both_ge10":comp_total,
            "within_one_ROI_white_and_yellow_fraction":round(comp_total/len(usable),8) if usable else 0.0,
            "per_photo_label_white_yellow_mixture":json.dumps(mixed,sort_keys=True),
            "independent_observers_cannot_rule_out_single_multicolour_flower":True,
            "genetic_or_population_colour_polymorphism_confirmed":False,
            **repeats,
        })
    if sum(x["source_photo_ids_selected"] for x in data)!=len(real) or sum(
       x["classifiable_ROI"] for x in data)!=roi_doc["n_classifiable_original_flower_photos"]:
        raise ValueError("Actual original ROI and source photo totals did not reconcile")
    outdir.mkdir(parents=True,exist_ok=True)
    table=outdir/f"wave{wave}_original_ROI_observer_geocell_mixed_colour_QC.csv"
    with table.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]),lineterminator="\n")
        w.writeheader();w.writerows(data)
    doc={
        "schema":"fcp_real_rgfca_original_photo_technical_observer_geo_organ_QC_v1",
        "status":"TECHNICAL_REAL_ORIGINAL_ROI_ONLY_NOT_POLYMORPHISM_ESTIMATE",
        "pilot_wave":wave,
        "n_original_ROI_photos":len(real),
        "n_classifiable_flower_ROIs":sum(x["classifiable_ROI"] for x in data),
        "n_taxa":len(data),
        "technical_original_source_receipt_sha256":sha256(source_doc_path),
        "technical_original_roi_receipt_sha256":sha256(roi_doc_path),
        "original_source_date_bulk_archive_membership_verified":False,
        "main_1000_RGFCA_P_poly_estimated":False,
        "post_pilot_structural_colour_review_required":True,
        "source_pilot_random_photo_draws_per_taxon":200,
        "pilot_resampling_fold_depths":[20,32,40,44,48],
        "resampling_does_NOT_create_independent_natural_populations":True,
        "actual_pairwise_source_photo_ID_overlap_by_draw_depth_included":True,
        "no_genetic_FCP_species_asserted":True,
        "CSV_sha256":sha256(table),
        "no_old_colour_measurements_reused":True,
    }
    (outdir/f"wave{wave}_original_ROI_observer_geocell_mixed_colour_QC.json").write_text(
        json.dumps(doc,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    return doc


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pilot-root",required=True,type=Path)
    p.add_argument("--original-rgfca-roster",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    x=p.parse_args(argv)
    print("REAL_RGFCA_ORIGINAL_PHOTO_OBSERVER_MIXED_ROI_QC",
          json.dumps(summarize(x.pilot_root,x.original_rgfca_roster,x.outdir),sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
