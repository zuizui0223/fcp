#!/usr/bin/env python3
"""Reanalyse real 12-taxon, 1034-photo ROI pilot GEOGRAPHY independently.

Both input waves are existing live-iNaturalist-API technical pilots with
created-before-2026-09-27 filters. They do NOT establish membership in the
immutable monthly iNaturalist 2026-09-27 bulk photo-ID snapshot. The analysis
never opens new photographs or past flower-colour measurement tables.

Photo clusters here describe geographic spread of the PHOTOGRAPHED SAMPLE,
not wild/native species range. Original positional accuracy and photo
observation dates were unavailable in these pilots: do not infer either.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

ROSTER_SHA="88aba565e1bdeef48ad3dc5eec0121c52448fda434aca849e542bb615f9d547c"
SOURCE_STATUS="SELECTED_ORIGINAL_PROVIDER_API_PHOTO_IDS_NO_IMAGE_PIXELS_YET"
SOURCE_KIND="iNaturalist_API_created_before_20260927_live_snapshot_NOT_dated_monthly_archive"
CLASSIFIED="classified_four_state_morph"
MORPHS={"white","yellow_orange","red_pink","blue_purple"}
EARTH_KM=6371.0088
WAVES=(
    (1,"rgfca_six_source_photo_candidates_before_colour.csv",
       "rgfca_six_source_precolour_receipt.json",
       "rgfca_original_new_photo_six_taxa_live_API_technical.csv",
       "rgfca_original_new_photo_six_taxa_live_API_technical_receipt.json",
       "fcp_rgfca_1000_enrolled_real_original_API_six_taxon_pilot_v1",
       "fcp_rgfca_new_original_ROI_v4_live_API_technical_pilot_v1",
       "REAL_NEW_ORIGINAL_PIXELS_TECHNICAL_PILOT_NOT_RGFCA_PPOLY",552,202),
    (2,"rgfca_six_wave2_source_photo_candidates_before_colour.csv",
       "rgfca_six_wave2_source_precolour_receipt.json",
       "rgfca_original_new_photo_six_taxa_wave2_live_API_technical.csv",
       "rgfca_original_new_photo_six_taxa_wave2_live_API_technical_receipt.json",
       "fcp_rgfca_1000_enrolled_real_original_API_disjoint_second_six_taxon_pilot_v1",
       "fcp_rgfca_new_original_ROI_v4_live_API_technical_second_disjoint_pilot_v1",
       "REAL_NEW_ORIGINAL_PIXELS_DISJOINT_WAVE2_TECHNICAL_NOT_RGFCA_PPOLY",482,203),
)
FIELDS=[
 "inat_taxon_id","species","source_opportunity_tier",
 "technical_original_photo_ids","technical_classifiable_photo_ROIs",
 "technical_classifiability_fraction","distinct_original_observations",
 "distinct_independent_observers","photographed_5degree_cells",
 "supported_photographed_5degree_cells",
 "max_supported_sample_geocell_centroid_km",
 "sample_spherical_centroid_p95_radius_km",
 "observed_technical_photographic_spread_category",
 "photo_labels_white","photo_labels_yellow_orange",
 "photo_labels_red_pink","photo_labels_blue_purple",
 "classifiable_same_ROI_white_yellow_both_ge10pct",
 "n_classifiable_photos_in_geographically_supported_cells",
 "original_position_uncertainty_available","true_native_range_verified",
 "monthly_20260927_photo_ID_membership_verified",
]


def sha(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for piece in iter(lambda:f.read(1<<20),b""):
            h.update(piece)
    return h.hexdigest()


def read(path:Path)->list[dict]:
    with path.open(newline="",encoding="utf-8") as file:
        return list(csv.DictReader(file))


def xyz(lat,lon):
    phi,lam=math.radians(lat),math.radians(lon)
    return (math.cos(phi)*math.cos(lam),math.cos(phi)*math.sin(lam),math.sin(phi))


def normalize(v):
    mag=math.sqrt(sum(x*x for x in v))
    if mag<1e-10:
        raise ValueError("Ambiguous globally distributed sample geographic centroid")
    return tuple(x/mag for x in v)


def great_circle(a,b):
    dot=max(-1.,min(1.,sum(x*y for x,y in zip(a,b))))
    return EARTH_KM*math.acos(dot)


def p95(values):
    y=sorted(values);pos=.95*(len(y)-1)
    i=int(pos);f=pos-i
    return y[i]+f*(y[min(len(y)-1,i+1)]-y[i])


def validate_wave(root:Path,signature:tuple)->list[dict]:
    (num,source_name,source_doc_name,measured_name,roi_doc_name,
     source_schema,roi_schema,roi_status,n_total,n_classifiable)=signature
    meta=root/"source_precolour"/source_name
    mreceipt_path=root/"source_precolour"/source_doc_name
    roi=root/"original_roi_results"/measured_name
    rreceipt_path=root/"original_roi_results"/roi_doc_name
    for path in (meta,mreceipt_path,roi,rreceipt_path):
        if not path.is_file():
            raise ValueError(f"Missing original frozen wave {num} input: {path.name}")
    m=json.loads(mreceipt_path.read_text(encoding="utf-8"))
    r=json.loads(rreceipt_path.read_text(encoding="utf-8"))
    if (m.get("schema")!=source_schema
        or m.get("status")!=SOURCE_STATUS
        or m.get("selected_metadata_sha256")!=sha(meta)
        or m.get("source_ROSTER_sha256")!=ROSTER_SHA
        or m.get("NOT_THE_IMMUTABLE_20260927_METADATA_CENSUS") is not True
        or m.get("EXACT_20260927_ARCHIVE_PHOTO_ID_MEMBERSHIP_VERIFIED") is not False
        or m.get("photo_selection_used_climate_soil_or_prior_photo_colours") is not False
        or m.get("photo_selection_used_literature_phenotype") is not False
        or m.get("image_bytes_opened") is not False
        or m.get("RGFCA_200_realization_photo_polymorphism_Ppoly_computed") is not False
        or r.get("schema")!=roi_schema or r.get("status")!=roi_status
        or r.get("source_kind")!=SOURCE_KIND
        or r.get("new_terminal_photo_csv_sha256")!=sha(roi)
        or r.get("source_image_ID_selection_sha256")!=sha(meta)
        or r.get("source_API_selection_receipt_sha256")!=sha(mreceipt_path)
        or r.get("exact_dated_20260927_source_membership_admitted") is not False
        or r.get("RGFCA_200x250_Ppoly_not_estimated") is not True
        or r.get("full_1000_species_ROI_study_complete") is not False
        or r.get("prior_FCP_measured_colours_loaded") is not False
        or m.get("n_selected_original_photo_IDs")!=n_total
        or r.get("selected_real_original_photo_IDs")!=n_total
        or r.get("original_pixels_ROI_v4_evaluated")!=n_total
        or r.get("original_pixels_downloaded")!=n_total
        or r.get("n_classifiable_original_flower_photos")!=n_classifiable):
        raise ValueError("Real original-photo ROI/source/model receipt mismatch or monthly-source forgery")
    selected=read(meta);pixels=read(roi)
    if len(selected)!=len(pixels)!=n_total:
        raise ValueError("Wrong number of selected ORIGINAL photo ID records")
    source_by_id={x["photo_id"]:x for x in selected}
    pixel_by_id={x["photo_id"]:x for x in pixels}
    if len(source_by_id)!=len(pixel_by_id)!=n_total or set(source_by_id)!=set(pixel_by_id):
        raise ValueError("Duplicate/replaced original source photo ID")
    for pid,x in pixel_by_id.items():
        y=source_by_id[pid]
        for key in (
            "inat_taxon_id","species","source_opportunity_tier","photo_id",
            "observation_id","observer_id","latitude","longitude",
            "photo_url_large","photo_license","screen_fold","sample_hash"):
            if x[key]!=y[key]:
                raise ValueError("Technical colour outcome/source original photograph identity switched")
        classified=x["measurement_status"]==CLASSIFIED
        if classified!=(x["morph"] in MORPHS):
            raise ValueError("Real frozen source ROI four-state terminal status invalid")
        try:
            lat,lon=float(y["latitude"]),float(y["longitude"])
        except ValueError as exc:
            raise ValueError("Invalid original photo geography") from exc
        if not(math.isfinite(lat) and math.isfinite(lon) and abs(lat)<=90 and abs(lon)<=180):
            raise ValueError("Invalid original source photo location")
    if sum(x["measurement_status"]==CLASSIFIED for x in pixels)!=n_classifiable:
        raise ValueError("Four-state original source photo classifiable count drift")
    if len({x["inat_taxon_id"] for x in pixels})!=6:
        raise ValueError("Technical six-species cohort was changed")
    return pixels


def summarize(rows:list[dict],outdir:Path)->dict:
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError("Do not overwrite original photo geography/colour technical audit")
    unique_ids={x["photo_id"] for x in rows}
    if len(rows)!=1034 or len(unique_ids)!=1034 or len({x["inat_taxon_id"] for x in rows})!=12:
        raise ValueError("Technical original photo rows overlap or cohort denominator drift")
    by=defaultdict(list)
    for r in rows:
        by[(r["inat_taxon_id"],r["species"],r["source_opportunity_tier"])].append(r)
    output=[]
    for (taxon,species,tier),pictures in sorted(by.items()):
        # Multiple photographed views of one observation are one geographic
        # record, not independent replicated plants or independent localities.
        cells=defaultdict(lambda:{"observation":set(),"observer":set(),"vectors":[],"classified":0})
        photos_by_obs={}
        colours=Counter()
        mixed=0
        for x in pictures:
            lat,lon=float(x["latitude"]),float(x["longitude"])
            cell=(min(35,math.floor((lat+90)/5)),math.floor(((lon+180)%360)/5))
            obs=x["observation_id"]
            if obs in photos_by_obs and photos_by_obs[obs]!=(lat,lon,x["observer_id"]):
                raise ValueError("One source observation has inconsistent location or observer")
            if obs not in photos_by_obs:
                photos_by_obs[obs]=(lat,lon,x["observer_id"])
                cells[cell]["observation"].add(obs)
                cells[cell]["observer"].add(x["observer_id"])
                cells[cell]["vectors"].append(xyz(lat,lon))
            if x["measurement_status"]==CLASSIFIED:
                colours[x["morph"]]+=1
                cells[cell]["classified"]+=1
                if float(x["flower_fraction_white"])>=.1 and float(x["flower_fraction_yellow"])>=.1:
                    mixed+=1
        supports={key:val for key,val in cells.items()
                  if len(val["observation"])>=3 and len(val["observer"])>=2}
        centroids={}
        for cell,rec in supports.items():
            pts=rec["vectors"]
            centroids[cell]=normalize(tuple(sum(p[i] for p in pts) for i in range(3)))
        max_dist=0.
        for i,c in enumerate(centroids):
            for d in list(centroids)[i+1:]:
                if len(supports[c]["observer"]|supports[d]["observer"])<4:continue
                max_dist=max(max_dist,great_circle(centroids[c],centroids[d]))
        positions=[xyz(lat,lon) for lat,lon,_ in photos_by_obs.values()]
        center=normalize(tuple(sum(p[i] for p in positions) for i in range(3)))
        radius95=p95([great_circle(center,p) for p in positions])
        obs_total=len(photos_by_obs)
        observers_total=len({x[2] for x in photos_by_obs.values()})
        if obs_total<12 or observers_total<6:
            cat="source_sample_insufficient_independent_observations"
        elif max_dist>=1000:
            cat="sample_supported_observed_geography_ge1000km_ACCURACY_UNKNOWN"
        elif radius95<=200:
            cat="sample_narrow_p95_le200km_ACCURACY_UNKNOWN"
        else:
            cat="sample_intermediate_or_unresolved_ACCURACY_UNKNOWN"
        cnt=sum(colours.values())
        output.append({
            "inat_taxon_id":taxon,"species":species,"source_opportunity_tier":tier,
            "technical_original_photo_ids":len(pictures),
            "technical_classifiable_photo_ROIs":cnt,
            "technical_classifiability_fraction":round(cnt/len(pictures),6),
            "distinct_original_observations":obs_total,
            "distinct_independent_observers":observers_total,
            "photographed_5degree_cells":len(cells),
            "supported_photographed_5degree_cells":len(supports),
            "max_supported_sample_geocell_centroid_km":round(max_dist,3),
            "sample_spherical_centroid_p95_radius_km":round(radius95,3),
            "observed_technical_photographic_spread_category":cat,
            "photo_labels_white":colours["white"],
            "photo_labels_yellow_orange":colours["yellow_orange"],
            "photo_labels_red_pink":colours["red_pink"],
            "photo_labels_blue_purple":colours["blue_purple"],
            "classifiable_same_ROI_white_yellow_both_ge10pct":mixed,
            "n_classifiable_photos_in_geographically_supported_cells":
                sum(q["classified"] for q in supports.values()),
            "original_position_uncertainty_available":False,
            "true_native_range_verified":False,
            "monthly_20260927_photo_ID_membership_verified":False,
        })
    total_class=sum(int(x["technical_classifiable_photo_ROIs"]) for x in output)
    if total_class!=405:
        raise ValueError("Real 1034 original photo classifiable totals inconsistent")
    classification=Counter(r["observed_technical_photographic_spread_category"] for r in output)
    outdir.mkdir(parents=True,exist_ok=True)
    table=outdir/"rgfca_12_REAL_original_photo_geography_vs_colour_technical.csv"
    with table.open("x",newline="",encoding="utf-8") as handle:
        w=csv.DictWriter(handle,fieldnames=FIELDS,lineterminator="\n")
        w.writeheader();w.writerows(output)
    summary={
        "schema":"fcp_rgfca_12_real_live_api_original_photo_geographic_sampling_QC_v1",
        "status":"VERIFIED_REAL_PIXELS_TECHNICAL_GEOGRAPHIC_SPAN_NOT_MONTHLY_RGFCA_FCP",
        "original_roi_pilot_actions_runs":[38832829706,38840579008],
        "original_roi_pilot_artifact_ids":[12018829689,12022221110],
        "data_source":"live_iNaturalist_observations_API_created_before_20260927_not_exact_monthly_archive",
        "original_rgfca_species_enrolment_sha256":ROSTER_SHA,
        "n_distinct_real_original_photo_IDs":len(unique_ids),
        "n_original_photo_classes_classifiable":total_class,
        "n_technically_evaluated_taxa":len(output),
        "n_sample_photo_supported_geo_ge1000km_unverified_position_accuracy":
            classification["sample_supported_observed_geography_ge1000km_ACCURACY_UNKNOWN"],
        "sample_geographic_categories":dict(classification),
        "not_original_20260927_monthly_source_IDs_admitted":True,
        "original_positional_accuracy_not_provided":True,
        "original_photo_observation_dates_not_provided":True,
        "no_new_photo_images_downloaded":True,
        "not_verified_wild_or_native_species_geographic_ranges":True,
        "no_genetic_polymorphism_claims":True,
        "no_RGFCA_main_1000_species_Ppoly_W_H2_or_environment_estimated":True,
        "same_flower_ROI_mixed_white_yellow_means_Ppoly_NOT_established":True,
        "original_real_photo_table_sha256":sha(table),
    }
    receipt=outdir/"rgfca_12_real_photo_geographic_technical_receipt.json"
    receipt.write_text(json.dumps(summary,indent=2,ensure_ascii=False,sort_keys=True)+"\n",encoding="utf-8")
    print("FCP_RGFCA_REAL12_GEO_PHOTO_QUALITY_QC",json.dumps({
        "n_real_original_photos":len(unique_ids),
        "n_classifiable":total_class,
        "n_source_species":len(output),
        "geographic_categories":dict(classification)},sort_keys=True),flush=True)
    return summary


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--wave1-root",required=True,type=Path)
    p.add_argument("--wave2-root",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    x=p.parse_args(argv)
    photos=validate_wave(x.wave1_root,WAVES[0])+validate_wave(x.wave2_root,WAVES[1])
    summarize(photos,x.outdir)


if __name__=="__main__":
    main()
