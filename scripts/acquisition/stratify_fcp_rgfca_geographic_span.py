#!/usr/bin/env python3
"""Classify OBSERVED geographic breadth before opening FCP/RGFCA flower pixels.

Uses the exact dated 2026-09-27 iNaturalist four-table RGFCA 1,000-taxon
source photograph-ID manifest (3,430,502 records). Keeps the original
outcome-blind 1,000-taxon cohort and 200x250 balanced schedule unchanged.

Broad geography here is DOCUMENTED GEOGRAPHIC PHOTOGRAPHIC SAMPLING EXTENT;
it is NOT verified native range, ecological niche breadth, or genetic
flower-colour polymorphism. Independent observations and independent observers
support both far-apart cells to avoid calling a single mislocated photo a
widespread species. Source counts are photo opportunities, not flower pixels.

The result is a pre-colour, taxon-only priority and matching ledger.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
from collections import Counter, defaultdict
from pathlib import Path

from scripts.acquisition.extract_rgfca_1000_original_monthly_photo_ids import (
    ARCHIVE_SHA, COUNT_SHA, CSV_FIELDS, RGFCA_ENROLLED_PHOTOS, ROSTER_SHA,
    SNAPSHOT,
)
from fcp_pipeline.photo_polymorphism_screen import sha256

SCHEMA="fcp_rgfca_original_1000_precolour_geographic_span_v1"
STATUS="SOURCE_ONLY_GEOGRAPHIC_PRIORITY_NO_FLOWER_PIXELS"
RECEIPT_SCHEMA="fcp_rgfca_targeted_original_20260927_photo_ID_census_v1"
RECEIPT_STATUS="COMPLETE_TARGETED_1000_MONTHLY_SOURCE_IDS_COUNT_MATCHED_NO_PHOTO_PIXELS"
ROSTER_FIELDS=("inat_taxon_id","species","source_photo_opportunities",
               "opportunity_stratum","enrolment_hash")
OUTPUT_FIELDS=(
    "inat_taxon_id","species","opportunity_stratum","source_photo_opportunities",
    "n_original_photo_metadata","n_distinct_observations","n_distinct_observers",
    "n_usable_geographic_observations","n_5degree_cells",
    "n_geographically_supported_cells","max_supported_cluster_separation_km",
    "geographic_radius_95_km","n_unknown_coordinate_accuracy",
    "n_rejected_coordinate_accuracy","geographic_sampling_category",
    "prioritize_broad_before_comparison","broad_vs_narrow_comparison_eligible",
    "is_native_range_confirmed","original_image_pixels_opened",
)
EARTH_KM=6371.0088
URL=re.compile(r"^https://inaturalist-open-data\.s3\.amazonaws\.com/photos/(\d+)/large\.(jpg|jpeg|png)$",re.I)
SEED="FCP_PRECOLOUR_GEOGRAPHIC_SPAN_20261016_V1"
BROAD_KM=1000.0
NARROW_P95_KM=200.0
MAX_POSITIONAL_ACCURACY_M=50000.0
MIN_GLOBAL_OBSERVATIONS=12
MIN_GLOBAL_OBSERVERS=6
MIN_CELL_OBSERVATIONS=3
MIN_CELL_OBSERVERS=2


def _hash(*parts):
    return hashlib.sha256("|".join([SEED,*map(str,parts)]).encode()).hexdigest()


def _xyz(lat,lon):
    p,l=math.radians(lat),math.radians(lon)
    return math.cos(p)*math.cos(l),math.cos(p)*math.sin(l),math.sin(p)


def _angle_km(a,b):
    dot=max(-1.0,min(1.0,sum(x*y for x,y in zip(a,b))))
    return EARTH_KM*math.acos(dot)


def _normalize(vec):
    s=math.sqrt(sum(x*x for x in vec))
    return tuple(x/s for x in vec) if s else (1.,0.,0.)


def _p95(values):
    if not values:return None
    values.sort()
    pos=.95*(len(values)-1)
    lo=int(pos);hi=min(lo+1,len(values)-1)
    return values[lo]+(values[hi]-values[lo])*(pos-lo)


def _cell(lat,lon):
    # -180 and +180 denote the same antimeridian location.
    return (min(35,int(math.floor((lat+90)/5))),
            int(math.floor(((lon+180)%360)/5)))


class SpeciesGeo:
    def __init__(self, tid, name, tier, total):
        self.tid,self.name,self.tier,self.total=tid,name,tier,total
        self.n_photo=0
        self.observations=set()
        self.observers=set()
        self.coordinates=[]
        self.cells=defaultdict(lambda:{"observer":set(),"obs":0,"xyz":[0.,0.,0.]})
        self.unknown_accuracy=0
        self.rejected_accuracy=0

    def add(self,row):
        self.n_photo+=1
        obs,watcher=row["observation_id"],row["observer_id"]
        if not obs or not watcher or not obs.strip() or not watcher.strip():
            raise ValueError("Immutable monthly source lacks observation/observer identifiers")
        if obs in self.observations:
            return # same observation photographed more than once: ONE locality
        self.observations.add(obs)
        self.observers.add(watcher)
        try:
            lat,lon=float(row["latitude"]),float(row["longitude"])
        except (ValueError,TypeError):
            self.rejected_accuracy+=1
            return
        if not math.isfinite(lat) or not math.isfinite(lon) or abs(lat)>90 or abs(lon)>180:
            self.rejected_accuracy+=1
            return
        raw=str(row.get("positional_accuracy_m","")).strip()
        if not raw:
            self.unknown_accuracy+=1
        else:
            try:
                accuracy=float(raw)
                if not math.isfinite(accuracy) or accuracy<0 or accuracy>MAX_POSITIONAL_ACCURACY_M:
                    self.rejected_accuracy+=1
                    return
            except ValueError:
                self.rejected_accuracy+=1
                return
        xyz=_xyz(lat,lon)
        self.coordinates.append(xyz)
        group=self.cells[_cell(lat,lon)]
        group["obs"]+=1
        group["observer"].add(watcher)
        for k in range(3):
            group["xyz"][k]+=xyz[k]

    def finish(self):
        if self.n_photo!=self.total:
            raise ValueError(f"RGFCA source photo ID count drifted: taxon={self.tid}, got={self.n_photo}, expected={self.total}")
        supported=[
            (key,_normalize(v["xyz"]),v["observer"],v["obs"])
            for key,v in self.cells.items()
            if v["obs"]>=MIN_CELL_OBSERVATIONS
            and len(v["observer"])>=MIN_CELL_OBSERVERS
        ]
        max_sep=0.0
        for i in range(len(supported)):
            _,center,watchers,_=supported[i]
            for j in range(i+1,len(supported)):
                _,next_center,next_watchers,_=supported[j]
                if len(watchers|next_watchers)<4:
                    continue
                max_sep=max(max_sep,_angle_km(center,next_center))
        if self.coordinates:
            sums=[sum(x[k] for x in self.coordinates) for k in range(3)]
            center=_normalize(sums)
            p95=_p95([_angle_km(center,x) for x in self.coordinates])
        else:
            p95=None
        n_obs=len(self.observations)
        n_users=len(self.observers)
        if n_obs<MIN_GLOBAL_OBSERVATIONS or n_users<MIN_GLOBAL_OBSERVERS or not self.coordinates:
            cat="geography_insufficient_independent_records"
        elif max_sep>=BROAD_KM:
            cat="broad_observed_1000km_supported"
        elif p95 is not None and p95<=NARROW_P95_KM:
            cat="narrow_observed_p95_le200km"
        else:
            cat="intermediate_or_unresolved_observed_span"
        return {
            "inat_taxon_id":self.tid,"species":self.name,
            "opportunity_stratum":self.tier,
            "source_photo_opportunities":self.total,
            "n_original_photo_metadata":self.n_photo,
            "n_distinct_observations":n_obs,
            "n_distinct_observers":n_users,
            "n_usable_geographic_observations":len(self.coordinates),
            "n_5degree_cells":len(self.cells),
            "n_geographically_supported_cells":len(supported),
            "max_supported_cluster_separation_km":round(max_sep,4),
            "geographic_radius_95_km":round(p95,4) if p95 is not None else "",
            "n_unknown_coordinate_accuracy":self.unknown_accuracy,
            "n_rejected_coordinate_accuracy":self.rejected_accuracy,
            "geographic_sampling_category":cat,
            "prioritize_broad_before_comparison":cat=="broad_observed_1000km_supported",
            "broad_vs_narrow_comparison_eligible":cat in {
                "broad_observed_1000km_supported","narrow_observed_p95_le200km"},
            "is_native_range_confirmed":False,
            "original_image_pixels_opened":False,
        }


def run(manifest:Path,source_receipt:Path,roster:Path,outdir:Path, *,
        synthetic:bool=False)->dict:
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError("Refuse to overwrite previously frozen geographical sampling priority")
    receipt=json.loads(source_receipt.read_text(encoding="utf-8"))
    if (receipt.get("schema")!=RECEIPT_SCHEMA
        or receipt.get("status")!=RECEIPT_STATUS
        or receipt.get("source_snapshot_date")!=SNAPSHOT
        or receipt.get("selected_all_photo_IDs_csv_sha256")!=sha256(manifest)
        or receipt.get("frozen_1000_source_enrolment_sha256")!=sha256(roster)
        or receipt.get("all_enrolled_species_match_independent_full_source_photo_ID_counts") is not True
        or receipt.get("original_image_pixels_opened") is not False
        or receipt.get("photos_NOT_classified_or_species_FCP_inferred") is not True
        or receipt.get("historical_flower_colour_CSV_imported") is not False
        or (not synthetic and (
            sha256(roster)!=ROSTER_SHA
            or receipt.get("source_archive_sha256")!=ARCHIVE_SHA
            or receipt.get("full_42111_photo_opportunity_counts_sha256")!=COUNT_SHA
            or receipt.get("enrolled_outcome_blind_taxa")!=1000
            or receipt.get("n_selected_original_photo_IDs")!=RGFCA_ENROLLED_PHOTOS))):
        raise ValueError("Missing correct untouched RGFCA dated-source ID/photo opportunity certificate")

    taxa={}
    with roster.open("r",encoding="utf-8",newline="") as source:
        reader=csv.DictReader(source)
        if tuple(reader.fieldnames or ())!=ROSTER_FIELDS:
            raise ValueError("Original RGFCA species-enrolment column contract changed")
        for rec in reader:
            key=rec["inat_taxon_id"]
            if not key.isdecimal() or int(key) in taxa or not rec["source_photo_opportunities"].isdecimal():
                raise ValueError("Duplicate or invalid RGFCA original taxon/source photo count")
            tid=int(key)
            taxa[tid]=SpeciesGeo(tid,rec["species"],rec["opportunity_stratum"],
                                 int(rec["source_photo_opportunities"]))
    if not taxa or (not synthetic and len(taxa)!=1000):
        raise ValueError("RGFCA original 1000-taxon cohort missing")
    observed=set()
    last=0
    count=0
    current=None
    results=[]
    with manifest.open("r",encoding="utf-8",newline="") as file:
        reader=csv.DictReader(file)
        if tuple(reader.fieldnames or ())!=CSV_FIELDS:
            raise ValueError("Original monthly full-taxonomy photo record fields changed")
        for rec in reader:
            key=rec["inat_taxon_id"]
            pid=rec["photo_id"]
            if (not key.isdecimal() or not pid.isdecimal()
                or int(key) not in taxa or rec["source_snapshot_date"]!=SNAPSHOT
                or rec["species"]!=taxa[int(key)].name
                or not (m:=URL.fullmatch(rec["photo_url_large"]))
                or m.group(1)!=pid
                or not rec.get("photo_license") or not rec.get("attribution")):
                raise ValueError("Original monthly source photo identity or copyright is not verified")
            tid=int(key)
            if tid!=current:
                if current is not None:
                    results.append(taxa[current].finish())
                if tid<=last or tid in observed:
                    raise ValueError("Source manifests not strictly sorted and unique by frozen taxon")
                last=tid
                current=tid
                observed.add(tid)
            taxa[tid].add(rec)
            count+=1
    if current is not None:results.append(taxa[current].finish())
    if (set(taxa)!=observed
        or count!=sum(s.total for s in taxa.values())
        or count!=receipt.get("n_selected_original_photo_IDs")):
        raise ValueError("Missing an original RGFCA taxon or original source photo ID")
    counts=Counter(x["geographic_sampling_category"] for x in results)
    # Every one of the original 1000 taxa is reported (not just broad taxa).
    def rank(r):
        category={"broad_observed_1000km_supported":0,
                  "narrow_observed_p95_le200km":1,
                  "intermediate_or_unresolved_observed_span":2,
                  "geography_insufficient_independent_records":3}
        return (category[r["geographic_sampling_category"]],
                -r["max_supported_cluster_separation_km"],
                _hash("priority",r["inat_taxon_id"]))
    ordered=sorted(results,key=rank)
    outdir.mkdir(parents=True,exist_ok=True)
    table=outdir/"rgfca_1000_original_source_observed_geographic_span.csv"
    priority=outdir/"rgfca_1000_precolour_broad_geography_first_priority.csv"
    for path,records in ((table,sorted(results,key=lambda r:r["inat_taxon_id"])),
                         (priority,ordered)):
        with path.open("w",encoding="utf-8",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=OUTPUT_FIELDS,lineterminator="\n")
            writer.writeheader()
            writer.writerows(records)
    summary={
        "schema":SCHEMA,"status":STATUS,
        "original_20260927_targeted_source_receipt_sha256":sha256(source_receipt),
        "original_20260927_targeted_manifest_sha256":sha256(manifest),
        "original_1000_species_roster_sha256":sha256(roster),
        "source_species":len(taxa),"source_photo_opportunities":count,
        "geographic_categories":dict(counts),
        "broad_minimum_cluster_separation_km":BROAD_KM,
        "narrow_95pct_centroid_radius_upper_km":NARROW_P95_KM,
        "minimum_independent_observations":MIN_GLOBAL_OBSERVATIONS,
        "minimum_distinct_observers":MIN_GLOBAL_OBSERVERS,
        "supported_5degree_cell_min_distinct_observations":MIN_CELL_OBSERVATIONS,
        "supported_5degree_cell_min_observers":MIN_CELL_OBSERVERS,
        "maximum_accepted_reported_position_accuracy_m":MAX_POSITIONAL_ACCURACY_M,
        "observed_photographic_distribution_not_verified_native_range":True,
        "existing_1000_taxon_enrolment_preserved":True,
        "existing_200x250_repetition_schedule_unchanged":True,
        "selection_before_any_original_photo_pixel":True,
        "geographic_priority_must_not_be_treated_as_random_sampling":True,
        "source_photo_selections_already_sealed_not_rewritten":True,
        "precolour_geographic_span_csv_sha256":sha256(table),
        "precolour_geographic_priority_csv_sha256":sha256(priority),
        "n_original_new_pixels_measured":0,
        "n_FCP_positive_species_classified":0,
        "climate_soil_W_or_H2_estimated":False,
    }
    report=outdir/"rgfca_1000_precolour_geographic_span_receipt.json"
    report.write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("FCP_RGFCA_ORIGINAL_GEO_SPAN_PRECOLOUR",json.dumps({
        "taxa":len(taxa),"source_photo_records":count,
        "class_counts":dict(counts),"no_pixels":True},sort_keys=True))
    return summary


def main(argv=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original-dated-source-photo-ids",type=Path,required=True)
    p.add_argument("--original-dated-source-receipt",type=Path,required=True)
    p.add_argument("--original-rgfca-roster",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args(argv)
    run(a.original_dated_source_photo_ids,a.original_dated_source_receipt,
        a.original_rgfca_roster,a.output)
    return 0

if __name__=="__main__":
    raise SystemExit(main())
