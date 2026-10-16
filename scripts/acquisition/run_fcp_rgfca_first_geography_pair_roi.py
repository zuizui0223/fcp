#!/usr/bin/env python3
"""Measure the FIRST source-supported broad/narrow taxon PAIR, not an arbitrary photo batch.

Only runs after independent complete 2026-09-27 source-ID, original image
selection and 128 worker packet hashes have been verified. Selects the first
metadata-eligible pre-colour broad vs narrow taxon pair, 4..16 original images
from each independent observer fold in EACH species (16..64 images total).
Source geography/identity never enters the exact frozen ROI-v4 evaluator.

The result is the first authentic MONTHLY-SOURCE pixel technical work unit,
not 1,000-taxon coverage, genetic polymorphism, a geographic effect estimate,
or a sample-independent randomized FCP inference.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from datetime import date
from collections import defaultdict, Counter
from pathlib import Path

import pandas as pd

from fcp_pipeline.photo_polymorphism_screen import sha256
from scripts.acquisition.build_fcp_uncapped_blind_batches import new_id
from scripts.acquisition.plan_fcp_photo_first_screen import SELECTED_FIELDS
from scripts.acquisition.prepare_first_rgfca_original_monthly_roi_unit import prepare
from scripts.acquisition.run_fcp_fresh_image_partition import run_batch
from scripts.acquisition.freeze_fcp_photo_first_blind_workers import WORKER,ACQUISITION
from scripts.acquisition.stratify_fcp_rgfca_geographic_span import _xyz, _angle_km

PHOTO_RE=re.compile(
    r"^https://inaturalist-open-data\.s3\.amazonaws\.com/photos/([0-9]+)/large\.(jpg|jpeg|png)$",
    re.I,
)
SOURCE_DATE="2026-09-27"
GEO_STATUS="SOURCE_ONLY_GEOGRAPHIC_PRIORITY_NO_FLOWER_PIXELS"
PHOTO_SCHEMA="fcp_original_rgfca_geographically_matched_real_roi_first_unit_v1"

SELECTION_JOIN=(
    "measurement_id","inat_taxon_id","species","photo_id",
    "observation_id","observer_id","latitude","longitude","screen_fold",
)
OUTPUT_COLUMNS=(*SELECTION_JOIN,"positional_accuracy_m","observed_on",
                "geographic_comparison_group","matched_pair_id")
PER_FOLD_MIN=4
PER_FOLD_TARGET=16


def read_doc(p:Path)->dict:
    v=json.loads(p.read_text(encoding="utf-8"))
    if not isinstance(v,dict):
        raise ValueError("Frozen original dated source receipt must be JSON object")
    return v


def read_csv(path:Path)->list[dict]:
    with path.open("r",newline="",encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def geographical_cell(row:dict)->tuple[int,int]:
    try:
        lat,lon=float(row["latitude"]),float(row["longitude"])
    except (TypeError,ValueError,KeyError) as exc:
        raise ValueError("Frozen original photo geographic coordinate invalid") from exc
    if (not math.isfinite(lat) or not math.isfinite(lon)
        or abs(lat)>90 or abs(lon)>180):
        raise ValueError("Original source photo location outside valid sphere")
    # Both +180 and -180 map to one longitude grid rather than two.
    return (min(35,int(math.floor((lat+90)/5))),
            int(math.floor(((lon+180)%360)/5)))


def geo_balanced_first_photos(pool:list[dict], target:int)->list[dict]:
    """First-stage original pixels geographically balanced WITHOUT flower labels.

    A pre-existing observer-disjoint source pool is partitioned into 5-degree
    sampling cells, then a deterministic round robin spreads the first
    few real pixels across localities instead of re-sorting all by photo ID
    hash alone. The original 64+64 source pool is NEVER reselected.
    """
    if not 1<=target<=16 or len(pool)<target:
        raise ValueError("Geographic first-pixel resource depth cannot be met")
    cells=defaultdict(list)
    for row in pool:
        cells[geographical_cell(row)].append(row)
    for rows in cells.values():
        rows.sort(key=lambda x:(x["selection_hash"],int(x["photo_id"])))
    selected=[]
    ordered=sorted(cells)
    while len(selected)<target:
        progressed=False
        for cell in ordered:
            if cells[cell]:
                selected.append(cells[cell].pop(0))
                progressed=True
                if len(selected)==target:
                    break
        if not progressed:
            raise ValueError("Source geographic-balanced first-pixel allocation incomplete")
    if len({r["photo_id"] for r in selected})!=len(selected):
        raise ValueError("Geographic first-pixel source ID reused")
    return selected


def supported_broad_photo_spread(selected_rows:list[dict])->dict:
    """Require selected pixels, not just the parent species, to span two real sites.

    A far-away single photo cannot satisfy the site support requirement. Each
    geographically separated 5-degree cell needs >=2 distinct source
    observations AND >=2 distinct observers. Position accuracy >50 km never
    supports the geographic criterion. This inspects source metadata only.
    """
    cells=defaultdict(lambda:{"observations":set(),"observers":set(),"xyz":[]})
    for row in selected_rows:
        raw=str(row.get("positional_accuracy_m","")).strip()
        # The subset accuracy is recovered from the independently SHA-frozen
        # monthly source photo-ID manifest. Do not allow missing or >50km
        # original positions to create a false wide-range first ROI batch.
        if not raw:
            continue
        try:
            accuracy=float(raw)
        except (TypeError,ValueError):
            continue
        if not math.isfinite(accuracy) or not 0<=accuracy<=50000:
            continue
        cell=geographical_cell(row)
        obs,observer=str(row["observation_id"]),str(row["observer_id"])
        g=cells[cell]
        if obs in g["observations"]:
            continue
        g["observations"].add(obs)
        g["observers"].add(observer)
        g["xyz"].append(_xyz(float(row["latitude"]),float(row["longitude"])))
    supported=[]
    for cell,source in cells.items():
        if len(source["observations"])<2 or len(source["observers"])<2:
            continue
        xyz=source["xyz"]
        mean=[sum(x[k] for x in xyz) for k in range(3)]
        length=math.sqrt(sum(x*x for x in mean))
        if length<1e-10:
            continue
        supported.append((cell,tuple(x/length for x in mean),source["observers"]))
    max_km=0.
    n_independent_distant_pairs=0
    for a in range(len(supported)):
        for b in range(a+1,len(supported)):
            # The original species-level broad-range certificate requires
            # at least FOUR distinct reporters across the two regions.
            # Never manufacture broad replication with the same two people
            # photographing both distant places.
            if len(supported[a][2]|supported[b][2])<4:
                continue
            d=_angle_km(supported[a][1],supported[b][1])
            max_km=max(max_km,d)
            n_independent_distant_pairs+=int(d>=1000)
    return {"n_supported_selected_cells":len(supported),
            "n_distant_cell_pairs_with_four_distinct_source_observers":
                n_independent_distant_pairs,
            "max_supported_selected_cell_centroid_distance_km":round(max_km,3),
            "geographic_accuracy_of_selected_subset_independently_verified":True,
            "selected_positions_from_exact_original_monthly_photo_IDs":True,
            "source_geo_classification_provenance_verified":True}


def reserve_two_supported_distant_sites(
    folds:dict[int,list[dict]], per_fold:int,
)->tuple[dict[int,list[dict]],dict] | None:
    """Start with independently supported distant sites BEFORE normal geo fill.

    The existing colour-blind preselected photos are never replaced. Simple
    round-robin over dozens of 5-degree cells can leave only ONE photo at each
    distant locality, failing a legitimate broad-vs-narrow comparison. Reserve
    two different original observations AND observers in two distant cells,
    while preserving each pre-existing observer-disjoint fold's photo budget.
    The chosen sites depend ONLY on original source geolocation and IDs.
    """
    if set(folds)!={0,1} or not 1<=per_fold<=16:
        raise ValueError("Geographically supported original source folds invalid")
    if any(len(folds[fold])<per_fold for fold in (0,1)):
        raise ValueError("Frozen source lacks original photos for both observer folds")
    candidates=defaultdict(list)
    for fold in (0,1):
        for row in folds[fold]:
            raw=str(row.get("positional_accuracy_m","")).strip()
            try:
                accuracy=float(raw)
            except (TypeError,ValueError):
                continue
            if not math.isfinite(accuracy) or not 0<=accuracy<=50000:
                continue
            candidates[geographical_cell(row)].append((fold,row))
    paired={}
    for cell,entries in candidates.items():
        entries=sorted(entries,key=lambda fr:(
            fr[1]["selection_hash"],int(fr[1]["photo_id"]),fr[0]))
        # Deterministic first admissible pair, never photo colour dependent.
        for i,(a,first) in enumerate(entries):
            for b,second in entries[i+1:]:
                if (first["observation_id"]!=second["observation_id"] and
                    first["observer_id"]!=second["observer_id"]):
                    vectors=[
                        _xyz(float(r["latitude"]),float(r["longitude"]))
                        for r in (first,second)
                    ]
                    xyz=[sum(v[k] for v in vectors) for k in range(3)]
                    size=math.sqrt(sum(v*v for v in xyz))
                    if size<1e-9:
                        continue
                    paired[cell]=(
                        ((a,first),(b,second)),
                        tuple(v/size for v in xyz)
                    )
                    break
            if cell in paired:
                break
    if len(paired)<2:
        return None
    pairs=[]
    for left in sorted(paired):
        for right in sorted(paired):
            if right<=left:continue
            p1,c1=paired[left]
            p2,c2=paired[right]
            km=_angle_km(c1,c2)
            if km<1000:continue
            reserve=[*p1,*p2]
            counts=Counter(fold for fold,_ in reserve)
            if any(counts[f]>per_fold for f in (0,1)):
                continue
            pairs.append((-km,left,right,reserve))
    for _,left,right,reserve in sorted(pairs):
        by_fold={0:[],1:[]}
        used=set()
        for fold,row in reserve:
            if row["photo_id"] in used:
                raise ValueError("Original photo ID reserved twice for sampled locality")
            used.add(row["photo_id"])
            by_fold[fold].append(row)
        for fold in (0,1):
            remaining=per_fold-len(by_fold[fold])
            if remaining>0:
                pool=[r for r in folds[fold] if r["photo_id"] not in used]
                if len(pool)<remaining:
                    continue
                by_fold[fold].extend(geo_balanced_first_photos(pool,remaining))
        if any(len(by_fold[f])!=per_fold for f in (0,1)):
            continue
        combined=by_fold[0]+by_fold[1]
        if len({r["photo_id"] for r in combined})!=2*per_fold:
            raise ValueError("Photo ID collision in geographically balanced pre-colour source")
        supported=supported_broad_photo_spread(combined)
        if (supported["n_supported_selected_cells"]>=2 and
            supported["max_supported_selected_cell_centroid_distance_km"]>=1000):
            supported["distant_support_regions_reserved_before_regular_sampling"]=True
            supported["reserved_geographic_cells"]=[list(left),list(right)]
            supported["observer_fold_source_photo_count_each"]=per_fold
            return by_fold,supported
    return None


def select_unit(source_root:Path,geo_root:Path,outdir:Path,*,
                _synthetic_source_gate:bool=False)->dict:
    """Create the independently source-signed, colour-blind first broad/narrow unit."""
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError("Cannot overwrite a frozen geographically matched ROI unit")
    full=source_root/"rgfca_1000_exact_source_photo_IDs"
    complete_manifest=full/"rgfca_1000_complete_original_source_photo_IDs.csv"
    source_doc=full/"rgfca_1000_targeted_monthly_photo_ID_receipt.json"
    selection=source_root/"rgfca_1000_photo_colour_blind_selection"
    sample_doc=selection/"photo_screen_source_selection_receipt.json"
    sample_csv=selection/"photo_screen_original_metadata_selected.csv"
    firewall=source_root/"rgfca_1000_original_photo_blind_worker_freeze"
    frozen_doc=firewall/"photo_screen_freeze.json"
    geo_doc=geo_root/"rgfca_1000_precolour_geographic_span_receipt.json"
    geo_taxa_csv=geo_root/"rgfca_1000_original_source_observed_geographic_span.csv"
    geo_pair_csv=geo_root/"rgfca_1000_original_source_broad_narrow_depthmatched_pairs.csv"
    for p in (source_doc,complete_manifest,sample_doc,sample_csv,frozen_doc,geo_doc,geo_taxa_csv,
              geo_pair_csv,firewall/"metadata_join_key.csv"):
        if not p.is_file():
            raise ValueError(f"Missing immutable 2026-09-27 pre-ROI source: {p.name}")
    # The full historical original source/128 worker receipts must be checked
    # before selecting or opening a single original image. The synthetic-only
    # override never exists on the public CLI.
    if not _synthetic_source_gate:
        prepare(source_root,outdir/"precolour_full_source_gate")
    full_doc=read_doc(source_doc)
    sample=read_doc(sample_doc)
    freeze=read_doc(frozen_doc)
    geo=read_doc(geo_doc)
    if (full_doc.get("source_snapshot_date")!=SOURCE_DATE
        or full_doc.get("selected_all_photo_IDs_csv_sha256")!=sha256(complete_manifest)
        or sample.get("source_snapshot_date")!=SOURCE_DATE
        or geo.get("status")!=GEO_STATUS
        or geo.get("source_species")!=(4 if _synthetic_source_gate else 1000)
        or geo.get("original_20260927_targeted_source_receipt_sha256")!=sha256(source_doc)
        or geo.get("original_20260927_targeted_manifest_sha256")!=
            full_doc.get("selected_all_photo_IDs_csv_sha256")
        or geo.get("precolour_geographic_span_csv_sha256")!=sha256(geo_taxa_csv)
        or geo.get("geographic_matching_pair_csv_sha256")!=sha256(geo_pair_csv)
        or sample.get("photo_selection_csv_sha256")!=sha256(sample_csv)
        or freeze.get("pre_colour_selection_csv_sha256")!=sha256(sample_csv)
        or freeze.get("source_photo_screen_selection_receipt_sha256")!=sha256(sample_doc)
        or freeze.get("image_worker_sees_only_three_fields") is not True
        or freeze.get("no_old_photo_colour_rows_loaded") is not True
        or geo.get("selection_before_any_original_photo_pixel") is not True
        or geo.get("n_original_new_pixels_measured")!=0
        or full_doc.get("original_image_pixels_opened") is not False
        or full_doc.get("historical_flower_colour_CSV_imported") is not False
        or sample.get("selected_sample_photo_pixels_measured") is not False
        or freeze.get("original_image_pixels_opened") is not False):
        raise ValueError("Geographic pair is not signed by the complete frozen ORIGINAL monthly photo ID universe")

    groups={}
    for row in read_csv(geo_taxa_csv):
        tid=row["inat_taxon_id"]
        if tid in groups:
            raise ValueError("Duplicated species geography source identity")
        groups[tid]=row
    pairs=read_csv(geo_pair_csv)
    if not pairs:
        raise ValueError("No source-supported broad/narrow comparison pairs; keep 1000 taxa and report no matched control")

    with sample_csv.open("r",newline="",encoding="utf-8") as f:
        reader=csv.DictReader(f)
        if tuple(reader.fieldnames or ())!=tuple(SELECTED_FIELDS):
            raise ValueError("Frozen original pre-colour photo selection columns changed")
        source_photo_ids=set()
        photo_rows=defaultdict(lambda:{0:[],1:[]})
        photo_count=0
        for row in reader:
            pid=row["photo_id"]
            if (not pid.isdecimal() or pid in source_photo_ids
                or row["source_snapshot_date"]!=SOURCE_DATE
                or not (url:=PHOTO_RE.fullmatch(row["photo_url_large"]))
                or url.group(1)!=pid
                or row["screen_fold"] not in ("0","1")
                or not row["observer_id"].strip()
                or not row["observation_id"].strip()
                or not row["photo_license"]
                or len(row["selection_hash"])!=16):
                raise ValueError("Photo ID, observer fold, source URL or license deviated from original source")
            source_photo_ids.add(pid)
            photo_count+=1
            photo_rows[row["inat_taxon_id"]][int(row["screen_fold"])].append(row)
    if photo_count!=sample.get("selected_photo_IDs"):
        raise ValueError("Original selected source-photo-ID count drifted")
    if set(photo_rows)-set(groups):
        raise ValueError("Original photo selection contains unrostered geographic taxon")

    # The blinded screen strips positional accuracy, but the exact original
    # 2026-09-27 monthly photo ledger retains it. Recover those ORIGINAL bytes
    # for every selected photo ID, verifying observation and coordinates, so
    # a lone poorly georeferenced photo cannot create false first-pixel range.
    selected_by_photo_id={
        row["photo_id"]:row
        for folds in photo_rows.values() for rows in folds.values() for row in rows
    }
    recovered=set()
    with complete_manifest.open("r",newline="",encoding="utf-8") as original_file:
        original=csv.DictReader(original_file)
        needed={"photo_id","observation_id","observer_id","inat_taxon_id",
                "latitude","longitude","positional_accuracy_m","observed_on",
                "source_snapshot_date"}
        if not needed.issubset(original.fieldnames or ()):
            raise ValueError("Original monthly source ID ledger lost georeference evidence")
        for record in original:
            pid=record["photo_id"]
            if pid not in selected_by_photo_id:
                continue
            if pid in recovered:
                raise ValueError("Original monthly source photo ID repeated")
            selected=selected_by_photo_id[pid]
            if (record["source_snapshot_date"]!=SOURCE_DATE or
                any(selected[k]!=record[k] for k in
                    ("inat_taxon_id","observation_id","observer_id")) or
                any(abs(float(selected[k])-float(record[k]))>1e-7
                    for k in ("latitude","longitude"))):
                raise ValueError("Selected original photographed site differs from immutable monthly record")
            selected["positional_accuracy_m"]=record["positional_accuracy_m"]
            # The pre-colour photo-selection CSV intentionally omitted the
            # original observation date. Recover it *only* from the original
            # SHA-sealed monthly photo ID universe, never the live API.
            selected["observed_on"]=record["observed_on"]
            recovered.add(pid)
    if recovered!=source_photo_ids:
        raise ValueError("Selected original photo ID absent from exact dated monthly archive")

    first_pair=None
    chosen=None
    each=None
    selected_geo_support=None
    for pair in pairs:
        b,n=pair["broad_taxon_id"],pair["narrow_taxon_id"]
        if (b not in groups or n not in groups
            or groups[b]["geographic_sampling_category"]!="broad_observed_1000km_supported"
            or groups[n]["geographic_sampling_category"]!="narrow_observed_p95_le200km"
            or groups[b]["opportunity_stratum"]!=groups[n]["opportunity_stratum"]
            or groups[b]["opportunity_stratum"]!=pair["opportunity_stratum"]
            or max(int(groups[b]["source_photo_opportunities"]),
                   int(groups[n]["source_photo_opportunities"])) >
                   4*min(int(groups[b]["source_photo_opportunities"]),
                         int(groups[n]["source_photo_opportunities"]))
            or pair["matched_before_any_original_photo_pixels"]!="True"):
            raise ValueError("Geographic comparison pair is ungrounded or unmatched by source opportunity")
        if not all(len(photo_rows[tid][fold])>=PER_FOLD_MIN
                   for tid in (b,n) for fold in (0,1)):
            continue
        depth=min(PER_FOLD_TARGET,*(len(photo_rows[tid][fold])
                                   for tid in (b,n) for fold in (0,1)))
        if not PER_FOLD_MIN<=depth<=PER_FOLD_TARGET:
            raise ValueError("Invalid immutable original ROI image compute batch bound")
        spread=reserve_two_supported_distant_sites(photo_rows[b],depth)
        if spread is None:
            # Full colour-blind 64+64 photo pool has inadequate independent
            # regional support even if the 1000-species metadata range is wide.
            continue
        broad_folds,supported=spread
        candidate=[]
        for fold in (0,1):
            candidate.extend((entry,"broad") for entry in broad_folds[fold])
            candidate.extend((entry,"narrow") for entry in
                geo_balanced_first_photos(photo_rows[n][fold],depth))
        first_pair=pair
        chosen=candidate
        each=depth
        selected_geo_support=supported
        break
    if first_pair is None:
        raise ValueError(
            "No matched broad/narrow pair has both observer folds AND "
            "two supported photographed regions >=1000 km apart")
    b,n=first_pair["broad_taxon_id"],first_pair["narrow_taxon_id"]
    if not 16<=len(chosen)<=64 or len({r["photo_id"] for r,_ in chosen})!=len(chosen):
        raise ValueError("Broad/narrow original photo unit improperly selected")
    geographic_cells_by_group={
        name:len({geographical_cell(row) for row,label in chosen if label==name})
        for name in ("broad","narrow")
    }

    # Join by stable blinded measurement_id to the *SHA-frozen* photo source
    # biological metadata. All original worker/acquisition files were verified
    # by prepare() above; nothing other than THREE original fields is sent
    # into the historical ROI-v4 source/model.
    wanted={new_id(r["photo_id"]):(r,label) for r,label in chosen}
    source_join={}
    with (firewall/"metadata_join_key.csv").open("r",newline="",encoding="utf-8") as f:
        reader=csv.DictReader(f)
        if not set(SELECTION_JOIN).issubset(reader.fieldnames or []):
            raise ValueError("Original frozen metadata/source biological join absent")
        for row in reader:
            mid=row["measurement_id"]
            if mid in wanted:
                if mid in source_join:
                    raise ValueError("Original source measurement key was duplicated")
                src,_=wanted[mid]
                for col in ("inat_taxon_id","species","photo_id",
                            "observation_id","observer_id","screen_fold"):
                    if row[col]!=src[col]:
                        raise ValueError("Original frozen source/photo observer identity changed")
                source_join[mid]=row
    if set(source_join)!=set(wanted):
        raise ValueError("Selected original photo IDs not found in the source-frozen 128-packet firewall")

    blind=[]
    acquisition=[]
    sidecar=[]
    for entry,label in chosen:
        mid=new_id(entry["photo_id"])
        filename=mid+".jpg"
        blind.append({"measurement_id":mid,"image_filename":filename,
                      "photo_license":entry["photo_license"]})
        acquisition.append({"measurement_id":mid,"image_filename":filename,
                            "photo_url_large":entry["photo_url_large"],
                            "photo_license":entry["photo_license"]})
        sidecar.append({**{k:source_join[mid][k] for k in SELECTION_JOIN},
                        "positional_accuracy_m":entry["positional_accuracy_m"],
                        "observed_on":entry["observed_on"],
                        "geographic_comparison_group":label,
                        "matched_pair_id":first_pair["pair_id"]})
    outdir.mkdir(parents=True,exist_ok=True)
    def save_rows(filename,cols,rows):
        path=outdir/filename
        with path.open("x",newline="",encoding="utf-8") as f:
            writer=csv.DictWriter(f,fieldnames=cols,lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        return sha256(path)
    blinded_sha=save_rows("first_matched_pair_original_blind_worker.csv",WORKER,blind)
    acquisition_sha=save_rows("first_matched_pair_original_photo_source_acquisition.csv",ACQUISITION,acquisition)
    source_sha=save_rows("first_matched_pair_original_source_join_key.csv",OUTPUT_COLUMNS,sidecar)
    summary={
        "schema":PHOTO_SCHEMA,
        "status":"ORIGINAL_MONTHLY_SOURCE_GEOGRAPHIC_PAIR_PREPARED_NO_PIXELS",
        "source_snapshot_date":SOURCE_DATE,
        "n_enrolled_species":geo["source_species"],
        "n_original_source_photo_ids_in_monthly_cohort":full_doc["n_selected_original_photo_IDs"],
        "n_photo_id_opportunities_selected_precolour":photo_count,
        "precolour_source_selection_sha256":sha256(sample_csv),
        "original_source_ID_receipt_sha256":sha256(source_doc),
        "frozen_worker_receipt_sha256":sha256(frozen_doc),
        "geographic_pair_csv_sha256":sha256(geo_pair_csv),
        "geographic_span_receipt_sha256":sha256(geo_doc),
        "matched_pair_id":first_pair["pair_id"],
        "broad_taxon_id":b,"narrow_taxon_id":n,
        "opportunity_stratum":first_pair["opportunity_stratum"],
        "per_taxon_per_observer_fold_photo_depth":each,
        "number_original_photo_IDs_in_first_unit":len(chosen),
        "n_broad_original_photos":2*each,
        "n_narrow_original_photos":2*each,
        "n_source_geographic_cells_selected_by_group":geographic_cells_by_group,
        "original_photo_first_unit_geographic_cell_round_robin":True,
        "broad_original_first_pixel_geographic_support":selected_geo_support,
        "source_location_accuracy_reverified_in_first_ROI_photo_subset":True,
        "original_photo_accuracy_kept_out_of_blind_image_evaluator":True,
        "original_observation_date_recovered_from_dated_monthly_source":True,
        "original_observation_date_kept_out_of_blind_image_evaluator":True,
        "original_blind_worker_sha256":blinded_sha,
        "original_source_acquisition_sha256":acquisition_sha,
        "original_biological_join_sha256":source_sha,
        "n_original_photo_pixels_opened":0,
        "selection_uses_no_flower_colour_or_environment":True,
        "source_image_IDs_from_verified_monthly_20260927":not _synthetic_source_gate,
        "new_colour_W_H2_or_geographic_association_estimated":False,
        "full_RGFCAROI1000_complete":False,
    }
    (outdir/"first_matched_pair_precolour_source_receipt.json").write_text(
        json.dumps(summary,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    return summary


class PairSourceIndex:
    def __init__(self,source_root:Path):
        self.worker=pd.read_csv(source_root/"first_matched_pair_original_blind_worker.csv",
                                dtype=str).fillna("")
        self.acq=pd.read_csv(source_root/"first_matched_pair_original_photo_source_acquisition.csv",
                             dtype=str).fillna("")
        if tuple(self.worker.columns)!=WORKER or tuple(self.acq.columns)!=ACQUISITION:
            raise ValueError("Original ROI photo worker is not 3-field blind")
        if (not 1<=len(self.worker)<=64
            or len(self.acq)!=len(self.worker)
            or len(set(self.worker.measurement_id))!=len(self.worker)
            or len(set(self.acq.measurement_id))!=len(self.acq)):
            raise ValueError("Original matched geographic ROI photo IDs are not unique")
        joined=self.worker.merge(self.acq,on=list(WORKER),how="left",validate="one_to_one")
        if joined.photo_url_large.isna().any():
            raise ValueError("Selected original source photo URL absent from frozen blind sample")
        self.acq=joined

    def get_batch(self,part,batch,size):
        if (part,batch,size)!=(0,0,64):
            raise ValueError("Only the first bounded original geographic pair is permitted")
        return self.worker,self.acq,1


def measure_unit(source_unit:Path,outdir:Path,*,method:Path,detector:Path,sam:Path,
                 loader=None,evaluator=None)->dict:
    receipt=read_doc(source_unit/"first_matched_pair_precolour_source_receipt.json")
    if receipt.get("schema")!=PHOTO_SCHEMA or receipt.get("status")!=\
       "ORIGINAL_MONTHLY_SOURCE_GEOGRAPHIC_PAIR_PREPARED_NO_PIXELS":
        raise ValueError("No valid original monthly photo pair receipt")
    for filename,key in (
        ("first_matched_pair_original_blind_worker.csv","original_blind_worker_sha256"),
        ("first_matched_pair_original_photo_source_acquisition.csv","original_source_acquisition_sha256"),
        ("first_matched_pair_original_source_join_key.csv","original_biological_join_sha256")):
        if sha256(source_unit/filename)!=receipt.get(key):
            raise ValueError("Original pre-colour pair metadata altered after selecting the ROI work unit")
    test_only=loader is not None or evaluator is not None
    if test_only and (loader is None or evaluator is None):
        raise ValueError("Synthetic evaluator and loader must be supplied together")
    if test_only != (receipt["source_image_IDs_from_verified_monthly_20260927"] is False):
        raise ValueError("Production original image pixels cannot use a synthetic provenance gate")
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError("Original photographic measurement output already sealed; no image reprocessing")
    idx=PairSourceIndex(source_unit)
    report=run_batch(source_unit,0,0,64,outdir,
                     source=method,detector=detector,sam=sam,
                     indexed_source=idx,loader=loader,evaluator=evaluator,
                     _test_only_skip_method_check=test_only)
    source_rows=read_csv(source_unit/"first_matched_pair_original_source_join_key.csv")
    terminals=pd.read_csv(outdir/"terminal_p000_b0000000.csv",dtype=str).fillna("")
    details=pd.DataFrame(source_rows,dtype=str)
    if len(terminals)!=len(details) or terminals.measurement_id.duplicated().any():
        raise ValueError("Measured original photo terminal IDs do not cover prepared pair")
    joined=details.merge(terminals,on="measurement_id",validate="one_to_one")
    if len(joined)!=len(details) or joined.measurement_id.duplicated().any():
        raise ValueError("Original individual photo terminal/source biological join failed")
    merged=outdir/"first_original_geographically_matched_photo_terminal.csv"
    if report["n_source_photo_ids"]!=receipt["number_original_photo_IDs_in_first_unit"]:
        raise ValueError("Original source photo acquisition/ROI count contradicted the source")
    if joined.measurement_status.eq("classified_four_state_morph").sum()!=report["n_classifiable"]:
        raise ValueError("Original ROI coarse flower status contradicted exact terminal receipt")
    # A white ray + yellow disc capitulum can produce multiple photographic
    # coarse states WITHOUT two biological flower-colour morphs. Report this
    # separately from among-individual colour differences in each species.
    white=pd.to_numeric(joined["flower_fraction_white"],errors="coerce")
    yellow=pd.to_numeric(joined["flower_fraction_yellow"],errors="coerce")
    classifiable=joined.measurement_status.eq("classified_four_state_morph")
    join_mixed=classifiable & white.ge(.10) & yellow.ge(.10)
    joined["same_roi_white_and_yellow_both_ge10pct"]=join_mixed
    # Broad vs narrow colour differences can be confounded by photographs
    # made in different seasons. The date is not used to pick pixels or
    # assign colours and never enters the ROI evaluation process. Keep its
    # presence and coverage visible in the post-colour comparison receipt.
    def source_month(raw):
        value=str(raw or "").strip()
        try:
            parsed=date.fromisoformat(value)
        except ValueError:
            return ""
        return parsed.strftime("%Y-%m")

    joined["source_observation_month"]=joined["observed_on"].map(source_month)
    joined.to_csv(merged,index=False)

    # Describe original-photo colour outcomes by sampled within-species
    # geographic cells. Do not promote isolated images or same-ROI organ
    # colours to geographic polymorphism or environmental selection.
    cell_fields=(
        "inat_taxon_id","species","geographic_comparison_group",
        "latitude_5degree_cell_index","longitude_5degree_cell_index",
        "n_measured_original_photo_ids","n_distinct_original_observations",
        "n_distinct_original_observers","n_original_photo_ROIs_classifiable",
        "n_original_photo_ROIs_unclassifiable",
        "n_same_ROI_white_and_yellow_mixed_ge10pct",
        "n_white_photo_labels","n_yellow_orange_photo_labels",
        "n_red_pink_photo_labels","n_blue_purple_photo_labels",
        "original_observation_months_json","original_photo_IDs_sha256",
        "is_genetic_or_population_FCP_test"
    )
    by_cell=defaultdict(list)
    for row in joined.to_dict("records"):
        loc=geographical_cell(row)
        by_cell[(str(row["inat_taxon_id"]),
                 str(row["geographic_comparison_group"]),loc)].append(row)
    site_records=[]
    for key,records in sorted(by_cell.items()):
        tid,group,cell=key
        if len({r["species"] for r in records})!=1:
            raise ValueError("One geographic original-photo cell mixed source species")
        n_good=sum(r["measurement_status"]=="classified_four_state_morph"
                   for r in records)
        row={
            "inat_taxon_id":tid,"species":records[0]["species"],
            "geographic_comparison_group":group,
            "latitude_5degree_cell_index":cell[0],
            "longitude_5degree_cell_index":cell[1],
            "n_measured_original_photo_ids":len(records),
            "n_distinct_original_observations":
                len({r["observation_id"] for r in records}),
            "n_distinct_original_observers":
                len({r["observer_id"] for r in records}),
            "n_original_photo_ROIs_classifiable":n_good,
            "n_original_photo_ROIs_unclassifiable":len(records)-n_good,
            "n_same_ROI_white_and_yellow_mixed_ge10pct":
                sum(bool(r["same_roi_white_and_yellow_both_ge10pct"])
                    for r in records),
            "n_white_photo_labels":sum(r["morph"]=="white" for r in records),
            "n_yellow_orange_photo_labels":
                sum(r["morph"]=="yellow_orange" for r in records),
            "n_red_pink_photo_labels":
                sum(r["morph"]=="red_pink" for r in records),
            "n_blue_purple_photo_labels":
                sum(r["morph"]=="blue_purple" for r in records),
            "original_observation_months_json":
                json.dumps(sorted({str(r["source_observation_month"])
                                  for r in records if r["source_observation_month"]}),
                           separators=(",",":")),
            "original_photo_IDs_sha256":hashlib.sha256(
                "\n".join(sorted({str(r["photo_id"]) for r in records})).encode()
            ).hexdigest(),
            "is_genetic_or_population_FCP_test":False,
        }
        if (sum(row[f"n_{x}_photo_labels"] for x in
                ("white","yellow_orange","red_pink","blue_purple"))!=n_good):
            raise ValueError("Source geographic cell coarse colour labels disagree with ROI statuses")
        site_records.append(row)
    if (sum(r["n_measured_original_photo_ids"] for r in site_records)!=len(joined)
        or sum(r["n_original_photo_ROIs_classifiable"] for r in site_records)
           !=int(joined.measurement_status.eq("classified_four_state_morph").sum())):
        raise ValueError("Source cell coverage omitted/reassigned original photographed sites")
    site_file=outdir/"first_original_geographic_5degree_cell_photo_colour_QC.csv"
    with site_file.open("x",newline="",encoding="utf-8") as file:
        writer=csv.DictWriter(file,fieldnames=cell_fields,lineterminator="\n")
        writer.writeheader()
        writer.writerows(site_records)
    groups={}
    for name,group in joined.groupby("geographic_comparison_group"):
        known={m for m in group.source_observation_month if m}
        groups[str(name)]={
            "n_original_images_measured":len(group),
            "n_newly_four_state_classifiable":
                int(group.measurement_status.eq("classified_four_state_morph").sum()),
            "n_original_provider_photo_acquisition_failed":
                int(group.measurement_status.eq("image_acquisition_failed").sum()),
            "n_same_original_roi_white_and_yellow_both_ge10pct":
                int(group.same_roi_white_and_yellow_both_ge10pct.sum()),
            "n_distinct_source_observers_in_measured_photos":
                int(group.observer_id.nunique()),
            "n_distinct_source_observations_in_measured_photos":
                int(group.observation_id.nunique()),
            "n_source_photos_with_valid_original_observation_month":
                int(group.source_observation_month.ne("").sum()),
            "n_distinct_original_observation_calendar_months":
                len(known),
            "n_distinct_original_observation_years":
                len({m[:4] for m in known}),
            "source_calendar_months_with_photographs_sorted":
                sorted(known),
            "n_distinct_geographic_5degree_cells_in_measured_photos":
                len({geographical_cell(r) for r in group.to_dict("records")}),
        }
    if set(groups)!={"broad","narrow"}:
        raise ValueError("Broad/narrow species group lost in the source-to-pixel join")
    shared=sorted(set(groups["broad"]["source_calendar_months_with_photographs_sorted"]) &
                  set(groups["narrow"]["source_calendar_months_with_photographs_sorted"]))
    result={
        "schema":"fcp_rgfca_genuine_original_monthly_geographic_pair_ROI64_v1",
        "status":("SYNTHETIC_ORIGINAL_GEOGRAPHIC_PAIR_TEST_NOT_BIOLOGY" if test_only
                  else "FIRST_GENUINE_MONTHLY_ORIGINAL_GEO_MATCHED_PIXELS_MEASURED_NOT_1000_SPECIES"),
        "source_snapshot_date":SOURCE_DATE,
        "precolour_source_unit_receipt_sha256":sha256(
            source_unit/"first_matched_pair_precolour_source_receipt.json"),
        "new_original_photo_pixels_processed":report["n_source_photo_ids"],
        "n_classifiable":report["n_classifiable"],
        "by_geographic_observed_source_group":groups,
        "shared_source_calendar_months_between_broad_and_narrow":
            shared,
        "geography_vs_season_sampling_confounding_reported":True,
        "original_date_not_used_for_precolour_photo_replacement":True,
        "original_geographic_5degree_cell_photo_QC_csv_sha256":sha256(site_file),
        "n_original_geographic_5degree_cells_with_photo_terminals":len(site_records),
        "geographic_site_photo_label_table_is_descriptive_not_FCP_test":True,
        "n_original_image_acquisition_failed":report["n_acquisition_failed"],
        "original_pixel_terminal_sha256":report["terminal_csv_sha256"],
        "original_geographic_source_join_terminal_sha256":sha256(merged),
        "original_unmodified_roi_method_fingerprint":report["method"],
        "full_RGFCAROI1000_complete":False,
        "confirmed_genetic_FCP_species":0,
        "measured_geographic_effect_estimates":0,
        "old_measured_colours_imported":False,
    }
    (outdir/"first_geo_pair_real_pixel_receipt.json").write_text(
        json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("FCP_ORIGINAL_GEO_PAIR_FIRST_PIXELS",json.dumps(
        {k:result[k] for k in ("status","new_original_photo_pixels_processed",
                               "n_classifiable","by_geographic_observed_source_group")},sort_keys=True),
          flush=True)
    return result


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dated-source-artifact-root",required=True,type=Path)
    p.add_argument("--geographic-matched-pair-result-root",required=True,type=Path)
    p.add_argument("--prepared-first-original-source-unit",required=True,type=Path)
    p.add_argument("--measurement-outdir",required=True,type=Path)
    p.add_argument("--original-measure-script",required=True,type=Path)
    p.add_argument("--original-detector",required=True,type=Path)
    p.add_argument("--original-efficient-sam",required=True,type=Path)
    args=p.parse_args(argv)
    source=select_unit(args.dated_source_artifact_root,
                       args.geographic_matched_pair_result_root,
                       args.prepared_first_original_source_unit)
    result=measure_unit(args.prepared_first_original_source_unit,args.measurement_outdir,
                       method=args.original_measure_script,
                       detector=args.original_detector,sam=args.original_efficient_sam)
    print(json.dumps({"source":source,"result":result},sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
