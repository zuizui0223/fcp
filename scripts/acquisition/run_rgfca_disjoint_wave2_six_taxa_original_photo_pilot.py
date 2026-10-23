#!/usr/bin/env python3
"""Second DISJOINT six-species real photo ROI pilot, not the 1000-taxon result.

The original 1000 species were chosen from the independently audited
2026-09-27 Open Dataset without looking at flower colour. This pilot uses
a SEPARATE, explicitly labeled bounded LIVE iNaturalist API source as
of 2026-10-16, restricted to observations created by 2026-09-27.
It MUST NOT be admitted as an exact 2026-09-27 archived source snapshot
or used to claim RGFCA biological population-level conclusions.

Six NEW taxa (ranks 3-4 within each original photo-opportunity tier) are
selected by the SAME immutable hash, disjoint from the first six taxa; up to two API pages/taxon and
48+48 original photos/taxon are drawn before image opening. The same
historical ROI-v4/EfficientSAM byte-frozen original worker is used, and
photo bytes are discarded after work. No other photo/colour data are read.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

import pandas as pd

from fcp_pipeline.fcp_allphotos import sha256_file
from fcp_pipeline.fcp_fresh_study import BIO, MORPHS, WORKER
from fcp_pipeline.photo_polymorphism_screen import observer_split
from scripts.acquisition import collect_fcp_uncapped_metadata as collector
from scripts.acquisition.run_fcp_fresh_image_partition import (
    PermanentOriginalPhotoError, image_loader, verify_frozen_method,
)

ROSTER_SHA = "88aba565e1bdeef48ad3dc5eec0121c52448fda434aca849e542bb615f9d547c"
SOURCE_DATE = "2026-09-27"
PILOT_SALT = "FCP_RGFCA_1000_PRECOLOUR_LIVE_API_PILOT_20261016_V1"
TIERS = ("shallow_128_255", "medium_256_1023", "deep_ge1024")
TAXA_PER_TIER = 2
PAGES_PER_TAXON = 2
PHOTOS_PER_FOLD = 48
SAMPLE_FIELDS = ("inat_taxon_id", "species", "source_opportunity_tier",
                 "photo_id", "observation_id", "observer_id",
                 "latitude", "longitude", "photo_url_large", "photo_license",
                 "provider_photo_url", "screen_fold", "sample_hash")
FINAL_FIELDS = (*SAMPLE_FIELDS, "measurement_id", "morph", "measurement_status",
                "global_classifiable", "image_bytes_sha256",
                *(f"flower_fraction_{x}" for x in BIO))


def sampling_hash(*values: object) -> str:
    return hashlib.sha256(("|".join([PILOT_SALT, *(str(v) for v in values)])).encode()).hexdigest()


def pilot_taxa(roster: Path, *, n_per_tier: int = TAXA_PER_TIER,
               _synthetic: bool = False) -> list[dict]:
    if not 1 <= n_per_tier <= TAXA_PER_TIER:
        raise ValueError("Cannot change independently limited RGFCA real pilot taxon number")
    if not _synthetic and sha256_file(roster) != ROSTER_SHA:
        raise ValueError("1,000-species real original-source RGFCA roster changed")
    selected = []
    seen = set()
    with roster.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        required = {"inat_taxon_id", "species", "source_photo_opportunities",
                    "opportunity_stratum"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("No admitted source-count-only RGFCA species cohort")
        rows = list(reader)
    if not _synthetic and len(rows) != 1000:
        raise ValueError("RGFCA production 1,000-species roster incomplete")
    for label in TIERS:
        candidates = [r for r in rows if r["opportunity_stratum"] == label]
        candidates.sort(key=lambda r: (sampling_hash("taxon", r["inat_taxon_id"]), int(r["inat_taxon_id"])))
        if len(candidates) < 2 * n_per_tier:
            raise ValueError("Source metadata tier is underpopulated for disjoint replication")
        for r in candidates[n_per_tier:2 * n_per_tier]:
            tid = str(r["inat_taxon_id"])
            if not tid.isdecimal() or int(tid) <= 0 or int(tid) in seen or not r["species"].strip():
                raise ValueError("Duplicate or malformed source taxon in pilot")
            if not _synthetic and int(r["source_photo_opportunities"]) < 128:
                raise ValueError("Pilot source taxon was not in original >=128 photo opportunity set")
            seen.add(int(tid))
            selected.append({
                "inat_taxon_id":int(tid), "species":r["species"],
                "source_opportunity_tier":label,
                "original_20260927_source_photo_count":int(r["source_photo_opportunities"]),
            })
    from scripts.acquisition.run_rgfca_real_1000_enrolled_six_photo_pilot import pilot_taxa as first_wave_taxa
    previous = {x["inat_taxon_id"] for x in first_wave_taxa(roster, _synthetic=_synthetic)}
    if previous & {x["inat_taxon_id"] for x in selected}:
        raise ValueError("Replication six species overlap original technical wave")
    return selected


def _trusted_s3_url(pid: str, provider_url: str) -> str | None:
    # The actual source API supplies the photo ID and extension. Only construct
    # the pre-existing official Open Data bucket large-image object schema when
    # the source URL itself proves the ID and extension. Reject unknown hosts.
    from urllib.parse import urlsplit
    u = urlsplit(provider_url)
    if (u.scheme != "https" or u.hostname not in
        {"inaturalist-open-data.s3.amazonaws.com", "static.inaturalist.org", "inaturalist-open-data.s3.us-east-1.amazonaws.com"}
        or u.query or u.fragment):
        return None
    match = re.fullmatch(r"/photos/([0-9]+)/(?:square|small|medium|large|thumb)\.([A-Za-z0-9]+)",u.path)
    if not match or match.group(1) != str(pid) or match.group(2).lower() not in ("jpg","jpeg","png"):
        return None
    ext = match.group(2)
    return f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.{ext}"


def select_photo_rows(taxon: dict, raw: list[dict], *,
                      per_fold: int = PHOTOS_PER_FOLD) -> list[dict]:
    tid = taxon["inat_taxon_id"]
    tier = taxon["source_opportunity_tier"]
    unique: dict[str,dict] = {}
    collisions = set()
    for obj in raw:
        pid = str(obj["photo_id"])
        if (int(obj["inat_taxon_id"]) != tid or obj["species"] != taxon["species"]
            or not pid.isdecimal()):
            raise ValueError("Source API photo swapped away from preselected taxon")
        source_url = str(obj.get("photo_url_large") or "")
        s3 = _trusted_s3_url(pid, source_url)
        if s3 is None:
            continue
        record = {
            "inat_taxon_id":tid, "species":taxon["species"],
            "source_opportunity_tier":tier, "photo_id":pid,
            "observation_id":str(obj["observation_id"]),
            "observer_id":str(obj["observer_id"]),
            "latitude":obj["latitude"],"longitude":obj["longitude"],
            "photo_url_large":s3,"photo_license":obj["photo_license"],
            "provider_photo_url":source_url,
        }
        if pid in unique and unique[pid] != record:
            collisions.add(pid)
        else:
            unique[pid] = record
    for pid in collisions:
        unique.pop(pid,None)
    chosen = []
    for fold in (0,1):
        items = sorted(
            (p for p in unique.values()
             if observer_split(tid,p["observer_id"]) == fold),
            key=lambda x: (sampling_hash("original_photo", tid, x["photo_id"]),int(x["photo_id"])))
        used_obs, observer_counts = set(),Counter()
        for obj in items:
            obs = obj["observation_id"]
            observer = obj["observer_id"]
            if obs in used_obs or observer_counts[observer] >= 8:
                continue
            used_obs.add(obs)
            observer_counts[observer] += 1
            obj = dict(obj)
            obj["screen_fold"] = fold
            obj["sample_hash"] = sampling_hash("original_photo",tid,obj["photo_id"])
            chosen.append(obj)
            if sum(x["screen_fold"] == fold for x in chosen) == per_fold:
                break
    if len({r["photo_id"] for r in chosen}) != len(chosen):
        raise ValueError("One original source photo ID repeated across folds")
    return chosen


def precolour_metadata(roster: Path, output: Path, *,
                       request=None, pages: int = PAGES_PER_TAXON,
                       _synthetic=False) -> dict:
    if pages != PAGES_PER_TAXON and not _synthetic:
        raise ValueError("Pilot API source-page exposure was predeclared as exactly two")
    if not 1 <= pages <= PAGES_PER_TAXON:
        raise ValueError("Pilot source page resource limit invalid")
    if output.exists() and any(output.iterdir()):
        raise ValueError("Source-only pilot metadata output already exists")
    output.mkdir(parents=True, exist_ok=True)
    taxa = pilot_taxa(roster,_synthetic=_synthetic)
    req = request or (lambda params: collector.network_fetch(params,pause_seconds=1.2))
    all_photos = []
    per_taxon = []
    for record in taxa:
        tid = record["inat_taxon_id"]
        cursor = 0
        candidates = []
        n_api_obs = 0
        for page in range(pages):
            params = collector.query_parameters(tid,SOURCE_DATE,cursor)
            response = req(params)
            results = response.get("results")
            if not isinstance(results,list) or len(results)>collector.PAGE_SIZE:
                raise ValueError("Provider returned malformed bounded API results")
            ids = [int(o["id"]) for o in results]
            if ids != sorted(set(ids)) or any(x <= cursor for x in ids):
                raise ValueError("Nonmonotonic real original API source cursor")
            for obs in results:
                candidates.extend(collector.observation_rows(
                    obs,species=record["species"],taxon_id=tid,
                    freeze=SOURCE_DATE,query_label="bounded_live_API_technical_pilot_NOT_dated_archive"))
            n_api_obs += len(ids)
            if ids:
                cursor = ids[-1]
            if len(results) < collector.PAGE_SIZE:
                break
        selected = select_photo_rows(record,candidates)
        per_taxon.append({**record,"api_observations_read":n_api_obs,
                          "candidate_licensed_photo_metadata_rows":len(candidates),
                          "selected_original_source_photo_IDs":len(selected),
                          "selected_fold0":sum(r["screen_fold"]==0 for r in selected),
                          "selected_fold1":sum(r["screen_fold"]==1 for r in selected),
                          "pagination_exhaustive_for_taxon":n_api_obs<pages*collector.PAGE_SIZE,
                          "new_colour_photos_measured":0})
        all_photos.extend(selected)
        print("RGFCA_REAL_API_PRECOLOUR_METADATA",json.dumps(per_taxon[-1],sort_keys=True),flush=True)
    selected = output/"rgfca_six_wave2_source_photo_candidates_before_colour.csv"
    with selected.open("w",newline="",encoding="utf-8") as f:
        w = csv.DictWriter(f,fieldnames=SAMPLE_FIELDS,lineterminator="\n")
        w.writeheader()
        w.writerows(all_photos)
    doc = {
        "schema":"fcp_rgfca_1000_enrolled_real_original_API_disjoint_second_six_taxon_pilot_v1",
        "status":"SELECTED_ORIGINAL_PROVIDER_API_PHOTO_IDS_NO_IMAGE_PIXELS_YET",
        "source_ROSTER_sha256":sha256_file(roster),
        "source_candidate_universe_species":1000 if not _synthetic else None,
        "pilot_taxa":len(taxa),"pages_per_taxon_upper_bound":pages,
        "maximum_source_API_requests":len(taxa)*pages,
        "source_api_created_d2_filter":SOURCE_DATE,
        "provider_API_read_live_at_run_time":True,
        "EXACT_20260927_ARCHIVE_PHOTO_ID_MEMBERSHIP_VERIFIED":False,
        "NOT_THE_IMMUTABLE_20260927_METADATA_CENSUS":True,
        "original_photo_selection_pre_colour":True,
        "photo_selection_used_literature_phenotype":False,
        "photo_selection_used_climate_soil_or_prior_photo_colours":False,
        "max_selected_photo_per_taxon":PHOTOS_PER_FOLD*2,
        "n_selected_original_photo_IDs":len(all_photos),
        "selected_metadata_sha256":sha256_file(selected),
        "image_bytes_opened":False,
        "original_frozen_ROI_v4_executed":False,
        "RGFCA_200_realization_photo_polymorphism_Ppoly_computed":False,
        "main_1000_taxon_original_dated_source_completed":False,
        "technical_only_disjoint_previous_six_taxa":True,
        "wave":2,
        "per_taxon":per_taxon,
    }
    (output/"rgfca_six_wave2_source_precolour_receipt.json").write_text(
        json.dumps(doc,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    return doc


def measure(source: Path, output: Path, original: Path, detector: Path,
            sam: Path) -> dict:
    if output.exists() and any(output.iterdir()):
        raise ValueError("Original ROI technical pilot output already exists")
    doc = json.loads((source/"rgfca_six_wave2_source_precolour_receipt.json").read_text(encoding="utf-8"))
    meta = source/"rgfca_six_wave2_source_photo_candidates_before_colour.csv"
    if (doc.get("status") != "SELECTED_ORIGINAL_PROVIDER_API_PHOTO_IDS_NO_IMAGE_PIXELS_YET"
        or doc.get("selected_metadata_sha256") != sha256_file(meta)
        or doc.get("EXACT_20260927_ARCHIVE_PHOTO_ID_MEMBERSHIP_VERIFIED") is not False
        or doc.get("original_photo_selection_pre_colour") is not True
        or doc.get("RGFCA_200_realization_photo_polymorphism_Ppoly_computed") is not False):
        raise ValueError("Unverified or pseudo-frozen photo selection cannot enter original ROI worker")
    with meta.open("r",newline="",encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if tuple(reader.fieldnames or ()) != SAMPLE_FIELDS:
            raise ValueError("Pilot source selection columns changed")
        candidates = list(reader)
    if len(candidates)!=doc["n_selected_original_photo_IDs"] or not candidates:
        raise ValueError("Empty/incomplete original source photo pilot")
    if len({x["photo_id"] for x in candidates}) != len(candidates):
        raise ValueError("Duplicate source photo identities cannot be remeasured")
    fingerprint = verify_frozen_method(original,detector,sam)
    output.mkdir(parents=True, exist_ok=True)
    failures = {}
    checksums = {}
    blind = []
    with tempfile.TemporaryDirectory(prefix="fcp_rgfca_six_photo_roi_") as scratch:
        images=Path(scratch)/"images"
        images.mkdir()
        for n,item in enumerate(candidates):
            pid = item["photo_id"]
            mid = "RGFCA-PILOT-"+hashlib.sha256((PILOT_SALT+"|"+pid).encode()).hexdigest().upper()[:24]
            item["measurement_id"] = mid
            name = mid+".jpg"
            try:
                data=image_loader(item["photo_url_large"])
                (images/name).write_bytes(data)
                checksums[mid]=hashlib.sha256(data).hexdigest()
                blind.append({"measurement_id":mid,"image_filename":name,
                              "photo_license":item["photo_license"]})
            except PermanentOriginalPhotoError as exc:
                failures[mid]=f"{type(exc).__name__}: {str(exc)[:150]}"
            if (n+1)%24==0 or n+1==len(candidates):
                print("FCP_RGFCA_REAL_ORIGINAL_PHOTO_DOWNLOAD",
                      json.dumps({"processed":n+1,"total":len(candidates),
                                  "downloaded":len(blind),"missing":len(failures)}),flush=True)
        raw = Path(scratch)/"original_roi_v4_blind.csv"
        if blind:
            frame = Path(scratch)/"blind_original_photo_worker.csv"
            with frame.open("w",newline="",encoding="utf-8") as f:
                writer=csv.DictWriter(f,fieldnames=list(WORKER),lineterminator="\n")
                writer.writeheader()
                writer.writerows(blind)
            subprocess.run([
                sys.executable,str(original),"--worker-manifest",str(frame),
                "--images-dir",str(images),"--detector-weight",str(detector),
                "--efficient-sam-dir",str(sam),"--output-csv",str(raw)],check=True)
            results=pd.read_csv(raw,dtype={"measurement_id":str}).fillna("")
        else:
            results=pd.DataFrame(columns=["measurement_id","morph","measurement_status",
                                         *(f"flower_fraction_{x}" for x in BIO)])
        if len(results)!=len(blind) or results.measurement_id.duplicated().any():
            raise ValueError("Original ROI-v4 lost or duplicated source-blind photo ids")
        outcomes = {x["measurement_id"]:x for x in results.to_dict(orient="records")}
        if set(outcomes)!={x["measurement_id"] for x in blind}:
            raise ValueError("Original ROI-v4 returned different source-selected photo identities")
        terminal = []
        for item in candidates:
            mid=item["measurement_id"]
            result=outcomes.get(mid)
            if result is None:
                if mid not in failures:
                    raise ValueError("Original selected photo missing from both ROI output and acquisition failures")
                result={"measurement_id":mid,"morph":"mixed_uncertain",
                        "measurement_status":"image_acquisition_failed",
                        **{f"flower_fraction_{x}":0.0 for x in BIO}}
            good=result["measurement_status"]=="classified_four_state_morph"
            if good!=(result["morph"] in MORPHS):
                raise ValueError("Original photo classified status inconsistent with frozen four coarse states")
            terminal.append({**item,"morph":result["morph"],
                             "measurement_status":result["measurement_status"],
                             "global_classifiable":str(bool(good)).lower(),
                             "image_bytes_sha256":checksums.get(mid,""),
                             **{f"flower_fraction_{x}":result[f"flower_fraction_{x}"] for x in BIO}})
    terminal_csv=output/"rgfca_original_new_photo_six_taxa_wave2_live_API_technical.csv"
    with terminal_csv.open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=FINAL_FIELDS,lineterminator="\n")
        writer.writeheader()
        writer.writerows(terminal)
    by_taxon={}
    for item in terminal:
        tid=str(item["inat_taxon_id"])
        rec=by_taxon.setdefault(tid,{"n_source_photos_in_technical_pilot":0,
                                     "n_classifiable":0,"n_observer_fold0":0,
                                     "n_observer_fold1":0,"classifiable_four_morph_counts":{}})
        rec["n_source_photos_in_technical_pilot"]+=1
        rec["n_classifiable"]+=item["global_classifiable"]=="true"
        rec["n_observer_fold0"]+=str(item["screen_fold"])=="0"
        rec["n_observer_fold1"]+=str(item["screen_fold"])=="1"
        if item["global_classifiable"]=="true":
            m=item["morph"]
            rec["classifiable_four_morph_counts"][m]=rec["classifiable_four_morph_counts"].get(m,0)+1
    summary={
        "schema":"fcp_rgfca_new_original_ROI_v4_live_API_technical_second_disjoint_pilot_v1",
        "status":"REAL_NEW_ORIGINAL_PIXELS_DISJOINT_WAVE2_TECHNICAL_NOT_RGFCA_PPOLY",
        "selected_real_original_photo_IDs":len(candidates),
        "original_pixels_downloaded":len(blind),
        "original_pixels_ROI_v4_evaluated":len(results),
        "n_classifiable_original_flower_photos":sum(x["global_classifiable"]=="true" for x in terminal),
        "number_taxa_selected":doc["pilot_taxa"],
        "source_API_selection_receipt_sha256":sha256_file(source/"rgfca_six_wave2_source_precolour_receipt.json"),
        "source_image_ID_selection_sha256":sha256_file(meta),
        "new_terminal_photo_csv_sha256":sha256_file(terminal_csv),
        "model_attestation":fingerprint,
        "prior_FCP_measured_colours_loaded":False,
        "exact_dated_20260927_source_membership_admitted":False,
        "source_kind":"iNaturalist_API_created_before_20260927_live_snapshot_NOT_dated_monthly_archive",
        "RGFCA_200x250_Ppoly_not_estimated":True,
        "full_1000_species_ROI_study_complete":False,
        "technical_only_disjoint_previous_six_taxa":True,
        "wave":2,
        "image_bytes_persisted":False,
        "by_taxon":by_taxon,
    }
    (output/"rgfca_original_new_photo_six_taxa_wave2_live_API_technical_receipt.json").write_text(
        json.dumps(summary,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("FCP_RGFCA_ACTUAL_NEW_ORIGINAL_ROI_PIXEL_PILOT",
          json.dumps({k:summary[k] for k in (
              "status","selected_real_original_photo_IDs",
              "original_pixels_downloaded","n_classifiable_original_flower_photos",
              "number_taxa_selected","by_taxon")},sort_keys=True),flush=True)
    return summary


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest="phase",required=True)
    p=sub.add_parser("source",help="Bounded live API; records source-photo IDs before opening pixels")
    p.add_argument("--rgfca-roster",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    p=sub.add_parser("measure",help="Original byte-pinned ROI-v4/EfficientSAM in blind photo worker")
    p.add_argument("--precolour-source",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    p.add_argument("--original-measure-script",required=True,type=Path)
    p.add_argument("--original-detector",required=True,type=Path)
    p.add_argument("--original-efficient-sam",required=True,type=Path)
    a=parser.parse_args(argv)
    if a.phase=="source":
        print(json.dumps(precolour_metadata(a.rgfca_roster,a.outdir),sort_keys=True))
    else:
        print(json.dumps(measure(a.precolour_source,a.outdir,a.original_measure_script,
                                 a.original_detector,a.original_efficient_sam),sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
