#!/usr/bin/env python3
"""Recover public point coordinates for ALL ORIGINAL FCP 42111+85337 photo records.

This is a metadata-only transport layer. No photo pixels or new colour labels.
Fixed photographed IDs only. Batch the public iNaturalist observation-ID API
in 100-ID requests, with rate limiting, frozen source audit and explicit
missingness. Exact source photos are never substituted.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import time
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError

import numpy as np
import pandas as pd

from audit_fcp_global42111_original_photo_coordinate_pilot_20261009 import (
    public_position,cell_id_for
)

SOURCE_RUN=37855050150
ORIGINAL_BREADTH=42111
ORIGINAL_TAXON_CELL=85337
BATCH=100
SHARDS=8
INTERVAL=1.2
TIMEOUT=45
ENDPOINT="https://api.inaturalist.org/v1/observations"
HEADERS={"User-Agent":"fcp-42111-original-taxoncells-georecovery/1.0","Accept":"application/json"}
SCHEMA="fcp_global42111_original_photo_full_public_coordinates_v1"


def sha(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""):h.update(block)
    return h.hexdigest()


def identity_tuple(row):
    return (int(row.inat_taxon_id),int(row.observation_id),int(row.photo_id))


def make_manifest(breadth:pd.DataFrame,cell:pd.DataFrame)->tuple[pd.DataFrame,dict]:
    fields=("inat_taxon_id","observation_id","photo_id")
    if len(breadth)!=ORIGINAL_BREADTH or len(cell)!=ORIGINAL_TAXON_CELL:
        raise RuntimeError("Previously measured source denominators drift")
    if not set(fields).issubset(breadth) or not set((*fields,"cell_id")).issubset(cell):
        raise RuntimeError("Original observation/photo identity missing")
    b=breadth[list(fields)].copy()
    b["cell_id"]=-1
    b["is_breadth_source"]=True
    b["is_taxon_cell_source"]=False
    c=cell[list(fields)+["cell_id"]].copy()
    c["is_breadth_source"]=False
    c["is_taxon_cell_source"]=True
    joint=pd.concat([b,c],ignore_index=True)
    for k in (*fields,"cell_id"):
        joint[k]=pd.to_numeric(joint[k],errors="raise").astype("int64")
    if (joint[list(fields)]<=0).any().any():
        raise RuntimeError("Frozen photo observation or taxon IDs nonpositive")
    nonnegative=joint.cell_id.ge(0)
    if not joint.loc[nonnegative,"cell_id"].between(0,161).all():
        raise RuntimeError("Taxon-cell source outside 162 global cells")
    # Every original photo is exactly one specimen taxon/observation identity.
    if joint.photo_id.nunique()!=joint[list(fields)].drop_duplicates().shape[0]:
        raise RuntimeError("Source photo ID maps to incompatible biological identity")
    if joint.groupby("observation_id").inat_taxon_id.nunique().gt(1).any():
        raise RuntimeError("Same original observed individual reused as different taxa")
    source=joint.groupby(list(fields),sort=True,as_index=False).agg(
        original_cell_id=("cell_id","max"),
        present_in_breadth=("is_breadth_source","max"),
        present_in_taxon_cell=("is_taxon_cell_source","max"),
    )
    if int(source.present_in_breadth.sum())!=ORIGINAL_BREADTH or int(source.present_in_taxon_cell.sum())!=ORIGINAL_TAXON_CELL:
        raise RuntimeError("Cannot reconstruct complete observed source usage")
    # Group all source photo rows into disjoint observation-ID retrieval shards.
    source["shard_index"]=source.observation_id.map(
        lambda i:int(hashlib.sha256(f"fcp42111_fixed_site_{int(i)}".encode()).hexdigest()[:8],16)%SHARDS
    ).astype(int)
    if not source.shard_index.between(0,SHARDS-1).all():
        raise RuntimeError("Invalid manifest shard rank")
    info={
        "schema":"fcp_global42111_full_original_photo_metadata_manifest_v1",
        "original_breadth_source_rows":ORIGINAL_BREADTH,
        "original_taxon_cell_source_rows":ORIGINAL_TAXON_CELL,
        "unique_original_photo_id_triples":len(source),
        "unique_original_observation_ids":int(source.observation_id.nunique()),
        "original_taxon_species":int(source.inat_taxon_id.nunique()),
        "shards":SHARDS,
        "source_colours_were_not_read":True,
        "new_photo_IDs_introduced":False,
        "no_photo_pixels_opened":True,
        "confirmatory_decisions_changed":False,
    }
    if info["original_taxon_species"]!=ORIGINAL_BREADTH:
        raise RuntimeError("Not every original FCP taxon retained")
    return source,info


def fetch_batch(ids:list[int])->dict:
    if not 1<=len(ids)<=BATCH or len(set(ids))!=len(ids):
        raise ValueError("Photo observation batch invalid")
    url=ENDPOINT+"?"+urlencode({"per_page":200,"id":",".join(str(i) for i in ids)})
    req=Request(url,headers=HEADERS)
    with urlopen(req,timeout=TIMEOUT) as response:
        obj=json.load(response)
    if not isinstance(obj,dict) or not isinstance(obj.get("results"),list):
        raise RuntimeError("Public observation source JSON missing results")
    return obj


def unpack(original:pd.DataFrame,returned:dict|None)->list[dict]:
    results=[]
    for r in original.itertuples(index=False):
        rec={
            "inat_taxon_id":int(r.inat_taxon_id),"observation_id":int(r.observation_id),
            "photo_id":int(r.photo_id),"original_cell_id":int(r.original_cell_id),
            "present_in_breadth":bool(r.present_in_breadth),
            "present_in_taxon_cell":bool(r.present_in_taxon_cell),
            "latitude":None,"longitude":None,"positional_accuracy_m":None,
            "site_geo_status":"SOURCE_OBSERVATION_NOT_RETURNED",
        }
        if returned is None:
            results.append(rec);continue
        tax=returned.get("taxon")
        if not isinstance(tax,dict) or int(tax.get("id") or 0)!=rec["inat_taxon_id"]:
            rec["site_geo_status"]="SOURCE_TAXON_IDENTITY_MISMATCH"
            results.append(rec);continue
        if not any(int(photo.get("id") or 0)==rec["photo_id"] for photo in returned.get("photos") or [] if isinstance(photo,dict)):
            rec["site_geo_status"]="EXACT_ORIGINAL_PHOTO_MISSING"
            results.append(rec);continue
        if bool(returned.get("captive")):
            rec["site_geo_status"]="CULTIVATED_NOT_WILD"
            results.append(rec);continue
        coord=public_position(returned)
        if coord is None:
            rec["site_geo_status"]="NO_UNOBSCURED_PUBLIC_COORDINATES"
            results.append(rec);continue
        lat,lon=coord
        if rec["original_cell_id"]>=0 and cell_id_for(lat,lon)!=rec["original_cell_id"]:
            rec["site_geo_status"]="CURRENT_PUBLIC_POSITION_OUTSIDE_HISTORICAL_EQUAL_AREA_CELL"
            results.append(rec);continue
        accuracy=returned.get("positional_accuracy")
        try: accuracy=float(accuracy)
        except (ValueError,TypeError):accuracy=None
        rec.update({
            "latitude":lat,"longitude":lon,
            "positional_accuracy_m":accuracy if accuracy is not None and np.isfinite(accuracy) else None,
            "site_geo_status":"VALID_PUBLIC_ORIGINAL_PHOTO_POINT",
        })
        results.append(rec)
    return results


def run_shard(source:pd.DataFrame,shard:int,outdir:Path,client=fetch_batch,sleep=time.sleep)->dict:
    if not 0<=shard<SHARDS:
        raise ValueError("Unknown original fixed shard")
    if not {"shard_index","observation_id","inat_taxon_id","photo_id","original_cell_id"}.issubset(source):
        raise RuntimeError("Source fixed manifest incomplete")
    subset=source.loc[source.shard_index==shard].copy()
    if len(subset)==0:
        raise RuntimeError("Frozen shard accidentally empty")
    byobs={int(obs):g for obs,g in subset.groupby("observation_id",sort=True)}
    ids=sorted(byobs)
    outdir.mkdir(parents=True,exist_ok=True)
    progress=outdir/f"original_geo_records_shard_{shard:02}.jsonl"
    if progress.exists():
        raise RuntimeError("Never silently re-fetch previously started source IDs")
    rows=[]
    failed_batches=0
    for offset in range(0,len(ids),BATCH):
        group=ids[offset:offset+BATCH]
        try:
            response=client(group)
            returned={int(z["id"]):z for z in response["results"] if isinstance(z,dict) and z.get("id")}
            if len(returned)>BATCH or not set(returned).issubset(set(group)):
                raise RuntimeError("Unrequested public observations returned")
            part=[r for i in group for r in unpack(byobs[i],returned.get(i))]
        except (HTTPError,URLError,TimeoutError,RuntimeError,ValueError) as err:
            failed_batches+=1
            part=[{
                **rec, "site_geo_status":"API_ERROR_UNRESOLVED",
                "error_class":type(err).__name__,
            } for i in group for rec in unpack(byobs[i],None)]
        rows.extend(part)
        with progress.open("a",encoding="utf-8") as f:
            for r in part:f.write(json.dumps(r,sort_keys=True)+"\n")
        if offset+BATCH<len(ids):sleep(INTERVAL)
    out=pd.DataFrame(rows)
    if len(out)!=len(subset) or out.photo_id.duplicated().any():
        raise RuntimeError("Shard lost or duplicated original photo identity")
    summary={
        "schema":SCHEMA,
        "type":"SHARD_METADATA_ONLY",
        "source_run_id":SOURCE_RUN,
        "shard_index":shard,"shards":SHARDS,
        "n_source_original_photos":len(subset),
        "n_original_observations_queried":len(ids),
        "n_100_ID_api_requests":int(np.ceil(len(ids)/BATCH)),
        "n_batch_errors":failed_batches,
        "n_exact_photo_public_geocoordinates":int(out.site_geo_status.eq("VALID_PUBLIC_ORIGINAL_PHOTO_POINT").sum()),
        "source_state_counts":{str(k):int(v) for k,v in out.site_geo_status.value_counts().items()},
        "colour_results_unopened":True,
        "images_unopened":True,
    }
    out.to_csv(outdir/f"original_geo_records_shard_{shard:02}.csv.gz",index=False,
               compression={"method":"gzip","compresslevel":9,"mtime":0})
    (outdir/f"receipt_shard_{shard:02}.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    print(json.dumps(summary,sort_keys=True),flush=True)
    return summary


def seal(manifest:pd.DataFrame,root:Path,outdir:Path)->dict:
    parts=[];receipts=[]
    for shard in range(SHARDS):
        rc=root/f"receipt_shard_{shard:02}.json"
        photo=root/f"original_geo_records_shard_{shard:02}.csv.gz"
        if not rc.exists() or not photo.exists():
            raise RuntimeError(f"Original metadata shard {shard} missing")
        z=json.loads(rc.read_text())
        if z["schema"]!=SCHEMA or z["shard_index"]!=shard or z["shards"]!=SHARDS:
            raise RuntimeError(f"Shard {shard} provenance changed")
        d=pd.read_csv(photo,low_memory=False)
        if len(d)!=z["n_source_original_photos"] or int(d.site_geo_status.eq("VALID_PUBLIC_ORIGINAL_PHOTO_POINT").sum())!=z["n_exact_photo_public_geocoordinates"]:
            raise RuntimeError(f"Shard {shard} receipt and observation table disagree")
        receipts.append(z);parts.append(d)
    out=pd.concat(parts,ignore_index=True)
    if len(out)!=len(manifest):
        raise RuntimeError("Source photograph ID denominator changed")
    if out.photo_id.duplicated().any() or set(out.photo_id.astype("int64"))!=set(manifest.photo_id.astype("int64")):
        raise RuntimeError("Some original measured photo IDs were replaced or omitted")
    expected=manifest[["photo_id","observation_id","inat_taxon_id","original_cell_id","present_in_breadth","present_in_taxon_cell"]]
    observed=out[list(expected.columns)]
    joined=observed.merge(expected,on="photo_id",validate="one_to_one",suffixes=("_a","_b"))
    for key in ("observation_id","inat_taxon_id","original_cell_id","present_in_breadth","present_in_taxon_cell"):
        if not (joined[key+"_a"].astype(str)==joined[key+"_b"].astype(str)).all():
            raise RuntimeError(f"Original {key} identity drift")
    counted=int(out.site_geo_status.eq("VALID_PUBLIC_ORIGINAL_PHOTO_POINT").sum())
    report={
        "schema":"fcp_global42111_original_public_geocoordinates_full_sealed_v1",
        "status":"ORIGINAL_PHOTO_PUBLIC_COORDINATE_METADATA_SEALED_NO_CLIMATE_TEST",
        "original_species_denominator":42111,
        "original_taxon_cell_denominator":85337,
        "n_original_distinct_photo_IDs":len(out),
        "n_distinct_original_observation_IDs":int(out.observation_id.nunique()),
        "n_orig_source_photo_public_point_geolocated":counted,
        "n_orig_source_photo_missing_or_unusable_point":len(out)-counted,
        "n_breadth_original_images_with_public_points":int((out.present_in_breadth & out.site_geo_status.eq("VALID_PUBLIC_ORIGINAL_PHOTO_POINT")).sum()),
        "n_taxon_cell_original_images_with_public_points":int((out.present_in_taxon_cell & out.site_geo_status.eq("VALID_PUBLIC_ORIGINAL_PHOTO_POINT")).sum()),
        "photo_metadata_status_counts":{str(k):int(v) for k,v in out.site_geo_status.value_counts().items()},
        "n_total_api_requests":sum(x["n_100_ID_api_requests"] for x in receipts),
        "n_total_batch_errors":sum(x["n_batch_errors"] for x in receipts),
        "source_shard_receipts":receipts,
        "never_infer_unclassified_flower_as_monomorphic":True,
        "photo_morph_labels_never_opened_for_selection":True,
        "new_image_pixels_opened":False,
        "actual_climate_soil_rasters_sampled":False,
        "confirmatory_decisions_changed":False,
    }
    outdir.mkdir(parents=True,exist_ok=True)
    out.to_csv(outdir/"all_original_measured_photo_public_coordinates.csv.gz",index=False,
               compression={"method":"gzip","compresslevel":9,"mtime":0})
    (outdir/"result.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps(report,indent=2,sort_keys=True),flush=True)
    return report


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--mode",required=True,choices=["manifest","shard","seal"])
    p.add_argument("--breadth",type=Path)
    p.add_argument("--taxon-cell",type=Path)
    p.add_argument("--manifest",type=Path)
    p.add_argument("--shard-index",type=int)
    p.add_argument("--shard-root",type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    a=p.parse_args()
    a.outdir.mkdir(parents=True,exist_ok=True)
    if a.mode=="manifest":
        if not a.breadth or not a.taxon_cell:raise ValueError("Original source required")
        s,r=make_manifest(pd.read_csv(a.breadth,low_memory=False),pd.read_csv(a.taxon_cell,low_memory=False))
        s.to_csv(a.outdir/"frozen_original_measured_photo_manifest.csv.gz",index=False,
                 compression={"method":"gzip","compresslevel":9,"mtime":0})
        (a.outdir/"result.json").write_text(json.dumps(r,indent=2,sort_keys=True)+"\n")
        print(json.dumps(r,sort_keys=True))
    elif a.mode=="shard":
        if not a.manifest or a.shard_index is None:raise ValueError("Frozen photo-ID manifest required")
        run_shard(pd.read_csv(a.manifest,low_memory=False),a.shard_index,a.outdir)
    else:
        if not a.manifest or not a.shard_root:raise ValueError("Unmodified manifest/shards required")
        seal(pd.read_csv(a.manifest,low_memory=False),a.shard_root,a.outdir)
if __name__=="__main__":
    main()
