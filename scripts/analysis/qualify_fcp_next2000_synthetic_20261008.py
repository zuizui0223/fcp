#!/usr/bin/env python3
"""Synthetic, outcome-blind end-to-end preflight for NEXT-2000 prospective FCP.

Feeds ONLY invented flower labels, localities and dates into a frozen historical
analysis engine. Does not download or inspect any of the newly selected photos.
Validates final JSON bytes and ZIP artifact packaging, to avoid a P500-style
late serialization failure.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

PREOPEN_MANIFEST_SHA256 = "e1301095533dfa755c1588cc550522bf31aabaef632083e265a138f0e3af6d2a"
PREOPEN_RECEIPT_SHA256 = None  # selected manifest SHA is the immutable auth key
HISTORICAL_ENGINE_COMMIT = "969c67d274868b5dc4e8c2c8c9fa5d1f1beb5527"
ENGINE_SOURCE_GIT_BLOB_SHA = "c7b277c249fc6121302b35ce2960e2003cfd0930"
ENGINE_RELATIVE_PATH = "scripts/analysis/run_polymorphism_specieswide_space_season_20261008.py"
N_SYNTHETIC_SPECIES = 12
N_SYNTHETIC_PHOTOS_PER_SPECIES = 80


def sha(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def open_engine(path: Path):
    spec=importlib.util.spec_from_file_location("fcp_frozen_geo_engine",path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Frozen spatial engine cannot be loaded")
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    for key in ("species_test","summarize","PERMUTATIONS","N_BOOTSTRAP"):
        if not hasattr(mod,key):
            raise RuntimeError(f"Frozen spatial engine lacks {key}")
    if int(mod.PERMUTATIONS)!=199 or int(mod.N_BOOTSTRAP)!=1999:
        raise RuntimeError("Historical calendar-null sampling constants have drifted")
    return mod


def synthetic_species(taxon_id:int, *, seasonal_only:bool=False,
                      one_observer:bool=False)->pd.DataFrame:
    n=N_SYNTHETIC_PHOTOS_PER_SPECIES
    lat=np.r_[np.full(40,10.0),np.full(40,12.0)]
    month=np.tile(np.r_[np.full(20,4),np.full(20,7)],2)
    labels=(np.where(month==4,"white","red_pink") if seasonal_only
            else np.array(["white"]*40+["red_pink"]*40))
    obs=(["the_same_observer"]*n if one_observer else
         [f"photographer_{k%10}" for k in range(n)])
    return pd.DataFrame({
        "inat_taxon_id":np.repeat(taxon_id,n),
        "species":[f"SyntheticGenus{taxon_id} example"]*n,
        "photo_id":np.arange(n)+taxon_id*n,
        "latitude":lat,
        "longitude":np.zeros(n),
        "morph":labels,
        "month":pd.Series(month,dtype="Int64"),
        "quarter":pd.Series((month-1)//3+1,dtype="Int64"),
        "year":pd.Series(np.full(n,2025),dtype="Int64"),
        "observer":obs,
    })


def preflight(engine, selected_manifest:Path, output_dir:Path) -> dict:
    if sha(selected_manifest)!=PREOPEN_MANIFEST_SHA256:
        raise RuntimeError("Prospective outcome-blind 2,000 species manifest hash drift")
    candidates=pd.read_csv(selected_manifest,sep="\t")
    if len(candidates)!=2000 or candidates.inat_taxon_id.nunique()!=2000:
        raise RuntimeError("A live prospective taxon identity gate failed")
    output_dir.mkdir(parents=True,exist_ok=True)
    checks={}
    allrows={}
    allnulls={}
    for policy in ("all","different_observer"):
        good=[];null=[]
        for i in range(N_SYNTHETIC_SPECIES):
            tid=900000000+i
            frame=synthetic_species(tid)
            result=engine.species_test(frame,"future_synthetic","month",policy)
            if result is None:
                raise RuntimeError(f"Synthetic spatial gradient unexpectedly not testable: {policy}")
            record,perm=result
            if not record["conditional_identifiable"] or record["excess_over_stratified_null"]<=0.2:
                raise RuntimeError("A true simulated spatial gradient did not survive month-null")
            if len(perm)!=199:
                raise RuntimeError("Null permutation count drift")
            good.append(record);null.append(perm)
        summary=engine.summarize(good,null,"future_synthetic","month",policy)
        if summary["n_species"]!=N_SYNTHETIC_SPECIES:
            raise RuntimeError("Synthetic species denominator drift")
        if not summary["estimable"] or summary["permutation_p_upper"]>0.05:
            raise RuntimeError("Positive synthetic geographic colour shift was not detected")
        if summary["species_bootstrap_95CI_excess"][0]<=0:
            raise RuntimeError("Synthetic whole-species bootstrap interval failed")
        checks[f"synthetic_spatial_clustering_{policy}"]=True
        allrows[policy]=good
        allnulls[policy]=null

    seasonal=synthetic_species(999999999,seasonal_only=True)
    season=engine.species_test(seasonal,"future_synthetic","month","all")
    if season is None or season[0]["conditional_identifiable"] or abs(season[0]["excess_over_stratified_null"])>1e-12:
        raise RuntimeError("Season-only false biological signal not held at null")
    checks["pure_season_switch_is_not_false_geography"]=True

    one=synthetic_species(888888888,one_observer=True)
    if engine.species_test(one,"future_synthetic","month","different_observer") is not None:
        raise RuntimeError("Same-observer contamination not excluded")
    checks["observer_disjoint_policy_blocks_single_photographer"]=True

    # Explicit reproducible source code presence, file input frozen by commit
    # in the workflow. No real colour measurements opened in any code path.
    reports={
        "schema":"fcp_next2000_synthetic_end_to_end_qualification_v1",
        "date_jst":"2026-10-08",
        "qualification_type":"SYNTHETIC_ONLY_NO_NEW_PHOTO_PIXEL_OPENING",
        "prospective_species_manifest_sha256":sha(selected_manifest),
        "selected_species_count":int(len(candidates)),
        "frozen_historical_engine_commit":HISTORICAL_ENGINE_COMMIT,
        "frozen_historical_engine_git_blob_sha":ENGINE_SOURCE_GIT_BLOB_SHA,
        "engine_permutations":int(engine.PERMUTATIONS),
        "engine_bootstrap":int(engine.N_BOOTSTRAP),
        "n_synthetic_species":N_SYNTHETIC_SPECIES,
        "n_synthetic_photo_records":N_SYNTHETIC_SPECIES*N_SYNTHETIC_PHOTOS_PER_SPECIES,
        "all_checks":checks,
        "synthetic_primary_effects":{
            policy:{
                "n_species":len(allrows[policy]),
                "mean_excess":float(np.mean([x["excess_over_stratified_null"] for x in allrows[policy]])),
                "result_is_biological_evidence":False}
            for policy in ("all","different_observer")
        },
        "outcome_firewall":{
            "new_photo_urls_fetched":False,
            "new_image_pixels_opened":False,
            "new_species_flower_colour_photos_classified":False,
            "real_new_colour_morph_outcome_used":False,
            "historical_posthoc_result_used_to_select_species":False,
            "new_prospective_confirmatory_verdict_computed":False
        },
        "technical_qualification_result":"PASS_SYNTHETIC_SERIALIZATION_AND_ARTIFACT_PACKAGING",
        "metadata_acquisition_authorized_by_this_receipt":False,
        "biological_pixel_opening_authorized_by_this_receipt":False,
        "confirmatory_decisions_changed":False
    }

    path=output_dir/"technical_result.json"
    path.write_text(json.dumps(reports,indent=2,sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")
    payload=json.loads(path.read_text(encoding="utf-8"))
    if payload["technical_qualification_result"]!=reports["technical_qualification_result"]:
        raise RuntimeError("Terminal serialized result did not validate after disk writing")
    zip_path=output_dir/"synthetic_technical_artifact.zip"
    with zipfile.ZipFile(zip_path,"w",compression=zipfile.ZIP_DEFLATED) as z:
        z.write(path,"technical_result.json")
    with zipfile.ZipFile(zip_path) as z:
        if z.namelist()!=["technical_result.json"]:
            raise RuntimeError("Artifact package missing terminal result")
        from_archive=json.loads(z.read("technical_result.json").decode("utf-8"))
    if from_archive!=payload:
        raise RuntimeError("Artifact archive data differ from validated disk result")
    reports["terminal_json_sha256"]=sha(path)
    reports["synthetic_artifact_zip_sha256"]=sha(zip_path)
    reports["artifact_reopened_and_validated"]=True
    (output_dir/"technical_execution_receipt.json").write_text(
        json.dumps(reports,indent=2,sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps(reports,indent=2,sort_keys=True,allow_nan=False))
    return reports


def main()->None:
    p=argparse.ArgumentParser()
    p.add_argument("--analysis-engine",required=True,type=Path)
    p.add_argument("--prospective-manifest",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    args=p.parse_args()
    m=open_engine(args.analysis_engine)
    preflight(m,args.prospective_manifest,args.outdir)


if __name__=="__main__":
    main()
