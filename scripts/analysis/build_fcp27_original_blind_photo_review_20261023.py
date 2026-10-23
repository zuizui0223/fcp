#!/usr/bin/env python3
"""FCP original-photo botanical review: 27 observed colour-mode candidates + 27 controls.

Selection is source-locked and deterministic. Reviewer packages with masked
identities never include inferred mode, morphology, species, photo IDs, dates,
observer IDs, locations or photographic palette features. Only a separate PI
mapping contains those fields. Human biological reviews are NOT inferred.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import os
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageOps
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "analysis"))
from fcp_full4475_hypothesis_census_20261021 import load_original, SOURCE_SHA
from fcp_photo_modes_geographic_shape_20261023 import mode_scores

MODE = "recurrent_photographic_mode_candidate"
CONTROL = "not_resolved_as_multimodal"
ALL_N = 447500
CLASS_N = 227107
MODE_N = 27
PAIR_N = 27
SELECTED_PER_SPECIES = 11
REVIEW_COLS = (
    "review_id", "review_order", "target_flower_visible",
    "petal_region_unambiguous", "observed_dominant_colour",
    "flower_stage", "lighting_white_balance",
    "multiple_flower_or_organ_confusion",
    "cultivation_or_taxon_concern", "review_confidence", "reviewer_notes"
)
SALT = "fcp27-mode-blind-photo-review-20261023-v1"


def need(test, message):
    if not test:
        raise RuntimeError("FCP27_BLIND_REVIEW: " + message)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for buf in iter(lambda: f.read(1 << 20), b""):
            h.update(buf)
    return h.hexdigest()


def hashed(*parts):
    return hashlib.sha256("|".join(map(str, (SALT, *parts))).encode()).hexdigest()


def match_species(mode, continuous, source):
    need(len(mode) == 1761 and mode.inat_taxon_id.is_unique, "modes artifact changed")
    need(len(continuous) == 1761 and continuous.taxon_id.is_unique,
         "continuous original species artifact changed")
    cohort = (source.groupby("inat_taxon_id", sort=False).fcp_source_cohort
              .first().rename("source_cohort"))
    join = mode.merge(continuous, left_on="inat_taxon_id",
                      right_on="taxon_id", how="inner", validate="one_to_one")
    need(join["species_x"].astype(str).eq(join["species_y"].astype(str)).all(),
         "source model vs continuous matched taxonomic name changed")
    join = (join.drop(columns="species_y").rename(columns={"species_x":"species"})
            .join(cohort, on="inat_taxon_id"))
    need(len(join) == 1761 and join.source_cohort.notna().all(),
         "source taxon join incomplete")
    positive = join.loc[join.mode_status == MODE].copy().sort_values("inat_taxon_id")
    candidates = join.loc[join.mode_status == CONTROL].copy().sort_values("inat_taxon_id")
    need(len(positive) == MODE_N and len(candidates) == 1691, "original source mode frequencies drift")
    # A negative control is not eligible merely because few geographic photo
    # rows made bimodality harder to detect. Match only controls passing the
    # same >=40 geo-valid photographic-mode detectability threshold.
    need(positive.n_spatial_classifiable.ge(40).all(), "positive mode source fit gate drift")
    candidates = candidates.loc[candidates.n_spatial_classifiable.ge(40)].copy()
    need(len(candidates) >= 100, "fewer than 100 comparable candidate controls")
    specs = [
        ("n_classified_original", 1., "log"),
        ("pairs_50km", 0.85, "log"),
        ("global_jsd", 1.3, "raw"),
        ("photographic_pc1_explained_variance", 0.9, "raw"),
    ]
    mat_a, mat_b = [], []
    for name,weight,transform in specs:
        a = positive[name].to_numpy(float)
        b = candidates[name].to_numpy(float)
        if transform == "log":
            a,b = np.log1p(a), np.log1p(b)
        lo,hi = np.nanquantile(np.concatenate([a,b]), [.10,.90])
        scale = max(float(hi-lo), .015)
        mat_a.append((a/scale)*weight)
        mat_b.append((b/scale)*weight)
    a = np.column_stack(mat_a)
    b = np.column_stack(mat_b)
    dist = np.sum((a[:,None,:]-b[None,:,:])**2, axis=2)
    dist += 0.8*(positive.photographic_coarse_label_count.to_numpy()[:,None]
                  !=candidates.photographic_coarse_label_count.to_numpy()[None,:])
    dist += 0.35*(positive.source_cohort.to_numpy()[:,None]
                   !=candidates.source_cohort.to_numpy()[None,:])
    row,col=linear_sum_assignment(dist)
    need(len(row) == PAIR_N and len(set(col)) == PAIR_N, "one-to-one control matching failed")
    positive = positive.iloc[row].copy()
    control = candidates.iloc[col].copy()
    positive["review_arm"] = "candidate"
    control["review_arm"] = "matched_control"
    positive["pair_no"] = np.arange(1,PAIR_N+1)
    control["pair_no"] = np.arange(1,PAIR_N+1)
    positive["match_distance"] = dist[row,col]
    control["match_distance"] = dist[row,col]
    combined = pd.concat([positive, control], ignore_index=True)
    need(len(combined) == 54 and combined.inat_taxon_id.is_unique, "wrong 54 matched taxon frame")
    return combined


def sample_stratum(group, selected, used_observers, target, quota, label):
    # Balanced across the same measured continuous-PC1 coordinate: the central
    # stratum is NOT assumed biologically intermediate or a third mode.
    hue = np.stack(group.palette_norm9.to_numpy()).astype(float)
    score, _ = mode_scores(np.sqrt(hue/hue.sum(axis=1)[:,None]))
    need(score is not None, "zero photo palette variation")
    position = float(np.quantile(score, target))
    g = group.copy().reset_index(drop=True)
    g["score"] = score
    g["rank_error"] = np.abs(score-position)
    g["_order"] = [hashed(label, int(v), "selection") for v in g.photo_id]
    g = g.sort_values(["rank_error","_order"],kind="stable")
    picked=[]
    for _, photo in g.iterrows():
        pid=int(photo.photo_id)
        oid=str(photo.observer_id)
        if pid in selected or (oid not in ("","nan","None") and oid in used_observers):
            continue
        picked.append((pid, label, float(photo.score)))
        selected.add(pid)
        if oid not in ("","nan","None"):
            used_observers.add(oid)
        if len(picked) == quota:
            break
    need(len(picked)==quota, f"insufficient distinct-photo/observer samples for {label}")
    return picked


def select_photo_rows(photo, classified, groups):
    flag = (photo.global_classifiable.fillna(False).astype(str).str.lower()
            .isin(["true","1","yes"]))
    need(int(flag.sum()) == CLASS_N, "classified photo denominator changed")
    rows, missing_summary = [], []
    for record in groups.itertuples(index=False):
        tid=int(record.inat_taxon_id)
        group = classified.loc[classified.inat_taxon_id == tid].copy()
        group = group.loc[(group.latitude.between(-90,90) &
                           group.longitude.between(-180,180) &
                           group.positional_accuracy_m.between(0,5000))]
        need(len(group) >= 40, f"taxon insufficient valid photos {tid}")
        selected=set()
        observers=set()
        strata=[]
        for quant, label in ((.12,"low_pc1"),(.50,"central_pc1"),(.88,"high_pc1")):
            strata.extend(sample_stratum(group,selected,observers,quant,3,label))
        unclassified=photo.loc[(photo.inat_taxon_id==tid)&(~flag)]
        unclassified=unclassified.sort_values("photo_id",kind="stable")
        unchosen=[]
        for obj in unclassified.itertuples(index=False):
            oid=str(obj.observer_id)
            if int(obj.photo_id) in selected:
                continue
            if oid in observers and oid not in ("","nan","None"):
                continue
            selected.add(int(obj.photo_id))
            if oid not in ("","nan","None"):
                observers.add(oid)
            unchosen.append((int(obj.photo_id),"unclassified_photo",None))
            if len(unchosen)==2:
                break
        strata.extend(unchosen)
        missing_summary.append({"taxon_id":tid,"review_arm":record.review_arm,
                                "n_original_unclassified":int(len(unclassified)),
                                "n_selected_unclassified":len(unchosen)})
        for pid,source_stratum,pc1 in strata:
            rows.append({"inat_taxon_id":tid,"species":record.species,
                         "pair_no":int(record.pair_no),
                         "review_arm":record.review_arm,
                         "photo_id":pid,"source_stratum":source_stratum,
                         "pc1_score":pc1})
    out=pd.DataFrame(rows)
    need(out.photo_id.is_unique and len(out) >= 54*9 and len(out) <= 54*11,
         "duplicate/invalid original reviewer photo selections")
    return out, pd.DataFrame(missing_summary)


def attach_original_meta(frame, photo_path):
    header = set(pd.read_csv(photo_path, nrows=0).columns)
    possible=("photo_url_large","photo_url","photo_license","photo_attribution",
              "observation_url","photo_attribution_text")
    optional=[v for v in possible if v in header]
    usecols=["photo_id"]+optional
    selected=set(frame.photo_id.astype(int))
    values=[]
    for chunk in pd.read_csv(photo_path, usecols=usecols, chunksize=30000, low_memory=False):
        part=chunk.loc[chunk.photo_id.isin(selected)]
        if len(part):
            values.append(part)
    found=pd.concat(values,ignore_index=True) if values else pd.DataFrame(columns=usecols)
    need(found.photo_id.is_unique,"source metadata photo ID duplication")
    out=frame.merge(found,on="photo_id",how="left",validate="one_to_one",indicator=True)
    need(out._merge.eq("both").all(),"review photos not in immutable original")
    return out.drop(columns="_merge"),optional


def blinded_packages(review, dest, download=False):
    # One anonymous species bundle per taxon; globally permuted between case/control.
    pairs=review[["inat_taxon_id","review_arm","pair_no"]].drop_duplicates().copy()
    pairs["bundle_id"]=["T"+hashed("taxon",int(x))[:10].upper() for x in pairs.inat_taxon_id]
    order=pairs.assign(_sort=pairs.inat_taxon_id.map(lambda x:hashed("bundleorder",int(x))))
    pairs=order.sort_values("_sort").drop(columns="_sort")
    rows=review.merge(pairs[["inat_taxon_id","bundle_id"]],on="inat_taxon_id",validate="many_to_one")
    rows["review_id"]=["I"+hashed("image",int(x))[:13].upper() for x in rows.photo_id]
    rows["_shuffle"]=[hashed("rowshuffle",str(x)) for x in rows.review_id]
    rows=rows.sort_values(["bundle_id","_shuffle"]).reset_index(drop=True)
    rows["review_order"]=np.arange(1,len(rows)+1)
    need(rows.review_id.is_unique and pairs.bundle_id.is_unique,"blind ID collision")
    sealed=dest/"pi_mapping_ONLY";sealed.mkdir(parents=True,exist_ok=True)
    rows.drop(columns="_shuffle").to_csv(sealed/"sealed_source_photo_mapping.csv",index=False)
    pairs.drop(columns="inat_taxon_id").to_csv(sealed/"taxon_pair_mapping.csv",index=False)
    public=dest/"blind_reviewer_ONLY";public.mkdir(parents=True,exist_ok=True)
    final=rows[["bundle_id","review_id","review_order"]].copy()
    for c in REVIEW_COLS:
        if c not in final:
            final[c]=""
    final=list(final.columns)
    for reviewer in (1,2):
        t=rows[["bundle_id","review_id"]].copy()
        t["_sort"]=[hashed("reviewer",reviewer,v) for v in t.review_id]
        t=t.sort_values("_sort").drop(columns="_sort").reset_index(drop=True)
        t["review_order"]=np.arange(1,len(t)+1)
        for key in REVIEW_COLS:
            if key not in t:
                t[key]=""
        # Reviewer sheets have no phenotype, taxon or source image ID.
        t[["bundle_id",*REVIEW_COLS]].to_csv(
            public/f"reviewer{reviewer}_blank_annotations.csv",index=False)
        t[["bundle_id","review_id","review_order"]].to_csv(
            public/f"reviewer{reviewer}_image_index.csv",index=False)
    return rows


LICENSED = {"cc0","cc-by","cc-by-sa","cc-by-nc","cc-by-nc-sa"}
SOURCE_BUCKET = "https://inaturalist-open-data.s3.amazonaws.com/photos"


def source_image_url(item):
    """Exact original URL, never an inferred photograph address."""
    for c in ("photo_url_large", "photo_url"):
        value=str(item.get(c, "")).strip()
        if value.startswith("https://"):
            host=(urlparse(value).hostname or "").lower()
            if (host.endswith(".inaturalist.org") or host=="inaturalist.org"
                or host.startswith("inaturalist-open-data.s3.")):
                return value
    return None


def photo_sources(item):
    licence=str(item.get("photo_license", "")).strip().lower()
    if licence not in LICENSED:
        return []
    pid=int(item["photo_id"])
    exact=source_image_url(item)
    candidates=([("authenticated_direct_source_url",exact)] if exact else [])
    # iNaturalist's documented licensed S3 photo ID/extension scheme.
    # Guessed paths are not evidence: count only a successful HTTP image decode
    # as a verified recovered original photo. Do not route unlicensed images.
    for ext in ("jpeg","jpg","png","webp"):
        url=f"{SOURCE_BUCKET}/{pid}/medium.{ext}"
        if url != exact:
            candidates.append(("source_id_s3_object_verified_by_http",url))
    return candidates


def try_image(item, output):
    options=photo_sources(item)
    if not options:
        return False,"ineligible_or_unknown_photo_license",None
    errors=[]
    for kind,url in options:
        for attempt in range(2):
            try:
                req=Request(url,headers={
                    "User-Agent":"zuizui0223-fcp-license-verified-blind-review/1.0",
                    "Accept":"image/*"})
                with urlopen(req,timeout=18) as r:
                    payload=r.read(7*1024*1024+1)
                if len(payload)>7*1024*1024 or len(payload)<1000:
                    raise ValueError("unexpected photo source size")
                with Image.open(io.BytesIO(payload)) as pic:
                    picture=ImageOps.exif_transpose(pic).convert("RGB")
                    picture.thumbnail((340,300),Image.Resampling.LANCZOS)
                    output.parent.mkdir(parents=True,exist_ok=True)
                    picture.save(output,format="JPEG",quality=82)
                return True,kind,url
            except Exception as e:
                errors.append(type(e).__name__)
                if getattr(e,"code",None) in (400,403,404,410):
                    break
                time.sleep(.3*(attempt+1))
    return False,"all_source_image_urls_unavailable:"+",".join(sorted(set(errors))),None


def make_sheets(rows, dest, fetch):
    blind=dest/"blind_reviewer_ONLY"
    files=blind/"images"
    files.mkdir(parents=True,exist_ok=True)
    if fetch:
        with ThreadPoolExecutor(max_workers=8) as pool:
            tasks={pool.submit(try_image,record,files/(record["review_id"]+".jpg")):record["review_id"]
                   for record in rows.to_dict(orient="records")}
            states={}
            for future in as_completed(tasks):
                states[tasks[future]]=future.result()
    else:
        states={r.review_id:(False,"image_download_not_requested",None) for r in rows.itertuples(index=False)}
    fields=[{"review_id":r.review_id,"rendered":bool(states[r.review_id][0]),
             "render_status":states[r.review_id][1],
             "verified_original_image_url":states[r.review_id][2]}
            for r in rows.itertuples(index=False)]
    pd.DataFrame(fields).to_csv(dest/"pi_mapping_ONLY"/"image_retrieval_audit.csv",index=False)
    for reviewer in (1,2):
        index=pd.read_csv(blind/f"reviewer{reviewer}_image_index.csv")
        per_page=20
        pages=blind/f"reviewer{reviewer}_contact_sheets"
        pages.mkdir(exist_ok=True)
        for start in range(0,len(index),per_page):
            subset=index.iloc[start:start+per_page]
            canvas=Image.new("RGB",(5*260,4*228),(245,245,245))
            draw=ImageDraw.Draw(canvas)
            for ix, row in enumerate(subset.itertuples(index=False)):
                x=(ix%5)*260;y=(ix//5)*228
                path=files/(row.review_id+".jpg")
                if path.exists():
                    with Image.open(path) as raw:
                        img=ImageOps.contain(raw,(248,192))
                        canvas.paste(img,(x+(260-img.width)//2,y+(192-img.height)//2))
                else:
                    draw.text((x+10,y+60),"IMAGE NOT\nAVAILABLE",fill=(80,80,80))
                draw.text((x+6,y+199),f"{row.review_order}: {row.review_id}",
                          fill=(25,25,25))
            canvas.save(pages/f"sheet_{start//per_page+1:02d}.jpg",quality=85)
    return fields


def synthetic():
    positives=pd.DataFrame({
        "inat_taxon_id":[1,2],"species":["A a","B b"],
        "mode_status":[MODE,MODE],
        "n_classified_original":[54,65],"pairs_50km":[130,120],
        "global_jsd":[.3,.4],"photographic_pc1_explained_variance":[.7,.6],
        "photographic_coarse_label_count":[2,3]})
    control=pd.DataFrame({
        "inat_taxon_id":[3,4],"species":["C c","D d"],
        "mode_status":[CONTROL,CONTROL],
        "n_classified_original":[55,64],"pairs_50km":[129,121],
        "global_jsd":[.31,.39],"photographic_pc1_explained_variance":[.69,.61],
        "photographic_coarse_label_count":[2,3]})
    score=np.full((60,9),1/3.)
    a,b=mode_scores(score)
    need(a is None,"truly identical-composition photo rows must be degenerate")
    hues=np.vstack([np.linspace(.1,.9,60),np.linspace(.9,.1,60)]).T
    nine=np.column_stack([hues,np.ones((60,7))*.05])
    pcs,_=mode_scores(np.sqrt(nine/nine.sum(axis=1)[:,None]))
    need(pcs is not None and np.isfinite(pcs).all(),"normal palette PC failed")
    mock=pd.DataFrame({"palette_norm9":list(nine/nine.sum(axis=1)[:,None]),
                       "photo_id":np.arange(60)+1000,
                       "observer_id":np.arange(60)+10})
    used=set();ob=set()
    collected=[]
    for q,label in ((.12,"low"),(.50,"middle"),(.88,"high")):
        collected+=sample_stratum(mock,used,ob,q,3,label)
    need(len(collected)==9 and len({x[0] for x in collected})==9,
         "blind-photo strata failed distinct photo selection")
    need(all(k not in REVIEW_COLS for k in
             ("species","mode_status","photo_id","latitude","longitude")),
         "private source field accidentally reviewer editable")
    print("FCP_27MODE_BLINDED_REVIEW_SYNTHETIC_PASS")


def main():
    parser=argparse.ArgumentParser()
    for v in ("photos","modes","continuous","out"):
        parser.add_argument("--"+v,type=Path)
    parser.add_argument("--self-test",action="store_true")
    parser.add_argument("--download-images",action="store_true")
    args=parser.parse_args()
    synthetic()
    if args.self_test:
        return
    need(all((args.photos,args.modes,args.continuous,args.out)),
         "original photos and both authenticated result artifacts required")
    source,classified=load_original(args.photos)
    need(len(source)==ALL_N and len(classified)==CLASS_N,"source drift")
    modes=pd.read_csv(args.modes)
    continuous=pd.read_csv(args.continuous)
    matched=match_species(modes,continuous,source)
    selected,missing=select_photo_rows(source,classified,matched)
    original,meta_cols=attach_original_meta(selected,args.photos)
    original["source_photo_page"]=[
        f"https://www.inaturalist.org/photos/{int(v)}" for v in original.photo_id
    ]
    args.out.mkdir(parents=True,exist_ok=True)
    rows=blinded_packages(original,args.out)
    states=make_sheets(rows,args.out,args.download_images)
    actual=int(sum(v["rendered"] for v in states))
    caseids=matched.loc[matched.review_arm=="candidate"].inat_taxon_id
    covariates=("n_classified_original","pairs_50km","global_jsd",
                "photographic_pc1_explained_variance")
    means={key:{
        "candidate_mean":float(matched.loc[matched.review_arm=="candidate",key].mean()),
        "control_mean":float(matched.loc[matched.review_arm=="matched_control",key].mean())}
        for key in covariates}
    manifest={
        "schema":"fcp_mode_blind_botanical_review_original_photos_v1",
        "source_photo_sha256":SOURCE_SHA,
        "source_original_photos":len(source),"classifiable_photos":len(classified),
        "previous_authenticated_original_mode_run":39923312111,
        "previous_authenticated_continuous_run":39918785955,
        "n_candidate_species":len(caseids),"n_one_to_one_matched_controls":PAIR_N,
        "total_selected_photo_rows":len(rows),
        "selected_petal_pc1_low_middle_high_each_species":3,
        "selected_unclassified_photo_rows":int(missing.n_selected_unclassified.sum()),
        "selected_reacquired_image_previews":actual,
        "direct_original_image_url_available":int(sum(source_image_url(r) is not None
                                                        for r in rows.to_dict("records"))),
        "documented_open_dataset_image_id_candidates":int(sum(
            bool(photo_sources(r)) for r in rows.to_dict("records"))),
        "optional_original_metadata_fields":meta_cols,
        "balance_summary":means,
        "image_retrieval_status":("photographic_previews_complete" if actual==len(rows)
                                  else "photographic_previews_incomplete"),
        "n_reviewers":2,
        "reviewers_are_blind_to_species_mode_date_site_observer_hue_score":True,
        "mapping_stored_separately_from_reviewer_directory":True,
        "public_repo_artifact_blind_privacy_not_cryptographic":True,
        "no_human_colour_or_botanical_judgements_yet":True,
        "reviewer_packets_are_blank_forms_not_validated_genetic_morphs":True,
        "note":"PI should distribute only blind_reviewer_ONLY; PI mapping is not an actual protected secret on a public CI system. Check source photo licences and attribution before sharing image sheets.",
    }
    (args.out/"review_manifest.json").write_text(json.dumps(manifest,indent=2,allow_nan=False)+"\n")
    print("FCP_27MODE_BLIND_REVIEW_READY",json.dumps(manifest,sort_keys=True))


if __name__=="__main__":
    main()
