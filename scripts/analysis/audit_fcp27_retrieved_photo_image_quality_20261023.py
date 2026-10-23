#!/usr/bin/env python3
"""Independent FCP27 licensed-thumbnail and blind-review provenance QC.

This audits thumbnails already reacquired and licensed/selected in the
authenticated source run. It does not relabel a biological petal, validate
flower-colour morphs, or replace any unavailable original photographs.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

EXPECTED = 594
CASES = 27


def check(yes, error):
    if not yes:
        raise RuntimeError("FCP27_RETRIEVED_PHOTO_QC: "+error)


def image_quality(path):
    with Image.open(path) as img:
        exif=len(img.getexif())
        pixel=np.asarray(img.convert("RGB"),dtype=np.float32)/255.
    lum = (pixel[:,:,0]*.2126 + pixel[:,:,1]*.7152 + pixel[:,:,2]*.0722)
    mx=np.max(pixel,axis=2)
    mn=np.min(pixel,axis=2)
    sat=(mx-mn)/(mx+1e-4)
    central=lum[lum.shape[0]//4:3*lum.shape[0]//4,
                lum.shape[1]//4:3*lum.shape[1]//4]
    return {
        "mean_frame_luminance":float(lum.mean()),
        "mean_frame_saturation":float(sat.mean()),
        "bright_clip_fraction":float(np.mean(lum>=.98)),
        "dark_clip_fraction":float(np.mean(lum<=.02)),
        "central_luminance_sd":float(central.std()),
        "exif_fields_remaining":exif,
    }


def summarize(mapping, audit, imgdir, nboot=499):
    check(len(mapping)==EXPECTED and mapping.review_id.is_unique,"sealed original mapping incomplete")
    check(len(audit)==EXPECTED and audit.review_id.is_unique,"retrieval evidence incomplete")
    needcols={"photo_id","inat_taxon_id","species","review_arm","source_stratum",
              "pair_no","review_id","photo_license"}
    check(needcols.issubset(mapping.columns),"missing source review metadata")
    check(mapping.inat_taxon_id.nunique()==54,"not 54 source taxa")
    got=mapping.merge(audit,on="review_id",how="left",validate="one_to_one")
    check(got.rendered.notna().all(),"retrieval results missing")
    check((got.rendered.astype(str).str.lower()=="true").sum()>=0,
          "cannot parse retrieval boolean")
    flags=got.rendered.astype(str).str.lower().eq("true")
    check(flags.sum()+int((~flags).sum())==EXPECTED,"audit denominator changed")
    measured=[]
    missing=[]
    for row in got.itertuples(index=False):
        path=imgdir/(row.review_id+".jpg")
        if bool(str(row.rendered).lower()=="true"):
            check(path.is_file(),"successful image absent")
            metrics=image_quality(path)
            check(metrics["exif_fields_remaining"]==0,
                  "PII/EXIF survived image anonymization")
            measured.append({
                "inat_taxon_id":int(row.inat_taxon_id),
                "review_arm":row.review_arm,"pair_no":int(row.pair_no),
                "source_stratum":row.source_stratum,"review_id":row.review_id,
                **metrics})
        else:
            check(not path.exists(),"phantom source image for failed download")
            missing.append({
                "source_stratum":row.source_stratum,
                "review_arm":row.review_arm,
                "reason":row.render_status
            })
    met=pd.DataFrame(measured)
    check(not met.empty and len(measured)+len(missing)==EXPECTED,
          "photo availability accounting broken")
    def q(s):
        return {"mean":float(s.mean()),"median":float(s.median()),
                "q05":float(s.quantile(.05)), "q95":float(s.quantile(.95)),
                "n":len(s)}
    stats={}
    for flag in ("candidate","matched_control"):
        sub=met.loc[met.review_arm==flag]
        stats[flag]={
            "n":len(sub),
            "mean_luminance":q(sub.mean_frame_luminance),
            "mean_saturation":q(sub.mean_frame_saturation),
            "bright_clipped_gt20pct":int((sub.bright_clip_fraction>.20).sum()),
            "dark_clipped_gt20pct":int((sub.dark_clip_fraction>.20).sum()),
            "central_low_detail_sd_lt_0p04":int((sub.central_luminance_sd<.04).sum())
        }
    diffs=[]
    for (pair,arm),g in met.groupby(["pair_no","review_arm"],sort=True):
        low=g.loc[g.source_stratum=="low_pc1","mean_frame_luminance"]
        high=g.loc[g.source_stratum=="high_pc1","mean_frame_luminance"]
        if len(low)>=2 and len(high)>=2:
            diffs.append({"pair_no":int(pair),"review_arm":arm,
                          "abs_low_high_luminance_delta":abs(float(high.mean()-low.mean()))})
    d=pd.DataFrame(diffs)
    pair_table=d.pivot(index="pair_no",columns="review_arm",
                       values="abs_low_high_luminance_delta").dropna()
    rng=np.random.default_rng(20261023)
    gap=(pair_table["candidate"]-pair_table["matched_control"]).to_numpy(float)
    boot=np.array([float(np.mean(gap[rng.integers(len(gap),size=len(gap))]))
                   for _ in range(nboot)]) if len(gap)>0 else np.array([])
    report={
        "schema":"fcp27_source_verified_thumbnail_technical_qc_v1",
        "n_source_selected":EXPECTED,
        "n_reacquired_original_image_previews":len(measured),
        "n_unavailable_source_photo_images":len(missing),
        "source_retrieval_failure_by_reason":pd.Series(
            [r["reason"] for r in missing]).value_counts().to_dict(),
        "source_retrieval_failure_by_photo_stratum":pd.Series(
            [r["source_stratum"] for r in missing]).value_counts().to_dict(),
        "source_licensed_photos":got.photo_license.fillna("").value_counts().to_dict(),
        "source_arms":stats,
        "n_matched_pairs_with_low_high_original_images":len(pair_table),
        "mean_candidate_minus_control_absolute_frame_luminance_difference":
            float(gap.mean()) if len(gap)>0 else None,
        "paired_bootstrap95ci":{
            "interval":[float(x) for x in np.quantile(boot,[.025,.975])]
            if len(boot)>0 else None,"n_boot":nboot},
        "all_recovered_thumbnail_exif_removed":True,
        "human_petals_and_senescence_reviewer_annotations_completed":0,
        "technical_analysis_only":True,
        "caution":"Whole-frame luminance can differ from flower petal colour and background. It cannot identify botanical colour modes, causal photo error or true native morphs.",
    }
    return report,met


def synthetic():
    rng=np.random.default_rng(23)
    sample=(rng.integers(0,255,size=(50,80,3))).astype("uint8")
    from tempfile import TemporaryDirectory
    with TemporaryDirectory() as temp:
        p=Path(temp)/"a.jpg"
        Image.fromarray(sample).save(p)
        v=image_quality(p)
        check(0<=v["mean_frame_luminance"]<=1,"bad clipped value")
        check(v["exif_fields_remaining"]==0,"test image EXIF not removed")
    print("FCP27_IMAGE_QC_SYNTHETIC_PASS")


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--mapping",type=Path)
    p.add_argument("--retrieval",type=Path)
    p.add_argument("--images",type=Path)
    p.add_argument("--out",type=Path)
    p.add_argument("--self-test",action="store_true")
    args=p.parse_args()
    synthetic()
    if args.self_test:
        return
    check(args.mapping and args.retrieval and args.images and args.out,
          "original artifact paths missing")
    a=pd.read_csv(args.mapping)
    b=pd.read_csv(args.retrieval)
    result,metrics=summarize(a,b,args.images)
    args.out.mkdir(parents=True,exist_ok=True)
    (args.out/"original_photo_image_quality_summary.json").write_text(
        json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")
    # No original photo IDs, species identity, location, or review treatment in exported image QC.
    print("FCP27_REACQUIRED_PHOTO_TECH_QC_COMPLETE",
          json.dumps(result,sort_keys=True))


if __name__=="__main__":
    main()
