#!/usr/bin/env python3
"""FCP mode-vs-space: blinded-photo visual pilot sensitivity, not confirmation.

The three flagged candidate taxa were selected *after viewing model-defined
low/central/high source photos* (NOT independent human adjudication). Removing
them is a qualitative stress test only; never call its nominal bootstrap CI
confirmatory or claim the remaining 24 have heritable morphs.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path

import numpy as np
import pandas as pd

MODE_SHA = "f799b8b0c0f6a250e5ff0c1ad15bf3150854f3099a6a9c03ad5a26711ba1e7df"
GEO_SHA = "6dfb1b8ef5b9bf60d8903b97f27257c4b2ddb26e45231b9870032a8b1b5cf528"
CASE = "recurrent_photographic_mode_candidate"
CONTROL = "not_resolved_as_multimodal"
FLAGGED = {
    "Argemone albiflora": "Photo-PC1 extremes show white open corollas; apparent density modes not a demonstrated petal-hue difference.",
    "Ochna serrulata": "Low PC1 includes red post-floral calyces/black fruits; high PC1 shows yellow corolla petals.",
    "Oxalis tenuifolia": "High PC1 includes red-pink unopened buds; central/low PC1 shows whitish open petals.",
}
REPS = 999
SEED = 20261023


def need(ok, detail):
    if not ok:
        raise RuntimeError("FCP27_POSTHOC_ORGAN_STAGE_AUDIT: "+detail)


def sha(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):
            h.update(b)
    return h.hexdigest()


def load(mode_path, geo_path, authenticated=True):
    if authenticated:
        need(sha(mode_path)==MODE_SHA and sha(geo_path)==GEO_SHA,
             "authenticated immutable original artifact bytes differ")
    m=pd.read_csv(mode_path,low_memory=False)
    g=pd.read_csv(geo_path,low_memory=False)
    need(len(m)==1761 and m.inat_taxon_id.is_unique,"mode taxa not exactly 1761")
    need(len(g)==1761 and g.taxon_id.is_unique,"geo taxa not exactly 1761")
    d=m.merge(g,left_on="inat_taxon_id",right_on="taxon_id",validate="one_to_one")
    need(len(d)==1761,"modal and spatial taxa mismatch")
    need(d.species_x.astype(str).eq(d.species_y.astype(str)).all(),
         "original species identity mismatch")
    need(np.max(np.abs(d.coarse_50km_depletion-d.coarse_50km_delta))<1e-10,
         "source coarse P2 mismatch after joining artifacts")
    need(d.mode_status.value_counts().to_dict()=={
        CONTROL:1691,CASE:27,"unstable_photo_mode_candidate":43},
        "original source mode statuses not fixed")
    return d.loc[d.mode_status.isin([CASE,CONTROL])].copy()


def design(frame, extra=False):
    z=np.column_stack([
        np.ones(len(frame)),
        frame.mode_status.eq(CASE).to_numpy(float),
        frame.global_jsd.to_numpy(float),
        np.log(frame.n_classified_original.to_numpy(float)),
        np.log1p(frame.pairs_50km.to_numpy(float))
    ])
    if extra:
        z=np.column_stack([
            z, frame.photographic_pc1_explained_variance.to_numpy(float),
            frame.photographic_coarse_label_count.to_numpy(float)
        ])
    y=frame.excess_jsd_50km_month.to_numpy(float)
    need(np.isfinite(z).all() and np.isfinite(y).all(),
         "source continuous estimand absent")
    return z,y


def fit(data, extra=False, reps=REPS):
    X,y=design(data,extra)
    need(np.linalg.matrix_rank(X)==X.shape[1],
         "covariate rank deficient")
    theta=float(np.linalg.lstsq(X,y,rcond=None)[0][1])
    blocks=[]
    # Exact genus-level cluster bootstrap with sufficient statistic reduction.
    for genus,subset in data.groupby("genus_x",sort=True):
        A,b=design(subset,extra)
        blocks.append((A.T@A,A.T@b))
    xx=np.stack([a for a,b in blocks])
    xy=np.stack([b for a,b in blocks])
    G=len(blocks)
    rng=np.random.default_rng(SEED+int(extra))
    samples=[]
    for _ in range(reps):
        cnt=np.bincount(rng.integers(G,size=G),minlength=G)
        A=np.einsum("g,gij->ij",cnt,xx)
        b=np.einsum("g,gi->i",cnt,xy)
        try:
            samples.append(float(np.linalg.solve(A,b)[1]))
        except np.linalg.LinAlgError:
            pass
    need(len(samples)>=int(.9*reps),"bootstrap rank failure")
    return {
        "beta_mode_minus_not_detected":theta,
        "genus_bootstrap_95ci":[float(x) for x in np.quantile(samples,[.025,.975])],
        "n_genus":G,"bootstrap_reps":len(samples),
        "adjustment":["species continuous JSD","log n classified photo",
                      "log(1+near pair count)"]
        +(["PC1 explained variance","coarse colour count"] if extra else [])
    }


def calculate(d,reps=REPS):
    flagged=d.loc[d.species_x.isin(FLAGGED)]
    need(len(flagged)==3 and flagged.mode_status.eq(CASE).all(),
         "inspected flags do not match original 27 mode cases")
    ans={}
    for name,removed in (("all_27_source_modes",[]),
                         ("excluding_3_pilot_flags",list(FLAGGED))):
        t=d.loc[~d.species_x.isin(removed)].copy()
        n_case=int(t.mode_status.eq(CASE).sum())
        n_control=int(t.mode_status.eq(CONTROL).sum())
        means=t.groupby("mode_status").excess_jsd_50km_month.mean()
        ans[name]={
            "n_mode_candidate":n_case,
            "n_controls":n_control,
            "mean_month_adjusted_50km_continuous_geo_by_group":{
                k:float(v) for k,v in means.items()},
            "unadjusted_mode_minus_control":
                float(means[CASE]-means[CONTROL]),
            "base":fit(t,False,reps),
            "extra_adjustment":fit(t,True,reps),
        }
    return {
        "schema":"fcp27_unblinded_ai_visual_pilot_organ_stage_source_conditional_v1",
        "source_taxon_frame":1761,
        "source_original_photo_run":39923312111,
        "source_continuous_run":39918785955,
        "visual_flag_examples":FLAGGED,
        "visually_reviewed_flags_NOT_independent_human":True,
        "uninspected_24_not_assumed_valid":True,
        "flag_selection_after_outcome_unregistered":True,
        "photo_strata_source":"original low/central/high PC1 review plates",
        "estimates":ans,
        "not_genetic_morphs_or_adaptation":True,
        "confidence_intervals_exploratory_not_post_selection_valid":True,
    }


def self_test():
    import tempfile
    example=pd.DataFrame({
        "mode_status":[CASE]*10+[CONTROL]*40,
        "global_jsd":np.linspace(.1,.6,50),
        "n_classified_original":np.tile([44,50,63,79,91],10),
        "pairs_50km":np.arange(50)+30,
        "photographic_pc1_explained_variance":np.tile([.4,.5,.65,.74,.9],10),
        "photographic_coarse_label_count":np.tile([1,2,3,4,2],10),
        "excess_jsd_50km_month":np.arange(50)/1000.,
        "genus_x":[f"genus{j//2}" for j in range(50)]
    })
    x,y=design(example)
    need(x.shape==(50,5) and len(y)==50,
         "synthetic regressor dimensions wrong")
    need(not np.isnan(x).any(),"null synthetic features")
    print("FCP27_ORGAN_STAGE_SENSITIVITY_SYNTHETIC_PASS")


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--modes",type=Path)
    p.add_argument("--geography",type=Path)
    p.add_argument("--out",type=Path)
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    self_test()
    if a.self_test:
        return
    need(a.modes and a.geography and a.out,"missing authenticated input paths")
    d=load(a.modes,a.geography)
    result=calculate(d)
    a.out.mkdir(parents=True,exist_ok=True)
    (a.out/"source_conditional_organ_stage_sensitivity.json").write_text(
        json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("FCP27_POSTHOC_ORGAN_STAGE_SENSITIVITY_COMPLETE",
          json.dumps(result,sort_keys=True,allow_nan=False))


if __name__=="__main__":
    main()
