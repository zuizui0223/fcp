#!/usr/bin/env python3
"""Split full 36-environment heldout multiclass prediction gain by original 4 photo hues.

Original FCP 42111 taxa / 100543 photo IDs, no hue reclassification.
Same genus-conditioned species-mean 5fold photo centroid-cell CV as parent
exact 128-block Shapley source. Brier class increments sum to full environment
improvement; NOT additive biochemical colour costs or selection coefficients.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from sklearn.model_selection import GroupKFold
from fcp_audit_36_preflight_20261011 import run as source,BLOCKS,CLASSES
from fcp_all36_global_region_cv_20261011 import aggregate_species,fit_predict
GEO=("abs_lat","lon_sin","lon_cos")
SCHEMA="fcp_original100543_joint36_photocolour_hue_brier_contribution_v1"
N_BOOT=499

def evaluate(raw,mode,region="world",strict=False):
    d=raw.copy()
    if strict:
        d=d[(d.flower_effective_pixels>=23826.)&
            (d.nuisance_pixel_fraction<=0.1375717742195135)]
    if region=="east":d=d[d.east]
    if region=="elsewhere":d=d[~d.east]
    sp=aggregate_species(d)
    feat=GEO+tuple(k for block,fs in BLOCKS.items()
                   if mode=="soil" or block!="soil" for k in fs)
    y=sp[["hue_"+x for x in CLASSES]].to_numpy(float)
    gain=np.full((len(sp),4),np.nan)
    folds=[]
    for k,(tr,te) in enumerate(GroupKFold(n_splits=5).split(sp,groups=sp.cell)):
        train,test=sp.iloc[tr],sp.iloc[te]
        good=train.genus.value_counts()
        choice=np.flatnonzero(test.genus.isin(good[good>=3].index).to_numpy())
        target=test.iloc[choice]
        if not len(target):continue
        assert not set(train.inat_taxon_id)&set(target.inat_taxon_id)
        pgeo=fit_predict(train,target,GEO)
        penv=fit_predict(train,target,feat)
        truth=y[te[choice]]
        gain[te[choice]]=(pgeo-truth)**2-(penv-truth)**2
        folds.append({"fold":k,"n_geoheldout":len(te),"n_eligible":len(choice)})
    keep=np.isfinite(gain).all(axis=1)
    if keep.sum()<200:raise RuntimeError("Insufficient heldout photo-taxonomic genus coverage")
    g=gain[keep]
    full=np.column_stack((g,g.sum(axis=1)))
    genus=sp.loc[keep,"genus"].to_numpy(str)
    _,grp=np.unique(genus,return_inverse=True)
    ng=int(grp.max())+1
    n=np.bincount(grp,minlength=ng)
    means=np.column_stack([np.bincount(grp,weights=full[:,j])/n for j in range(5)])
    rng=np.random.default_rng(2026101113+len(region)+(0 if mode=="soil" else 100)+(200 if strict else 0))
    bs=[]
    for start in range(0,N_BOOT,50):
        ids=rng.integers(0,ng,size=(min(50,N_BOOT-start),ng))
        bs.append(means[ids].mean(axis=1))
    ci=np.quantile(np.concatenate(bs),[.025,.975],axis=0).T
    keys=[*CLASSES,"total_multiclass_Brier"]
    x={key:{"genus_equal_gain":float(means[:,j].mean()),
            "genus_cluster_95CI":ci[j].tolist(),
            "species_equal_gain":float(full[:,j].mean())} for j,key in enumerate(keys)}
    if not np.isclose(sum(x[k]["genus_equal_gain"] for k in CLASSES),
                      x["total_multiclass_Brier"]["genus_equal_gain"],atol=1e-12):
        raise RuntimeError("Original four flower hue error gain does not reassemble to full Brier")
    return {"n_source_test_species":int(keep.sum()),"n_test_genera":ng,
            "n_source_photo_rows":len(d),"source_folds":folds,
            "hue_Brier_contributions":x,"four_class_gain_reassembles_full":True}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--mode",choices=("soil","climate","both"),default="both")
    p.add_argument("--regions",choices=("world","east","elsewhere","all"),default="all")
    p.add_argument("--quality",action="store_true")
    p.add_argument("--outdir",required=True,type=Path)
    a=p.parse_args()
    u,denom,masks=source()
    modes=["soil","climate"] if a.mode=="both" else [a.mode]
    regions=["world","east","elsewhere"] if a.regions=="all" else [a.regions]
    cases={}
    for mode in modes:
        for region in regions:
            key=mode+"_"+region+("_quality" if a.quality else "")
            cases[key]=evaluate(u.loc[masks[mode]],mode,region,a.quality)
    result={"schema":SCHEMA,"original_photo_ids":100543,
            "original_classified_photo_ids":45953,"source_complete_photo_counts":denom,
            "source_colour_measurements_untouched":True,"source_nongenetic_photo_fitness_information":False,
            "cases":cases,
            "nonclaims":["Original photo white labels are exposure-sensitive and not pigment molecular identities.",
                "Class-specific Brier gains sum to full multiclass prediction, not separate causal fitness or colour transition components.",
                "Original photo centroid-cell heldout, not complete phylogenetic/spatial causal control."]}
    a.outdir.mkdir(exist_ok=True,parents=True)
    (a.outdir/"result.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({k:{"species":v["n_source_test_species"],
                        "genus_equal_Brier":{name:round(z["genus_equal_gain"],9)
                        for name,z in v["hue_Brier_contributions"].items()}}
                     for k,v in cases.items()},indent=2))
if __name__=="__main__":main()
