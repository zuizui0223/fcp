#!/usr/bin/env python3
"""Cross-fit source photographs to audit apparent white-state thermal modulation.

Original FCP 447,500 measured photos. All observed photos from a photographer
within a species stay in ONE of two deterministic, outcome-blind folds.
Group baseline white-classified photo abundance on fold A, estimate thermal /
white 9-colour photo regression on entirely different photographers' fold B.
Swap fold roles and re-run, plus one-photo-per-observer outcome sensitivity.

No genotype, independently reviewed petals, pigment cost or selection measured.
The photographic baseline and out-of-fold response are independent OBSERVERS,
not independent cameras/platforms; the source classification is incomplete.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"/"analysis"))
from fcp_full4475_hypothesis_census_20261021 import (
    SOURCE_SHA,BIO_ARCHIVE_SHA,TOTAL_ROWS,TOTAL_TAXA,CLASSIFIABLE_ROWS,
    load_original,sample_climate,all_coords_ok,PHOTO_LABELS
)
from fcp_all4475_seven_abiotic_axes_20261021 import add_srad,SRAD_ZIP_SHA,AXES
from fcp4475_nine_colour_seven_climate_covariance_20261024 import design_geo_month
from fcp4475_shared_thermal_moisture_photo_hues_20261024 import (
    compose_macros,DOMAINS
)

SEED=20261024
NULL_REPS=199
GENUS_REPS=999
MIN_BASELINE_CLASS=10
MIN_TEST_CLASS=18
MIN_BASELINE_OBSERVERS=8
MIN_TEST_OBSERVERS=10
MIN_TEST_CELLS=3
RESID_FRACTION=.05
COLS=("none","rare","mixed","nearly_all")
COHORTS=("discovery","validation","prospective","wave4","new_U100")
PANEL_ORDER=("forward","reverse")
TEST_MODE=("all_heldout","one_photo_per_observer")
PRIMARY_THRESHOLD=3
STRICT_THRESHOLD=5


def ensure(ok,detail):
    if not ok:raise RuntimeError("FCP_CROSSFIT_THERMAL_PHOTO: "+detail)


def sha1_text(s):
    return int.from_bytes(hashlib.sha256(str(s).encode()).digest()[:8],"big")


def observer_folds(source_taxon):
    """Photo-COLOUR-BLIND fold allocation, all source photos not only classified."""
    obs=source_taxon.observer_id.dropna().unique()
    tokens=sorted(obs,key=lambda u:(sha1_text(str(int(source_taxon.inat_taxon_id.iloc[0]))+
                             "|"+str(u)+"|"+str(SEED)),str(u)))
    return {token:i%2 for i,token in enumerate(tokens)}


def category(nwhite,nnonwhite,threshold=PRIMARY_THRESHOLD):
    """Based only on baseline fold photos, NEVER the heldout test photos."""
    if nwhite==0:return "none"
    if nwhite<threshold:return "rare"
    if nnonwhite>=threshold:return "mixed"
    return "nearly_all"


def withheld_slope(photo_subset,taxon,fold,thin=False):
    """One species' heldout observational source white/thermal partial slope.

    Colour/photo month expected covariance subtracted EXACTLY. Nulls shuffle
    white ROI within same calendar month; no cross-month/source-venue shuffle.
    """
    value={
      "heldout_n_classified":len(photo_subset),
      "heldout_n_complete":0,"heldout_n_observers":0,"heldout_n_cells":0,
      "thermal_unique_var_frac":None,"heldout_BIO5_q90_q10_C":None,
      "status":"unavailable",
    }
    if thin:
        photo_subset=(photo_subset.sort_values(["observer_id","photo_id"],kind="stable")
                     .drop_duplicates("observer_id",keep="first"))
    g=photo_subset
    geo=all_coords_ok(g)
    valid=geo&g.date.notna().to_numpy()&np.isfinite(g[list(AXES)].to_numpy(float)).all(axis=1)
    g=g.loc[valid].reset_index(drop=True)
    n=len(g)
    nob=g.observer_id.nunique()
    value["heldout_n_complete"]=int(n)
    value["heldout_n_observers"]=int(nob)
    if n<MIN_TEST_CLASS or nob<MIN_TEST_OBSERVERS:
        value["status"]="insufficient_heldout_photo_or_observers"
        return value,None,None
    cells=np.column_stack([
        np.floor((90-g.latitude.to_numpy(float))*6),
        np.floor((g.longitude.to_numpy(float)+180)*6),
    ]).astype(int)
    ncell=len(np.unique(cells,axis=0))
    value["heldout_n_cells"]=int(ncell)
    if ncell<MIN_TEST_CELLS:
        value["status"]="insufficient_heldout_raster_cells"
        return value,None,None
    env=g[list(AXES)].to_numpy(float)
    q5=np.quantile(env[:,AXES.index("BIO5")],[.1,.9])
    value["heldout_BIO5_q90_q10_C"]=float(q5[1]-q5[0])
    means=env.mean(axis=0)
    sd=env.std(axis=0)
    temp_ok=sd>np.maximum(1e-7,2e-10*np.maximum(np.abs(means),1.))
    if not temp_ok[0] or not temp_ok[2]:
        value["status"]="no_variable_temperature_in_heldout"
        return value,None,None
    z=(env-means)/np.where(temp_ok,sd,1.)
    mac=compose_macros(z)
    design,month=design_geo_month(g)
    rank=np.linalg.matrix_rank(design,tol=1e-9)
    value["heldout_geo_month_rank"]=int(rank)
    if n-rank<8:
        value["status"]="insufficient_degrees_after_geo_month"
        return value,None,None
    residual=mac-design@np.linalg.lstsq(design,mac,rcond=1e-8)[0]
    warm=residual[:,0]
    controls=residual[:,1:]
    x=warm-controls@np.linalg.lstsq(controls,warm,rcond=1e-9)[0]
    v=float(np.mean(x*x))
    raw=float(np.mean(mac[:,0]**2))
    retained=v/max(raw,1e-12)
    value["thermal_unique_var_frac"]=float(retained)
    if raw<1e-9 or retained<RESID_FRACTION or not np.isfinite(v):
        value["status"]="collinear_temperature_after_heldout_controls"
        return value,None,None
    palette=np.stack(g.palette_norm9.to_numpy()).astype(float)
    ensure(palette.shape==(n,9) and np.allclose(palette.sum(axis=1),1,atol=1e-10,rtol=0),
           "source original heldout nine-hue vector lost")
    y=palette[:,0]
    xc=x.copy()
    groups=[]
    for m in np.unique(month):
        idx=np.flatnonzero(month==m)
        xc[idx]-=xc[idx].mean()
        if len(idx)>1:groups.append(idx)
    # Identical to observed covariance less the exact calendar-month-matched
    # climate permutation expectation, divided by the unpermuted climate variance.
    slope=float(np.dot(xc,y)/(n*v))
    rng=np.random.default_rng(SEED+int(taxon)*31+int(fold)*71+int(thin)*101)
    null=np.zeros(NULL_REPS,float)
    for idx in groups:
        orders=np.argsort(rng.random((NULL_REPS,len(idx))),axis=1)
        null+=(y[idx][orders]@xc[idx])/(n*v)
    value["status"]="crossfitted_heldout_thermal_white_estimable"
    value["heldout_thermal_white_slope"]=slope
    value["heldout_null_mean"]=float(null.mean())
    ensure(np.isfinite(slope) and np.isfinite(null).all(),
           "heldout source/photo temperature slope not finite")
    return value,slope,null


def photo_split(source_raw,source_classified):
    """Return all 4,475 source taxa x 2 directions, observer-disjoint photos."""
    out=[]
    allobs=source_raw.groupby("inat_taxon_id",sort=True)
    n_class_source=0
    for count,(taxon,cl) in enumerate(source_classified.groupby("inat_taxon_id",sort=True)):
        raw=allobs.get_group(taxon)
        obsfold=observer_folds(raw)
        s=cl.loc[cl.observer_id.notna()].copy()
        s["crossfit_fold"]=s.observer_id.map(obsfold)
        ensure(s.crossfit_fold.isin([0,1]).all(),"crossfit observer photo overlap/missing fold")
        # A heldout photographer must not appear in the opposing baseline.
        both=(set(s.loc[s.crossfit_fold==0,"observer_id"]) &
              set(s.loc[s.crossfit_fold==1,"observer_id"]))
        ensure(not both,"heldout baseline photographer overlap")
        n_class_source+=len(cl)
        for direction in (0,1):
            base=s.loc[s.crossfit_fold==direction]
            test=s.loc[s.crossfit_fold==1-direction]
            nwhite=int((base.morph=="white").sum())
            nnw=int(len(base)-nwhite)
            row={
              "inat_taxon_id":int(taxon),"species":str(cl.species.iloc[0]),
              "genus":str(cl.species.iloc[0]).split()[0],
              "source_cohort":str(cl.fcp_source_cohort.iloc[0]),
              "direction":"forward" if direction==0 else "reverse",
              "baseline_fold":direction,
              "test_fold":1-direction,
              "n_classified_original":len(cl),
              "n_original_photo":100,
              "n_observer_missing_in_source":int(raw.observer_id.isna().sum()),
              "baseline_n_classified":len(base),
              "baseline_n_observers":int(base.observer_id.nunique()),
              "baseline_n_white":nwhite,
              "baseline_n_nonwhite":nnw,
              "baseline_white_fraction":float(
                  np.mean([float(v[0]) for v in base.palette_norm9])) if len(base) else None,
              "baseline_group3":category(nwhite,nnw,PRIMARY_THRESHOLD),
              "baseline_group5":category(nwhite,nnw,STRICT_THRESHOLD),
              "baseline_is_adequately_sampled":bool(
                  len(base)>=MIN_BASELINE_CLASS and
                  base.observer_id.nunique()>=MIN_BASELINE_OBSERVERS),
              "test_n_classified":len(test),
              "test_n_observers":int(test.observer_id.nunique())
            }
            for thin in (False,True):
                record=row.copy()
                record["test_photo_mode"]="all_heldout" if not thin else "one_photo_per_observer"
                z,estimate,null=withheld_slope(test,taxon,direction,thin)
                record.update(z)
                # Label source-original taxon row even if baseline unsupported.
                if not row["baseline_is_adequately_sampled"]:
                    record["status"]="baseline_photo_opportunity_insufficient"
                    record["heldout_thermal_white_slope"]=None
                    null=None
                record["_null"]=null
                out.append(record)
        if (count+1)%1000==0:
            print("FCP_CROSSFIT_SOURCE_TAXA",count+1,flush=True)
    ensure(n_class_source==CLASSIFIABLE_ROWS,"original nine-photo-classified denominator changed")
    d=pd.DataFrame(out)
    ensure(len(d)==TOTAL_TAXA*2*2 and
           d.groupby(["direction","test_photo_mode"]).size().eq(TOTAL_TAXA).all(),
           "original 4,475 source taxa not retained in four crossfit panels")
    return d


def group_bootstrap(frame,group_col,seed=SEED,rep=GENUS_REPS):
    d=frame.loc[frame.heldout_thermal_white_slope.notna()].copy()
    categories=COLS
    data={}
    rng=np.random.default_rng(seed)
    for cat in categories:
        x=d.loc[d[group_col]==cat]
        y=x.heldout_thermal_white_slope.to_numpy(float)
        codes,genus=pd.factorize(x.genus,sort=True)
        G=len(genus)
        if len(y)<20 or G<15:
            data[cat]={"status":"not_estimable","n_taxa":len(y),"n_genera":G}
            continue
        sums=np.zeros(G,float)
        counts=np.bincount(codes,minlength=G)
        np.add.at(sums,codes,y)
        draws=rng.integers(G,size=(rep,G))
        boots=sums[draws].sum(axis=1)/counts[draws].sum(axis=1)
        data[cat]={
            "status":"source_observer_disjoint_crossfit",
            "n_taxa":int(len(y)),"n_genera":int(G),
            "mean":float(y.mean()),"median":float(np.median(y)),
            "fraction_positive":float(np.mean(y>0)),
            "genus_bootstrap_marginal95":[float(v) for v in np.quantile(boots,[.025,.975])]
        }
    return data


def regression(frame,group_col,seed=SEED):
    """Regress heldout source-photo slopes on INDEPENDENT observer baseline groups."""
    d=frame.loc[frame.heldout_thermal_white_slope.notna()].copy()
    covars=[
        "baseline_n_classified","heldout_n_complete","heldout_n_cells",
        "heldout_BIO5_q90_q10_C","thermal_unique_var_frac"
    ]
    d=d.loc[np.isfinite(d[covars].to_numpy(float)).all(axis=1)].copy()
    if len(d)<100:
        return {"status":"insufficient_independent_photo_opportunity","n_taxa":len(d)}
    labels=pd.get_dummies(d[group_col],dtype=float).reindex(columns=list(COLS)[1:],fill_value=0.)
    cohort=pd.get_dummies(d.source_cohort,dtype=float).reindex(
        columns=list(COHORTS)[1:],fill_value=0.)
    extra=d[covars].to_numpy(float)
    extra[:,2]=np.log1p(extra[:,2])
    extra[:,3]=np.log1p(extra[:,3])
    sd=extra.std(axis=0)
    extra=(extra-extra.mean(axis=0))/np.maximum(sd,1e-12)
    X=np.column_stack((np.ones(len(d)),labels.to_numpy(float),extra,cohort.to_numpy(float)))
    names=["intercept",*["crossfit_"+g+"_vs_none" for g in COLS[1:]],
           *["z_"+c for c in covars],
           *["source_cohort_"+k for k in COHORTS[1:]]]
    p=X.shape[1]
    if np.linalg.matrix_rank(X)!=p:
        return {"status":"photographic_group_design_rank_deficient","n_taxa":len(d)}
    y=d.heldout_thermal_white_slope.to_numpy(float)
    estimates=np.linalg.lstsq(X,y,rcond=1e-8)[0]
    codes,genus=pd.factorize(d.genus,sort=True)
    G=len(genus)
    cross=np.zeros((G,p,p),float)
    rhs=np.zeros((G,p),float)
    np.add.at(cross,codes,np.einsum("ni,nj->nij",X,X))
    np.add.at(rhs,codes,X*y[:,None])
    picks=np.random.default_rng(seed).integers(G,size=(GENUS_REPS,G))
    weights=np.zeros((GENUS_REPS,G),float)
    np.add.at(weights,(np.repeat(np.arange(GENUS_REPS),G),picks.ravel()),1)
    gram=np.einsum("bg,gij->bij",weights,cross)
    target=weights@rhs
    ensure(gram.shape==(GENUS_REPS,p,p) and target.shape==(GENUS_REPS,p),
           "genus bootstrap OLS shape contract broken")
    bs=[]
    for a,b in zip(gram,target):
        try:bs.append(np.linalg.solve(a,b))
        except np.linalg.LinAlgError:pass
    ensure(len(bs)>.90*GENUS_REPS,
           "crossfit group regression unstable under genus resampling")
    boots=np.stack(bs)
    cis=np.quantile(boots,[.025,.975],axis=0).T
    # Group response null computed from internal source-month matched shuffles,
    # holding SOURCE fold-baseline groups and photo opportunity covariates fixed.
    nullY=np.stack(d["_null"].to_numpy()).astype(float)
    ensure(nullY.shape==(len(d),NULL_REPS),"crossfit source null/photo response invalid")
    nullCoef=np.linalg.lstsq(X,nullY,rcond=1e-8)[0]
    results={}
    for i,name in enumerate(names):
        observed=float(estimates[i])
        q=np.quantile(boots[:,i],[.025,.975])
        null=nullCoef[i]
        p2=float((1+np.count_nonzero(np.abs(null)>=abs(observed)))/(NULL_REPS+1))
        results[name]={
            "coefficient":observed,
            "genus_bootstrap95":[float(v) for v in q],
            "two_sided_month_matched_original_source_null_p_199":p2,
            "genus95_excludes_zero":bool(q[0]>0 or q[1]<0)
        }
    return {
        "status":"observer_disjoint_source_photo_regression",
        "n_taxa":int(len(d)),"n_genera":int(G),
        "baseline_category":"no independent baseline photo labelled white",
        "n_monthly_photo_null_replicates":NULL_REPS,
        "coef":results,
        "genus_bootstrap_draws":len(bs),
        "scope":"noncausal original nine-hue PHOTO and thermal climate relation, same source platform"
    }


def one_panel(d,direction,mode):
    sub=d.loc[(d.direction==direction)&(d.test_photo_mode==mode)].copy()
    ensure(len(sub)==TOTAL_TAXA and sub.inat_taxon_id.is_unique,
           "missing original taxa from crossfit panel")
    out={
        "n_original_source_species":TOTAL_TAXA,
        "n_independent_baseline_supported":int(sub.baseline_is_adequately_sampled.sum()),
        "n_heldout_thermal_estimable":int(sub.heldout_thermal_white_slope.notna().sum()),
        "status_counts":sub.status.value_counts().to_dict(),
        "n_valid_both3":int(((sub.baseline_group3=="mixed")&sub.heldout_thermal_white_slope.notna()).sum()),
        "three_image_baseline_group":{
            "results":group_bootstrap(sub,"baseline_group3",SEED+int(direction=="reverse")*83+int(mode=="one_photo_per_observer")*31),
            "regression":regression(sub,"baseline_group3",SEED+int(direction=="reverse")*41+int(mode=="one_photo_per_observer")*23)
        },
        "five_image_baseline_group":{
            "results":group_bootstrap(sub,"baseline_group5",SEED+103+int(direction=="reverse")*83+int(mode=="one_photo_per_observer")*31),
            "regression":regression(sub,"baseline_group5",SEED+103+int(direction=="reverse")*41+int(mode=="one_photo_per_observer")*23)
        }
    }
    # Stratified baseline photo-white fractions, computed on an independent set.
    seen=sub.loc[sub.heldout_thermal_white_slope.notna()&
                 sub.baseline_white_fraction.notna()].copy()
    if len(seen)>=100:
        buckets=pd.qcut(seen.baseline_white_fraction,q=4,labels=False,duplicates="drop")
        labels=["Q"+str(int(v)+1) for v in buckets]
        seen["_q"]=labels
        out["independent_baseline_white_fraction_quartiles"]={
            "bounds":[float(x) for x in np.unique(np.quantile(seen.baseline_white_fraction,[0,.25,.5,.75,1]))],
            "n_taxa":len(seen),
            "means":[{"group":q,"n_taxa":len(y),
                      "mean_heldout_white_warmth_slope":float(y.heldout_thermal_white_slope.mean())}
                     for q,y in seen.groupby("_q",sort=True)]
        }
    return out


def synthetic():
    fake=pd.DataFrame({
        "inat_taxon_id":[1001]*4,
        "observer_id":[1,1,2,3],
        "photo_id":[1,2,3,4]
    })
    a=observer_folds(fake)
    ensure(a[1] in (0,1) and a[2] in (0,1) and a[3] in (0,1),
           "whole-observer photo split not deterministic")
    ensure(a==observer_folds(fake),"observer hash split nonrepeatable")
    ensure(category(0,4)=="none" and category(2,5)=="rare"
           and category(3,6)=="mixed" and category(4,2)=="nearly_all",
           "crossfit photo-white four categories invalid")
    # Known PHOTO climate effect independent of observation month:
    x=np.array([-2.,-1.,1.,2.,-2.,-1.,1.,2.])
    y=np.array([.1,.2,.4,.5,.1,.2,.4,.5])
    month=np.array([5,5,5,5,6,6,6,6])
    centered=x.copy()
    for m in np.unique(month):
        idx=month==m
        centered[idx]-=centered[idx].mean()
    exact=float(centered@y/len(y))
    expected=0
    for m in np.unique(month):
        idx=(month==m)
        expected+=(x[idx].sum()*y[idx].mean())/len(y)
    ensure(np.isclose(exact,x@y/len(y)-expected,atol=1e-12) and exact>0,
           "photo-calendar-month-conditional covariance identity broken")
    # NumPy genus regression batch/vector support fixture.
    M=np.tile(np.eye(12)[None],(7,1,1))
    b=np.ones((7,12))
    ans=np.linalg.solve(M,b[...,None])[...,0]
    ensure(ans.shape==(7,12),"genus regression solve batch dimensions broken")
    print("FCP_CROSSFITTED_447500_PHOTO_OBSERVER_INDEPENDENCE_SELFTEST_PASS")


def main():
    p=argparse.ArgumentParser()
    for arg in ("photos","worldclim-bio","worldclim-srad","out"):
        p.add_argument("--"+arg,type=Path)
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    synthetic()
    if a.self_test:return
    ensure(a.photos and a.worldclim_bio and a.worldclim_srad and a.out,
           "three original SHA-pinned measured photo and climate inputs required")
    full,d=load_original(a.photos)
    d=sample_climate(d,a.worldclim_bio)
    d=add_srad(d,a.worldclim_srad)
    ensure(len(full)==TOTAL_ROWS and len(d)==CLASSIFIABLE_ROWS
           and d.inat_taxon_id.nunique()==TOTAL_TAXA,
           "source 447500 original photo/4475 species/227107 classifiable count drift")
    rows=photo_split(full,d)
    a.out.mkdir(parents=True,exist_ok=True)
    to_save=rows.drop(columns=["_null"])
    to_save.to_csv(a.out/"fcp4475_source_all_taxa_twofold_observer_crossfit_thermal_white.csv",
                   index=False,float_format="%.12g")
    panels={}
    for direction in PANEL_ORDER:
        for mode in TEST_MODE:
            panel=direction+"__"+mode
            panels[panel]=one_panel(rows,direction,mode)
            print("FCP_CROSSFIT_PANEL",panel,panels[panel]["n_heldout_thermal_estimable"],
                  flush=True)
    # Note: no shared photo observer IDs or per-photo coordinates in derived
    # artifact; source test identifiers remain at taxon level.
    summary={
      "schema":"fcp447500_observer_disjoint_crossfit_white_baseline_thermal_slope_v1",
      "date_jst":"2026-10-24",
      "source_original_photo_sha256":SOURCE_SHA,
      "original_BIO_worldclim_sha256":BIO_ARCHIVE_SHA,
      "original_SRAD_worldclim_sha256":SRAD_ZIP_SHA,
      "source_original_photos":TOTAL_ROWS,
      "source_taxa":TOTAL_TAXA,
      "source_classified_nine_photo_hue":CLASSIFIABLE_ROWS,
      "source_photo_colour_unclassified":TOTAL_ROWS-CLASSIFIABLE_ROWS,
      "photo_fold_assignment":"SHA256 taxon_id|observer_id|seed ranking on ALL ORIGINAL photos before photo-colour or climate eligibility filters, alternating original observers",
      "baseline_and_heldout_photo_observers_strictly_disjoint":True,
      "fold_direction":"forward baseline0 -> heldout1; reverse baseline1 -> heldout0",
      "baseline_min_classified_photos":MIN_BASELINE_CLASS,
      "baseline_min_independent_observers":MIN_BASELINE_OBSERVERS,
      "heldout_min_climate_photo":MIN_TEST_CLASS,
      "heldout_min_observers":MIN_TEST_OBSERVERS,
      "heldout_min_distinct_source_worldclim_cells":MIN_TEST_CELLS,
      "heldout_min_unabsorbed_temperature_residual_fraction":RESID_FRACTION,
      "heldout_test_climate_model":"same five fixed source climate domains as preceding original analysis: shared BIO1/BIO5 thermal z, BIO12/BIO14 shared precip, BIO4, BIO15 and solar; heldout-only climatological normalization and broad XYZ+calendar month climate adjustment",
      "heldout_original_photo_temperature_slope":"nine-palette white versus heldout uniquely residualized standardized macrotemperature after exact month-conditional null; observational conditional regression, not temperature causality",
      "original_photo_group_thresholds":[PRIMARY_THRESHOLD,STRICT_THRESHOLD],
      "photographed_white_baseline_not_a_genetic_morph":True,
      "both_photo_halves_same_platform_not_independent_colour_assay":True,
      "photo_background_and_flower_organs_not_independently_validated":True,
      "folds_allow_original_same_species_multiple_photos_per_observer_only_inside_one_fold":True,
      "month_matched_photographic_null_per_taxon":NULL_REPS,
      "genus_bootstrap_repetitions":GENUS_REPS,
      "all_4475_source_taxa_present_in_each_of_four_panels":True,
      "source_biological_conclusions_deliberately_bounded":True,
      "panels":panels,
    }
    (a.out/"fcp4475_observer_disjoint_crossfit_white_thermal_results.json").write_text(
        json.dumps(summary,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("FCP_CROSSFIT_WHITE_THERMAL_COMPLETE",json.dumps({
      panel:{"n":q["n_heldout_thermal_estimable"],
             "three_class_groups":{k:{"n":v.get("n_taxa"),
                                     "mean":v.get("mean"),
                                     "CI":v.get("genus_bootstrap_marginal95")}
                                   for k,v in q["three_image_baseline_group"]["results"].items()},
             "adjusted":q["three_image_baseline_group"]["regression"].get("coef",{}).get(
                 "crossfit_mixed_vs_none"),
             "strict5":q["five_image_baseline_group"]["regression"].get("coef",{}).get(
                 "crossfit_mixed_vs_none")}
      for panel,q in panels.items()
    },sort_keys=True,allow_nan=False))


if __name__=="__main__":
    main()
