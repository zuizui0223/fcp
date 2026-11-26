#!/usr/bin/env python3
"""Own FCP photo census: within-species white label vs other-taxon colour context.

One model. No external species-trait literature, no reconstructed phenotypes,
no alternative radius search. Context is equal-species photographed appearance,
not abundance, actual co-flowering events, interaction or reproductive fitness.
"""
from __future__ import annotations
import argparse, hashlib, json, math
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from fcp_full4475_hypothesis_census_20261021 import (
    load_original, sample_climate, all_coords_ok, SOURCE_SHA, TOTAL_TAXA,
    TOTAL_ROWS, CLASSIFIABLE_ROWS,
)
from fcp_all4475_seven_abiotic_axes_20261021 import add_srad, AXES, SRAD_ZIP_SHA
from fcp4475_shared_thermal_moisture_photo_hues_20261024 import compose_macros
from fcp4475_nine_colour_seven_climate_covariance_20261024 import design_geo_month

SEED=20261126
RADIUS=50.
MIN_REFS=5
MIN_OUTSIDE=2
MIN_FOCAL=25
MIN_CELLS=5
MIN_DF=8
MIN_RESID_SD=.03
MIN_RETAIN=.05
BOOT=5000

def require(ok,msg):
    if not bool(ok): raise RuntimeError(msg)

def xyz(lat,lon):
    la=np.deg2rad(np.asarray(lat,float));lo=np.deg2rad(np.asarray(lon,float))
    c=np.cos(la)
    return np.column_stack((c*np.cos(lo),c*np.sin(lo),np.sin(la)))

def fold(genus):
    return hashlib.sha256(("fcp-context-v1:"+genus).encode()).digest()[0]%2

def context_values(codes,white,observers,local_idx,focal_obs,tot,tw,obs_index):
    """Local and OUTSIDE-50km within-month means over identical reference taxa.

    Remove every record of the focal photographer, both local and outside.
    Each retained reference species needs >=1 local and >=2 outside photos.
    Reference observations are already thinned to one per taxon/observer.
    """
    li=np.asarray(local_idx,dtype=int)
    li=li[observers[li]!=focal_obs]
    if len(li)<MIN_REFS:return None
    spp,inv=np.unique(codes[li],return_inverse=True)
    n=np.bincount(inv); w=np.bincount(inv,weights=white[li])
    # Only reference records of this observer, at most one per taxon.
    oi=obs_index.get(focal_obs)
    if oi is not None:
        pos=np.searchsorted(spp,codes[oi])
        good=pos<len(spp)
        pp=pos[good];oo=oi[good]
        good2=spp[pp]==codes[oo]
        pp=pp[good2];oo=oo[good2]
        ob_n=np.bincount(pp,minlength=len(spp))
        ob_w=np.bincount(pp,weights=white[oo],minlength=len(spp))
    else:
        ob_n=np.zeros(len(spp));ob_w=np.zeros(len(spp))
    outside=tot[spp]-n-ob_n
    outsidew=tw[spp]-w-ob_w
    require(np.all(outside>=0) and np.all(outsidew>=-1e-8) and
            np.all(outsidew<=outside+1e-8),"outside-count partition error")
    good=outside>=MIN_OUTSIDE
    if good.sum()<MIN_REFS:return None
    keep=good[inv]
    if len(np.unique(observers[li[keep]]))<MIN_REFS:return None
    local=float(np.mean(w[good]/n[good]))
    baseline=float(np.mean(outsidew[good]/outside[good]))
    return local,baseline,local-baseline,int(good.sum()),int(n[good].sum())

def compute_context(d):
    coords=xyz(d.latitude,d.longitude)
    out=np.full((len(d),5),np.nan)
    month=d.month.to_numpy(int); folds=d.fold.to_numpy(int)
    tax=d.tax_code.to_numpy(int); obs=d.obs_code.to_numpy(int)
    w=d.white.to_numpy(float)
    chord=2*np.sin(RADIUS/(2*6371.0088))
    for mo in range(1,13):
        for fd in (0,1):
            q=np.flatnonzero((month==mo)&(folds==fd))
            ref=np.flatnonzero((month==mo)&(folds!=fd))
            if not len(q) or not len(ref):continue
            rt=tax[ref];rw=w[ref];ro=obs[ref]
            tot=np.bincount(rt,minlength=TOTAL_TAXA)
            tw=np.bincount(rt,weights=rw,minlength=TOTAL_TAXA)
            order=np.argsort(ro,kind="stable")
            cuts=np.flatnonzero(np.diff(ro[order]))+1
            obs_index={int(ro[a[0]]):a for a in np.split(order,cuts) if len(a)}
            tree=cKDTree(coords[ref])
            for off in range(0,len(q),256):
                ix=q[off:off+256]
                pairs=tree.query_ball_point(coords[ix],chord)
                for i,li in zip(ix,pairs):
                    v=context_values(rt,rw,ro,li,int(obs[i]),tot,tw,obs_index)
                    if v is not None:out[i]=v
            print("context",mo,fd,"focal",len(q),"supported",
                  int(np.isfinite(out[q,0]).sum()),flush=True)
    return out

def resid(x,C):
    return x-C@np.linalg.lstsq(C,x,rcond=1e-8)[0]

def fit_taxon(g):
    row={"inat_taxon_id":int(g.inat_taxon_id.iloc[0]),"species":str(g.species.iloc[0]),
         "genus":str(g.genus.iloc[0]),"fold":int(g.fold.iloc[0]),
         "n_observer_thinned_geo_dated":len(g),
         "n_context":int(g.local_white.notna().sum()),
         "n_complete":0,"status":"too_few_focal_photos"}
    ok=g.local_white.notna()&np.isfinite(g[list(AXES)].to_numpy(float)).all(axis=1)
    g=g.loc[ok].reset_index(drop=True); n=len(g)
    row["n_complete"]=n
    if n<MIN_FOCAL:return row
    cells=np.column_stack((np.floor((90-g.latitude.to_numpy())*6),
                           np.floor((g.longitude.to_numpy()+180)*6)))
    nc=len(np.unique(cells,axis=0));row["n_cells"]=nc
    if nc<MIN_CELLS:
        row["status"]="too_few_cells";return row
    env=g[list(AXES)].to_numpy(float)
    sd=env.std(axis=0);mu=env.mean(axis=0)
    Z=(env-mu)/np.where(sd>1e-8,sd,1)
    G,_=design_geo_month(g)
    # Include focal climate, geography, month, and context sampling opportunity.
    C=np.column_stack((G,compose_macros(Z),np.log1p(g.n_ref_species),
                      np.log1p(g.n_ref_photos)))
    X=g[["outside_white","local_departure"]].to_numpy(float)
    y=g.white.to_numpy(float)
    rank=np.linalg.matrix_rank(np.column_stack((C,X)),tol=1e-7)
    df=n-rank;row["df"]=int(df)
    if df<MIN_DF:
        row["status"]="insufficient_residual_df";return row
    # Require both distinct context axes to retain interpretable variation.
    rs=[]
    for j in (0,1):
        rr=resid(X[:,j],np.column_stack((C,X[:,1-j])))
        var=float(np.mean(rr*rr))
        raw=float(np.var(X[:,j]))
        rs.append((math.sqrt(var),var/max(raw,1e-15)))
    row.update(outside_resid_sd=rs[0][0],departure_resid_sd=rs[1][0],
               outside_retained=rs[0][1],departure_retained=rs[1][1])
    if any(a<MIN_RESID_SD or b<MIN_RETAIN for a,b in rs):
        row["status"]="context_axes_not_jointly_identifiable";return row
    R=resid(X,C);Y=resid(y,C)
    beta=np.linalg.lstsq(R,Y,rcond=1e-8)[0]
    require(np.isfinite(beta).all(),"nonfinite slope")
    row.update(status="estimable",n_white=int(y.sum()),
               n_nonwhite=int(n-y.sum()),outside_beta_pp_per10pp=float(10*beta[0]),
               departure_beta_pp_per10pp=float(10*beta[1]),
               mean_local_white=float(g.local_white.mean()),
               mean_outside_white=float(g.outside_white.mean()),
               mean_n_reference_species=float(g.n_ref_species.mean()))
    return row

def summary(rows):
    if not len(rows):return {"n_species":0}
    a=rows[["outside_beta_pp_per10pp","departure_beta_pp_per10pp"]].to_numpy(float)
    groups=rows.groupby("genus",sort=True)[["outside_beta_pp_per10pp","departure_beta_pp_per10pp"]]
    sums=groups.sum().to_numpy(float);counts=groups.size().to_numpy(int)
    rng=np.random.default_rng(SEED)
    sims=[]
    for _ in range(BOOT):
        ix=rng.integers(len(sums),size=len(sums))
        sims.append(sums[ix].sum(axis=0)/counts[ix].sum())
    ci=np.quantile(np.asarray(sims),[.025,.975],axis=0)
    z={"n_species":int(len(rows)),"n_genera":int(len(sums)),
       "n_focal_photos":int(rows.n_complete.sum()),
       "n_species_with_both_labels":int(((rows.n_white>0)&(rows.n_nonwhite>0)).sum())}
    for k,j in [("outside_species_baseline",0),("local_departure_primary",1)]:
        z[k]={"mean_pp_per10pp":float(a[:,j].mean()),
              "median_pp_per10pp":float(np.median(a[:,j])),
              "conditional_genus_bootstrap95":[float(ci[0,j]),float(ci[1,j])],
              "fraction_positive":float((a[:,j]>0).mean()),
              "fraction_negative":float((a[:,j]<0).mean())}
    return z

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--photos",type=Path,required=True)
    ap.add_argument("--bio",type=Path,required=True)
    ap.add_argument("--srad",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    t,d=load_original(a.photos)
    require(len(t)==447500 and t.inat_taxon_id.nunique()==4475 and len(d)==227107,"source counts")
    register=t[["inat_taxon_id","species"]].drop_duplicates("inat_taxon_id").copy()
    register["genus"]=register.species.str.split().str[0]
    register["fold"]=register.genus.map(fold)
    d["date"]=pd.to_datetime(d.observed_on,format="mixed",errors="coerce",utc=True)
    valid=all_coords_ok(d)&d.date.notna().to_numpy()&d.observer_id.notna().to_numpy()
    d=d.loc[valid].sort_values("photo_id",kind="stable").drop_duplicates(
        ["inat_taxon_id","observer_id"]).reset_index(drop=True)
    d["genus"]=d.species.str.split().str[0]
    d["fold"]=d.genus.map(fold)
    d["month"]=d.date.dt.month.astype(int)
    d["white"]=d.morph.eq("white").astype(float)
    ids={int(t):i for i,t in enumerate(sorted(register.inat_taxon_id))}
    d["tax_code"]=d.inat_taxon_id.map(ids).astype(int)
    d["obs_code"]=pd.factorize(d.observer_id,sort=True)[0]
    require(d.groupby("genus").fold.nunique().eq(1).all(),"genus split leakage")
    print("observer-thinned geo-dated photos",len(d),"taxa",d.inat_taxon_id.nunique(),flush=True)
    # Test exact counting on a constructed case, not scientific evidence.
    arr=np.arange(5).repeat(4);ww=np.tile([1.,0.,0.,1.],5);oo=np.arange(20)
    tt=np.bincount(arr);tw=np.bincount(arr,weights=ww)
    z=context_values(arr,ww,oo,np.arange(0,20,4),-1,tt,tw,{})
    require(z is not None and np.allclose(z[:3],[1,1/3,2/3]),"context counting test")
    vals=compute_context(d)
    cols=["local_white","outside_white","local_departure","n_ref_species","n_ref_photos"]
    d[cols]=vals
    valid_context=d.local_white.notna()
    require(np.allclose(d.loc[valid_context,"local_white"],
                        d.loc[valid_context,"outside_white"]+d.loc[valid_context,"local_departure"]),
            "local vs outside identity")
    print("context supported",int(valid_context.sum()),"taxa",
          d.loc[valid_context,"inat_taxon_id"].nunique(),flush=True)
    d=sample_climate(d,a.bio)
    d=add_srad(d,a.srad)
    rec=[fit_taxon(g) for _,g in d.groupby("inat_taxon_id",sort=True)]
    fits=pd.DataFrame(rec)
    ledger=register.merge(fits.drop(columns=["species","genus","fold"]),on="inat_taxon_id",how="left",validate="one_to_one")
    ledger["status"]=ledger.status.fillna("no_valid_observer_thinned_photo")
    require(len(ledger)==4475 and ledger.inat_taxon_id.is_unique,"ledger coverage")
    # Source-derived context records, without exposing photographer identity.
    d[["inat_taxon_id","photo_id","fold","month","white",*cols]].to_csv(
        a.out/"fcp_context_photo_support.csv.gz",index=False,compression="gzip")
    ledger.to_csv(a.out/"fcp_context_species_ledger.csv",index=False)
    good=ledger.loc[ledger.status.eq("estimable")].copy()
    result={"schema":"fcp_own_photo_other_species_colour_context_v1",
            "status":"ORIGINAL_SOURCE_HASH_VERIFIED",
            "source_sha256":SOURCE_SHA,"srad_sha256":SRAD_ZIP_SHA,
            "registry":{"source_photos":len(t),"species":len(register),"classified":227107,
                        "observer_thinned_geo_dated_photos":len(d),
                        "n_context_supported_photos":int(valid_context.sum()),
                        "n_context_supported_species":int(d.loc[valid_context,"inat_taxon_id"].nunique()),
                        "status_counts":{str(k):int(v) for k,v in ledger.status.value_counts().items()}},
            "design":{"radius_km":RADIUS,"calendar_month_same_not_same_year":True,
                      "genus_split_two_complementary_reference_folds":True,
                      "minimum_reference_species":MIN_REFS,"minimum_reference_outside_photos_same_month":MIN_OUTSIDE,
                      "minimum_focal_photos":MIN_FOCAL,"minimum_cells":MIN_CELLS,"minimum_df":MIN_DF,
                      "minimum_context_residual_sd":MIN_RESID_SD,"minimum_retained_variance":MIN_RETAIN,
                      "focal_observer_removed_from_context":True,
                      "reference_weighting":"equal other species",
                      "response":"binary original white photo category",
                      "predictors":"outside same-month reference species baseline and local-minus-outside white share; same reference species weights",
                      "controls":"five original climate domains, geography XYZ, calendar month, reference species and photo counts",
                      "units":"focal white probability percentage points per +10 percentage points context",
                      "bootstrap":"5000 whole-genus resamples conditional on fixed context; geographic and reciprocal-fold dependence not modeled"},
            "pooled":summary(good),"complementary_folds":{str(k):summary(good.loc[good.fold.eq(k)]) for k in (0,1)},
            "limitations":["Own archived iNaturalist original-photo measurement database, not independently counted plant or flower abundance.",
                           "Context is same calendar month across years and 50-km region, not actual synchronous co-flowering or shared pollinator neighborhood.",
                           "Disjoint reference genera and focal photographer removal prevent direct self inclusion, but complementary folds reuse the same archive in swapped roles.",
                           "Genus intervals condition on constructed context and do not propagate full reference-sample uncertainty or shared geographic dependence.",
                           "Associations may reflect shared environment, habitat, observer preferences and measurement artifacts; neither attraction nor rare-morph fitness is identified.",
                           "Original photo label not verified hereditary morph or petal spectrum."]}
    (a.out/"fcp_macro_colour_context_result.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("FCP_MACRO_CONTEXT_RESULT_BEGIN",flush=True)
    print(json.dumps(result,indent=2,allow_nan=False),flush=True)
    print("FCP_MACRO_CONTEXT_RESULT_END",flush=True)

if __name__=="__main__":main()
