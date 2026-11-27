#!/usr/bin/env python3
"""Resolve the five environmental distance axes on the SAME FCP photo pairs.

Frozen 447,500-photo source, previously identified 1,626-species mixing frame.
One six-predictor rank model within species, calendar-month intercepts removed.
A climate domain is not called dominant merely because its marginal interval
excludes zero: compare its coefficient with the other four directly.
No extra thresholds, external trait curation, or fitness imputation.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata
import fcp_geography_climate_colour_mixing_20261127 as previous
from fcp_full4475_hypothesis_census_20261021 import (
    load_original, sample_climate, all_coords_ok, SOURCE_SHA, BIO_ARCHIVE_SHA,
)
from fcp_all4475_seven_abiotic_axes_20261021 import add_srad, AXES, SRAD_ZIP_SHA
from fcp_spatial_environment_architecture_20261127 import climate_domains, DOMAINS
from fcp_white_chromatic_spatial_partition_20261127 import xyz

SEED=20261127
B=5000
PRED=["geography",*DOMAINS]
TERMS=[k+"_beta_pp_per10pct" for k in PRED]
CONTRAST="thermal_minus_other_four_mean_pp_per10pct"
MIN_RESIDUAL_FRACTION=.05

def check(ok,msg):
    if not bool(ok):
        raise RuntimeError("FCP_FIVE_ENVIRONMENTAL_DISTANCE_DOMAINS: "+msg)

def pair_stats(g):
    """Use EXACTLY the previous eligible month rule and pair graph."""
    old=previous.month_stats(g)
    if old is None:return None
    n=len(g)
    i,j=np.triu_indices(n,1)
    ids=g.observer_id.to_numpy()
    obs=g.observation_id.to_numpy()
    keep=(ids[i]!=ids[j])&(obs[i]!=obs[j])
    i,j=i[keep],j[keep]
    m=len(i)
    check(m==old["m"],"prior pair count drift")
    p=xyz(g.latitude,g.longitude)
    geo=2*6371.0088*np.arcsin(np.clip(np.linalg.norm(p[i]-p[j],axis=1)/2,0,1))
    clim=np.abs(g[DOMAINS].to_numpy(float)[i]-g[DOMAINS].to_numpy(float)[j])
    raw=np.column_stack((geo,clim))
    X=np.column_stack([rankdata(raw[:,k],method="average")/m
                       for k in range(raw.shape[1])])
    X-=X.mean(0)
    w=g.white.to_numpy(int)
    y=(w[i]!=w[j]).astype(float)
    y-=y.mean()
    return {"xx":X.T@X,"xy":X.T@y,"m":m,"n":n}

def fit(g,old):
    r={"n_complete_photos":len(g),"status":"outside_prior_1626_frame"}
    if old["status"]!="neighbour_overlap_estimable":return r
    r["prior_geo_rank"]=old.get("geo_rank_beta_pp_per10pct")
    r["prior_combined_climate_rank"]=old.get("climate_rank_beta_pp_per10pct")
    r["n_used_photos"]=old["n_used_photos"]
    a=[x for _,m in g.groupby("month",sort=True) if (x:=pair_stats(m)) is not None]
    check(sum(x["m"] for x in a)==old["n_all_pairs"],"original graph not reproduced")
    check(sum(x["n"] for x in a)==old["n_used_photos"],"original photo support changed")
    xx=sum(x["xx"] for x in a);xy=sum(x["xy"] for x in a)
    diag=np.diag(xx)
    r["n_all_pairs"]=int(old["n_all_pairs"])
    if np.any(diag<=1e-12):
        r["status"]="constant_distance_domain";return r
    scale=np.sqrt(diag)
    corr=xx/np.outer(scale,scale)
    corr=(corr+corr.T)/2
    if np.linalg.matrix_rank(corr,tol=1e-9)<len(PRED):
        r["status"]="rank_deficient_distance_domains";return r
    inv=np.linalg.inv(corr)
    residual=1/np.diag(inv)
    check(np.all(np.isfinite(residual)) and np.all(residual>=0),
          "nonfinite residual distance opportunity")
    for k,v in zip(PRED,residual):r[k+"_unique_rank_variance_fraction"]=float(v)
    if np.any(residual<MIN_RESIDUAL_FRACTION):
        r["status"]="one_or_more_domains_not_separately_identifiable";return r
    beta=np.linalg.solve(corr,xy/scale)/scale
    check(np.all(np.isfinite(beta)),"nonfinite distance rank model")
    r.update(status="all_six_axes_estimable",**dict(zip(TERMS,map(float,10*beta))))
    r[CONTRAST]=r[TERMS[1]]-float(np.mean([r[k] for k in TERMS[2:]]))
    return r

def summarize(frame):
    cols=TERMS+[CONTRAST]
    a=frame[cols].to_numpy(float)
    check(np.all(np.isfinite(a)),"missing fitted slope")
    gg=frame.groupby("genus",sort=True)
    sums=gg[cols].sum().to_numpy(float); ns=gg.size().to_numpy(int)
    rng=np.random.default_rng(SEED)
    draws=np.empty((B,len(cols)))
    for k in range(B):
        idx=rng.integers(len(ns),size=len(ns))
        draws[k]=sums[idx].sum(0)/ns[idx].sum()
    mu=a.mean(0)
    interval=np.quantile(draws,[.025,.975],axis=0)
    # Simultaneous family: geography, five unique climate axes and one
    # prespecified thermal-minus-other-climate contrast; same bootstrap draws.
    se=draws.std(0,ddof=1)
    check(np.all(se>0),"degenerate interval")
    critical=float(np.quantile(np.max(np.abs((draws-mu)/se),axis=1),.95))
    res=[]
    for i,c in enumerate(cols):
        res.append({"term":c,"mean_pp_per10pct":float(mu[i]),
                    "median":float(np.median(a[:,i])),
                    "fraction_positive":float((a[:,i]>0).mean()),
                    "conditional_genus_bootstrap95":interval[:,i].tolist(),
                    "simultaneous95":[float(mu[i]-critical*se[i]),float(mu[i]+critical*se[i])]})
    return {"n_species":len(frame),"n_genera":len(ns),
            "n_photos":int(frame.n_used_photos.sum()),
            "n_pair_opportunities":int(frame.n_all_pairs.sum()),
            "coefficients":res,"simultaneous_critical":critical}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--photos",type=Path,required=True)
    p.add_argument("--bio",type=Path,required=True)
    p.add_argument("--srad",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    previous.self_test()
    t,d=load_original(a.photos)
    check(len(t)==447500 and t.inat_taxon_id.nunique()==4475 and len(d)==227107,
          "original source denominators")
    reg=t[["inat_taxon_id","species","fcp_source_cohort"]].drop_duplicates("inat_taxon_id").copy()
    reg["genus"]=reg.species.str.split().str[0]
    d["date"]=pd.to_datetime(d.observed_on,format="mixed",errors="coerce",utc=True)
    valid=all_coords_ok(d)&d.date.notna().to_numpy()&d.observer_id.notna().to_numpy()&d.observation_id.notna().to_numpy()
    q=d.loc[valid].sort_values("photo_id",kind="stable").drop_duplicates(
        ["inat_taxon_id","observer_id"]).copy()
    check(len(q)==218709,"prior photo filter drift")
    q["month"]=q.date.dt.month.astype(int)
    q["white"]=q.morph.eq("white").astype(int)
    q=sample_climate(q,a.bio);q=add_srad(q,a.srad)
    q=q.loc[np.isfinite(q[list(AXES)].to_numpy(float)).all(1)].copy()
    check(len(q)==218682,"prior climate frame drift")
    scaling=climate_domains(q)
    rows=[]; oldrows=[]
    for tid,g in q.groupby("inat_taxon_id",sort=True):
        old=previous.summarize_taxon(g)
        if old["status"]=="neighbour_overlap_estimable":
            oldrows.append(old)
        rows.append({"inat_taxon_id":int(tid),**fit(g,old)})
    prev=pd.DataFrame(oldrows)
    check(len(prev)==1626 and prev.n_used_photos.sum()==69516,
          "prior mixing sample not reproduced")
    check(abs(prev.geo_rank_beta_pp_per10pct.mean()-.28464234138831496)<1e-10
          and abs(prev.climate_rank_beta_pp_per10pct.mean()-.15304256775243996)<1e-10,
          "prior joint rank coefficients not reproduced")
    ledger=reg.merge(pd.DataFrame(rows),on="inat_taxon_id",how="left",validate="one_to_one")
    ledger["status"]=ledger.status.fillna("no_climate_complete_photo")
    check(len(ledger)==4475 and ledger.inat_taxon_id.is_unique,"registry support")
    good=ledger.loc[ledger.status.eq("all_six_axes_estimable")].copy()
    check(len(good)>30,"insufficient comparative six-axis support")
    result={"schema":"fcp_own_photo_five_climate_distance_axes_v1",
            "status":"ORIGINAL_PHOTO_CLIMATE_HASH_VERIFIED",
            "photo_sha256":SOURCE_SHA,"bio_sha256":BIO_ARCHIVE_SHA,
            "srad_sha256":SRAD_ZIP_SHA,
            "source":{"n_photos":len(t),"n_species":len(reg),"classified":len(d),
                      "climate_complete_thinned":len(q)},
            "previous_joint_rank_reproduced":{"species":len(prev),
              "photos":int(prev.n_used_photos.sum()),
              "geographic_mean":float(prev.geo_rank_beta_pp_per10pct.mean()),
              "combined_climate_mean":float(prev.climate_rank_beta_pp_per10pct.mean())},
            "design":{"response":"opposite white/nonwhite photograph labels",
                       "predictors":PRED,
                       "distance":"geographic great-circle and each absolute climate-domain difference, ranked within species-calendar-month on identical pair opportunities",
                       "slope_units":"probability percentage points per +10 percentile rank, conditional on all other distance ranks",
                       "minimum_unique_rank_variance_fraction":MIN_RESIDUAL_FRACTION,
                       "support":"previous 1626 species restricted to all six axes jointly identifiable",
                       "primary":"joint five climate domain coefficients; direct thermal minus mean of other four slopes",
                       "inference":"species-equal mean, 5000 whole-genus bootstrap, marginal and simultaneous intervals across seven reported estimates"},
            "status_counts":{str(k):int(v) for k,v in ledger.status.value_counts().items()},
            "matched_domain_model":summarize(good),
            "limitations":["Conditional photographic associations, not measured inherited petal morphs or fitness.",
              "Distance ranks represent relative contrasts among sampled observations; they do not compare effects per physical unit.",
              "Separate domain slopes and prior combined climate slope are different estimands.",
              "A positive thermal slope denotes more opposite-colour pairs at larger temperature differences; it does not indicate which colour is favoured.",
              "Fixed source pair graphs share photo vertices and species share spatial settings; whole-genus intervals do not propagate all photographic or spatial uncertainty.",
              "A joint domain coefficient is not a causal climate effect, and the thermal contrast tests mean slopes rather than every possible pairwise difference."]}
    ledger.to_csv(a.out/"fcp_five_domain_distance_species_ledger.csv",index=False)
    (a.out/"fcp_five_domain_distance_result.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("FCP_FIVE_DOMAIN_DISTANCE_RESULT_BEGIN",flush=True)
    print(json.dumps(result,indent=2,allow_nan=False),flush=True)
    print("FCP_FIVE_DOMAIN_DISTANCE_RESULT_END",flush=True)

if __name__=="__main__":main()
