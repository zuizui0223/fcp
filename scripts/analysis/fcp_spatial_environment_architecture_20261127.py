#!/usr/bin/env python3
"""Spatial environmental organization and within-species photographic colour replacement.

One comparative ecological test using the frozen 447,500-photo archive and its
existing WorldClim rasters. No literature trait enrichment, threshold search,
missing-colour imputation, or frequency-dependent fitness inference.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from fcp_full4475_hypothesis_census_20261021 import (
    load_original, sample_climate, all_coords_ok, SOURCE_SHA, BIO_ARCHIVE_SHA,
)
from fcp_all4475_seven_abiotic_axes_20261021 import add_srad, AXES, SRAD_ZIP_SHA
from fcp_white_chromatic_spatial_partition_20261127 import (
    xyz, one_taxon, LABELS, MIN_COMMON_PAIRS,
)

SEED=20261127
B=5000
DOMAINS=["thermal","precipitation","temperature_seasonality",
         "precipitation_seasonality","solar"]
EXPECTED_N=774
EXPECTED_DELTA=.01950300365085596

def check(test,msg):
    if not bool(test):
        raise RuntimeError("FCP_ENVIRONMENT_ARCHITECTURE: "+msg)

def climate_domains(q):
    # Common, colour-outcome-blind scaling across the thinned source photos.
    # First combine standardized BIO1/BIO5 and BIO12/BIO14, then scale each
    # of five domains to unit SD. Thus thermal and moisture each have one vote.
    z=q[list(AXES)].to_numpy(float)
    m=z.mean(0);s=z.std(0)
    check(np.isfinite(z).all() and np.all(s>0),"invalid source climate")
    z=(z-m)/s
    v=np.column_stack(((z[:,0]+z[:,2])/2,(z[:,3]+z[:,4])/2,
                        z[:,1],z[:,5],z[:,6]))
    vm=v.mean(0);vs=v.std(0)
    check(np.all(vs>0),"constant climate domain")
    v=(v-vm)/vs
    q[DOMAINS]=v
    return {"raw_axes":list(AXES),"raw_source_mean":m.tolist(),
            "raw_source_sd":s.tolist(),"domain_mean_before_final_scale":vm.tolist(),
            "domain_sd_before_final_scale":vs.tolist()}

def month_summary(g):
    n=len(g)
    if n<3:return None
    a,b=np.triu_indices(n,1)
    keep=(g.observer_id.to_numpy()[a]!=g.observer_id.to_numpy()[b])&(
          g.observation_id.to_numpy()[a]!=g.observation_id.to_numpy()[b])
    a,b=a[keep],b[keep]
    if len(a)==0:return None
    p=xyz(g.latitude.to_numpy(float),g.longitude.to_numpy(float))
    km=2*6371.0088*np.arcsin(np.clip(np.linalg.norm(p[a]-p[b],axis=1)/2,0,1))
    near=km<=50.;far=km>=250.
    nn=int(near.sum());nf=int(far.sum());K=min(nn,nf)
    if K==0:return None
    labels=g.label.to_numpy(int)
    mismatch=(labels[a]!=labels[b]).astype(float)
    c=np.bincount(labels,minlength=4)
    E=1-float(np.sum(c*(c-1)))/(n*(n-1))
    env=g[DOMAINS].to_numpy(float)
    envdist=np.mean((env[a]-env[b])**2,axis=1)
    # All original same-month observer-separated pairs define total breadth,
    # not merely the two extreme geographic bands.
    result=np.array([mismatch[near].mean(),mismatch[far].mean(),E,
                     envdist[near].mean(),envdist[far].mean(),envdist.mean(),
                     km[near].mean(),km[far].mean()],float)
    check(np.isfinite(result).all(),"nonfinite pair means")
    return K,nn,nf,result

def species(g):
    parts=[v for _,m in g.groupby("month",sort=True)
           if (v:=month_summary(m)) is not None]
    row={"n_climate_complete":len(g),"n_matched_months":len(parts),
         "status":"no_common_near_far_month"}
    if not parts:return row
    K=sum(x[0] for x in parts)
    L=sum(x[0]*x[3] for x in parts)/K
    names=["colour_near","colour_far","colour_expected","env_near",
           "env_far","env_breadth","near_mean_km","far_mean_km"]
    row.update(dict(zip(names,map(float,L))))
    row.update(common_pair_opportunity=int(K),
               near_pairs=sum(x[1] for x in parts),
               far_pairs=sum(x[2] for x in parts),
               colour_replacement_pp=100*float(L[1]-L[0]))
    row.update({f"mean_{k}":float(g[k].mean()) for k in DOMAINS})
    row["mean_abs_latitude"]=float(g.latitude.abs().mean())
    row["status"]="too_few_common_pair_opportunities"
    if K<MIN_COMMON_PAIRS:return row
    # Bounded environmental allocation contrast; undefined if both distance
    # bands have identical climatic values. It measures spatial partition,
    # while the separate env_breadth variable measures total magnitude.
    if L[3]+L[4]<=1e-12:
        row["status"]="no_sampled_climate_difference";return row
    row["env_spatial_partition"]=float((L[4]-L[3])/(L[4]+L[3]))
    check(-1.00000000001<=row["env_spatial_partition"]<=1.00000000001,
          "climate organization bound failed")
    row["status"]="estimable"
    return row

def regress(rows):
    d=rows.copy()
    y=d.colour_replacement_pp.to_numpy(float)
    d["log_env_breadth"]=np.log1p(d.env_breadth)
    d["log_near_distance"]=np.log1p(d.near_mean_km)
    d["log_far_distance"]=np.log1p(d.far_mean_km)
    d["log_common_pairs"]=np.log1p(d.common_pair_opportunity)
    d["log_photos"]=np.log1p(d.n_climate_complete)
    # Primary term is environmental spatial partition; breadth is separate.
    # Species-equal adjusted descriptive association, not a mechanism fit.
    features=["env_spatial_partition","log_env_breadth","colour_expected",
              "log_near_distance","log_far_distance","log_common_pairs",
              "log_photos","n_matched_months","mean_abs_latitude"]+[
              "mean_"+k for k in DOMAINS]
    v=d[features].to_numpy(float)
    mu=v.mean(0);sd=v.std(0)
    check(np.all(np.isfinite(v)) and np.all(sd>1e-10),"constant or invalid model predictors")
    z=(v-mu)/sd
    coh=pd.get_dummies(d.fcp_source_cohort.astype(str),drop_first=True,dtype=float)
    X=np.column_stack((np.ones(len(d)),z,coh.to_numpy(float)))
    names=["intercept"]+features+["cohort_"+str(k) for k in coh.columns]
    k=X.shape[1]
    check(len(d)>k+30 and np.linalg.matrix_rank(X)==k,"insufficient model rank")
    beta=np.linalg.lstsq(X,y,rcond=1e-10)[0]
    pred=X@beta
    # Diagnose residual support of the planned partition predictor,
    # rather than searching alternative covariates after seeing results.
    oth=np.delete(X,1,axis=1)
    pr=X[:,1]-oth@np.linalg.lstsq(oth,X[:,1],rcond=1e-10)[0]
    residual_fraction=float(np.mean(pr**2))
    groups=d.genus.to_numpy(str)
    gg=np.unique(groups)
    XX=[];XY=[]
    for g in gg:
        xx=X[groups==g];yy=y[groups==g]
        XX.append(xx.T@xx);XY.append(xx.T@yy)
    XX=np.asarray(XX);XY=np.asarray(XY)
    rng=np.random.default_rng(SEED)
    draws=[]
    rejected=0
    for _ in range(B):
        idx=rng.integers(len(gg),size=len(gg))
        a=XX[idx].sum(0);b=XY[idx].sum(0)
        if np.linalg.matrix_rank(a)<k:
            rejected+=1;continue
        draws.append(np.linalg.solve(a,b))
    check(len(draws)>=.99*B,"bootstrap model rank unstable")
    ci=np.quantile(np.asarray(draws),[.025,.975],axis=0)
    coef=[{"term":name,"estimate_pp":float(beta[j]),
           "conditional_genus_bootstrap95":[float(ci[0,j]),float(ci[1,j])]}
          for j,name in enumerate(names)]
    quart=np.quantile(v[:,0],[.25,.5,.75])
    # IQR contrast on original environmental partition scale, other variables fixed.
    iqr_sd=(quart[2]-quart[0])/sd[0]
    return {"n_species":len(d),"n_genera":len(gg),
            "n_focal_climate_complete_photos":int(d.n_climate_complete.sum()),
            "n_species_with_expected_colour_mismatch_zero":int(d.colour_expected.eq(0).sum()),
            "formula":"colour_replacement_pp ~ standardized(env_spatial_partition + log1p(env_breadth) + expected_colour_mismatch + log1p(near/far mean distance) + log1p(common pairs) + log1p(photos) + matched months + mean absolute latitude + five mean climate domains) + acquisition cohort",
            "coefficients":coef,"primary":coef[1],"breadth":coef[2],
            "primary_residual_variance_fraction":residual_fraction,
            "design_rank":int(k),"bootstrap_accepted":len(draws),
            "bootstrap_rejected":rejected,
            "primary_partition_q25_median_q75":quart.tolist(),
            "primary_partition_sd":float(sd[0]),
            "primary_iqr_contrast_pp":{"estimate":float(beta[1]*iqr_sd),
                "conditional_genus_bootstrap95":[float(ci[0,1]*iqr_sd),float(ci[1,1]*iqr_sd)]},
            "mean_colour_replacement_pp":float(y.mean()),
            "mean_colour_near":float(d.colour_near.mean()),
            "mean_colour_far":float(d.colour_far.mean()),
            "r2_in_sample":float(1-np.sum((y-pred)**2)/np.sum((y-y.mean())**2)),
            "raw_partition_breadth_correlation":float(np.corrcoef(v[:,0],v[:,1])[0,1])}

def self_test():
    g=pd.DataFrame({"latitude":[0.]*6,"longitude":[0,.01,.02,4,4.01,4.02],
        "label":[0,0,0,1,1,1],"observer_id":range(6),
        "observation_id":range(6),"month":[5]*6})
    for k in DOMAINS:g[k]=[0,0,0,1,1,1]
    K,nn,nf,x=month_summary(g)
    check((K,nn,nf)==(6,6,9),"synthetic matched pair geometry")
    check(np.allclose(x[:6],[0,1,.6,0,1,.6]),"synthetic environmental or colour split")
    row=species(g)
    check(row["colour_replacement_pp"]==100.,"synthetic scale contrast")

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--photos",type=Path,required=True)
    p.add_argument("--bio",type=Path,required=True)
    p.add_argument("--srad",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    self_test()
    t,d=load_original(args.photos)
    check(len(t)==447500 and t.inat_taxon_id.nunique()==4475 and len(d)==227107,"source denominators")
    reg=t[["inat_taxon_id","species","fcp_source_cohort"]].drop_duplicates("inat_taxon_id").copy()
    reg["genus"]=reg.species.str.split().str[0]
    d["date"]=pd.to_datetime(d.observed_on,format="mixed",errors="coerce",utc=True)
    valid=all_coords_ok(d)&d.date.notna().to_numpy()&d.observer_id.notna().to_numpy()&d.observation_id.notna().to_numpy()
    q=d.loc[valid].sort_values("photo_id",kind="stable").drop_duplicates(
        ["inat_taxon_id","observer_id"]).copy()
    q["month"]=q.date.dt.month.astype(int)
    q["label"]=pd.Categorical(q.morph,categories=LABELS).codes
    check(len(q)==218709 and (q.label>=0).all(),"prior thinned photographic frame drift")
    # Verify the previous 774-taxon comparison BEFORE adding climate eligibility.
    old=[]
    for tid,g in q.groupby("inat_taxon_id",sort=True):
        r=one_taxon(g)
        if r["status"] in ("both_axes_estimable","total_contrast_estimable_axes_insufficient"):
            old.append({"inat_taxon_id":int(tid),"original_replacement":r["raw_all_far_minus_near"]})
    old=pd.DataFrame(old)
    check(len(old)==EXPECTED_N and abs(old.original_replacement.mean()-EXPECTED_DELTA)<1e-10,
          "previous matched-band colour result did not reproduce")
    print("Previous exact matched geographic sample reproduced",len(old),float(old.original_replacement.mean()),flush=True)
    q=sample_climate(q,args.bio)
    q=add_srad(q,args.srad)
    validclim=np.isfinite(q[list(AXES)].to_numpy(float)).all(axis=1)
    q=q.loc[validclim].copy()
    scales=climate_domains(q)
    a=[]
    for tid,g in q.groupby("inat_taxon_id",sort=True):
        a.append({"inat_taxon_id":int(tid),**species(g)})
    ledger=reg.merge(pd.DataFrame(a),on="inat_taxon_id",how="left",validate="one_to_one")
    ledger=ledger.merge(old,on="inat_taxon_id",how="left",validate="one_to_one")
    ledger["status"]=ledger.status.fillna("no_climate_complete_photo")
    check(len(ledger)==4475 and ledger.inat_taxon_id.is_unique,"registry coverage")
    good=ledger.loc[ledger.status.eq("estimable")].copy()
    check(good.original_replacement.notna().all(),"climate filtering gained geographic opportunity")
    print("Eligible environmental-architecture taxa",len(good),flush=True)
    result={"schema":"fcp_original_photo_spatial_environment_organization_v1",
      "status":"HASH_VERIFIED_ORIGINAL_PHOTO_AND_WORLDCLIM",
      "source_sha256":SOURCE_SHA,"bio_sha256":BIO_ARCHIVE_SHA,"srad_sha256":SRAD_ZIP_SHA,
      "source_original_photos":447500,"source_species":4475,"source_classified":227107,
      "prior_matched_species_reproduced":len(old),
      "prior_colour_replacement_reproduced_pp":100*float(old.original_replacement.mean()),
      "climate_complete_thinned_photos":len(q),
      "status_counts":{str(k):int(v) for k,v in ledger.status.value_counts().items()},
      "design":{"near_km_le":50,"far_km_ge":250,"same_calendar_month_across_years":True,
         "pair_weights":"same month weight min(n_near,n_far) on both distance bands",
         "minimum_matched_pair_opportunity":MIN_COMMON_PAIRS,
         "environment_metric":"mean squared distance in five globally standardized fixed climate domains",
         "environment_partition":"(far mean env distance - near mean env distance)/(far + near); bounded [-1,1]",
         "environment_breadth":"same-month mean climate squared distance over all observer-separated pairs, matched-month weights",
         "response":"100*(far minus near probability of any original 4-state colour mismatch)",
         "primary":"adjusted species-equal slope per one SD environmental spatial partition, with breadth, climate means, observation opportunities and geography controlled",
         "uncertainty":"whole genus bootstrap, 5000 draws; fixed observed pair graphs, climate and colours"},
      "global_climate_scaling":scales,"model":regress(good),
      "limitations":["Original four-class photographic labels, not validated individual petal morphs.",
        "Same calendar month does not establish synchronous flowering or independent plant individuals.",
        "WorldClim normals describe spatial climate only, not interannual fluctuations during selection.",
        "Pair-distance bands and normalized environmental contrast are descriptive organization measures, not neutral genetic structure.",
        "Species-level controls do not remove all habitat, geography, shared species ranges, observer and measurement biases.",
        "Conditional genus bootstrap does not account for the full original photo sampling or cross-genus spatial dependence.",
        "This is an observational macroecological association, not a causal test of migration, frequency dependence, adaptation or fitness."]}
    ledger.to_csv(args.out/"fcp_environment_architecture_species_ledger.csv",index=False)
    (args.out/"fcp_environment_architecture_result.json").write_text(
        json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("FCP_ENV_ARCH_RESULT_BEGIN",flush=True)
    print(json.dumps(result,indent=2,allow_nan=False),flush=True)
    print("FCP_ENV_ARCH_RESULT_END",flush=True)

if __name__=="__main__":main()
