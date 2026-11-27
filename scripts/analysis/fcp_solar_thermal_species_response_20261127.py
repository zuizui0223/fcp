#!/usr/bin/env python3
"""Species-specific solar response: a direct between-observer prediction test.

Test whether signed white-fraction/solar slopes are reproducible across
disjoint photographers, rather than treating the magnitude of noisy single
fits as evidence for opposing biological responses. Thermal is a matched
benchmark. Source photo measurements, no assumed genotypes or fitness.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from fcp_full4475_hypothesis_census_20261021 import (
    load_original,sample_climate,all_coords_ok,SOURCE_SHA,BIO_ARCHIVE_SHA,
)
from fcp_all4475_seven_abiotic_axes_20261021 import add_srad,AXES,SRAD_ZIP_SHA
from fcp4475_shared_thermal_moisture_photo_hues_20261024 import compose_macros
from fcp4475_nine_colour_seven_climate_covariance_20261024 import design_geo_month
from fcp4475_crossfitted_observer_independent_thermal_white_20261024 import observer_folds

SEED=20261127
B=5000
TARGETS={"thermal":0,"solar":4}
MIN_PHOTOS=18
MIN_CELLS=3
MIN_RESID_DF=8
MIN_RETAIN=.05
# Common climate units across folds. Require residual support of at least
# 0.1 pooled within-species climate SD, not just a ratio of two tiny variances.
MIN_UNIQUE_SD=.10

def check(ok,msg):
    if not bool(ok):raise RuntimeError("FCP_SOLAR_SPECIES_RESPONSE: "+msg)

def common_domains(g):
    e=g[list(AXES)].to_numpy(float)
    mu=e.mean(0);sd=e.std(0)
    valid=sd>np.maximum(1e-7,2e-10*np.maximum(1.,np.abs(mu)))
    z=(e-mu)/np.where(valid,sd,1.)
    a=compose_macros(z)
    a=a-a.mean(0)
    s=a.std(0)
    valid=s>1e-8
    a=a/np.where(valid,s,1.)
    check(np.isfinite(a).all(),"invalid climate units")
    return a

def project(x,C):
    return x-C@np.linalg.lstsq(C,x,rcond=1e-8)[0]

def fold_fit(g,M):
    r={"n":len(g),"status_thermal":"too_few_photos",
       "status_solar":"too_few_photos"}
    if len(g)<MIN_PHOTOS:return r
    cells=np.column_stack((np.floor((90-g.latitude.to_numpy(float))*6),
                           np.floor((g.longitude.to_numpy(float)+180)*6)))
    nc=len(np.unique(cells,axis=0));r["cells"]=nc
    if nc<MIN_CELLS:
        r.update(status_thermal="too_few_cells",status_solar="too_few_cells")
        return r
    G,_=design_geo_month(g)
    R=project(M,G)
    X=np.column_stack((G,M))
    rank=np.linalg.matrix_rank(X,tol=1e-8)
    r["df"]=len(g)-int(rank)
    if r["df"]<MIN_RESID_DF:
        r.update(status_thermal="too_few_residual_df",status_solar="too_few_residual_df")
        return r
    y=np.stack(g.palette_norm9.to_numpy())[:,0]
    check(np.isfinite(y).all() and np.all((y>=0)&(y<=1)),
          "invalid measured white fraction")
    for name,j in TARGETS.items():
        u=project(R[:,j],np.delete(R,j,axis=1))
        uv=float(np.mean(u*u)); raw=float(np.var(M[:,j]))
        r[name+"_unique_sd"]=float(np.sqrt(max(uv,0)))
        r[name+"_retained"]=uv/max(raw,1e-15)
        if raw<=1e-10 or uv<MIN_UNIQUE_SD**2 or uv/raw<MIN_RETAIN:
            r["status_"+name]="not_identifiable_in_fold";continue
        # Since month fixed effects are in G, residual covariance already
        # subtracts the exact within-month shuffle expectation.
        beta=100*float(u@y/(u@u))
        check(np.isfinite(beta),"nonfinite heldout slope")
        r[name+"_slope_pp"]=beta
        r["status_"+name]="estimable"
    return r

def pearson(x,y):
    a=np.asarray(x,float)-np.mean(x);b=np.asarray(y,float)-np.mean(y)
    den=np.sqrt((a@a)*(b@b))
    return float(a@b/den) if den>1e-15 else np.nan

def metrics(a):
    """a columns: solar A,B and thermal A,B; same species, fixed climate units."""
    ans=[]
    for i in (0,2):
        x,y=a[:,i],a[:,i+1]
        # Rank correlation tests whether the species' *signed* solar response
        # predicts the other fold's signed response, robust to extreme slopes.
        r=pearson(rankdata(x),rankdata(y))
        raw=pearson(x,y)
        center=float(np.mean((x-x.mean())*(y-y.mean())))
        sx,sy=np.sign(x),np.sign(y)
        # Chance-corrected sign association; zeros have sign zero.
        sign_cov=float(np.mean(sx*sy)-sx.mean()*sy.mean())
        pred=float((np.mean(sx*y)+np.mean(sy*x))/2)
        # Remove constant directional expectation; any positive value is
        # out-of-fold species-specific sign information, not pooled mean.
        pred-=float((sx.mean()*y.mean()+sy.mean()*x.mean())/2)
        ans.extend([r,raw,center,sign_cov,pred,float(x.mean()),float(y.mean()),
                    float((x.mean()+y.mean())/2),
                    float(np.mean((sx==sy)&(sx!=0)))])
    ans.append(ans[0]-ans[9])
    return np.asarray(ans,float)

NAMES=[]
for name in ("solar","thermal"):
    NAMES += [name+"_"+k for k in (
        "signed_slope_spearman","signed_slope_pearson","centered_crossfold_covariance_pp2",
        "centered_sign_agreement","centered_sign_prediction_pp",
        "fold_A_mean_pp","fold_B_mean_pp","paired_mean_pp","same_nonzero_sign_fraction")]
NAMES.append("solar_minus_thermal_spearman")

def summary(d):
    cols=["solar_slope_pp_A","solar_slope_pp_B","thermal_slope_pp_A","thermal_slope_pp_B"]
    a=d[cols].to_numpy(float)
    check(np.isfinite(a).all(),"missing matched response")
    val=metrics(a)
    check(np.isfinite(val).all(),"degenerate observed response metric")
    ix=[np.asarray(v,int) for v in d.reset_index(drop=True).groupby("genus",sort=True).indices.values()]
    rng=np.random.default_rng(SEED)
    boot=[]
    for _ in range(B):
        draw=rng.integers(len(ix),size=len(ix))
        r=metrics(a[np.concatenate([ix[i] for i in draw])])
        if np.isfinite(r).all():boot.append(r)
    boot=np.asarray(boot)
    check(len(boot)>=.99*B,"unstable comparative bootstrap")
    ci=np.quantile(boot,[.025,.975],axis=0)
    primary_idx=[0,9,18] # two repeatabilities and their direct difference
    se=boot[:,primary_idx].std(0,ddof=1)
    check(np.all(se>0),"degenerate rank interval")
    crit=float(np.quantile(np.max(np.abs((boot[:,primary_idx]-val[primary_idx])/se),axis=1),.95))
    estimates={}
    for i,k in enumerate(NAMES):
        estimates[k]={"estimate":float(val[i]),"conditional_genus_bootstrap95":ci[:,i].tolist()}
        if i in primary_idx:
            j=primary_idx.index(i)
            estimates[k]["simultaneous95_three_rank_contrasts"]=[
                float(val[i]-crit*se[j]),float(val[i]+crit*se[j])]
    return {"n_species":len(d),"n_genera":len(ix),
            "n_photo_A":int(d.n_A.sum()),"n_photo_B":int(d.n_B.sum()),
            "n_solar_both_positive":int(((a[:,0]>0)&(a[:,1]>0)).sum()),
            "n_solar_both_negative":int(((a[:,0]<0)&(a[:,1]<0)).sum()),
            "n_solar_opposite_sign":int((a[:,0]*a[:,1]<0).sum()),
            "bootstrap_accepted":len(boot),"estimates":estimates}

def self_test():
    a=np.array([[-2,-3,-2,1],[-1,-1,-1,-1],[1,2,1,-2],[2,4,2,2]],float)
    z=metrics(a)
    check(abs(z[0]-1)<1e-12 and z[2]>0 and z[3]>0 and z[4]>0,
          "known bidirectional repeatability not recovered")
    # A shared positive mean alone must not be interpreted as heterogeneity.
    b=a.copy();b[:,0]+=10;b[:,1]+=10
    zz=metrics(b)
    check(abs(z[0]-zz[0])<1e-12 and abs(zz[3])<1e-12 and abs(zz[4])<1e-12,
          "constant sign is not species-specific information")

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--photos",type=Path,required=True)
    p.add_argument("--bio",type=Path,required=True)
    p.add_argument("--srad",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    self_test()
    t,d=load_original(a.photos)
    check(len(t)==447500 and len(d)==227107 and t.inat_taxon_id.nunique()==4475,
          "source denominators")
    reg=t[["inat_taxon_id","species","fcp_source_cohort"]].drop_duplicates("inat_taxon_id").copy()
    reg["genus"]=reg.species.str.split().str[0]
    d["date"]=pd.to_datetime(d.observed_on,format="mixed",errors="coerce",utc=True)
    valid=all_coords_ok(d)&d.date.notna().to_numpy()&d.observer_id.notna().to_numpy()&d.observation_id.notna().to_numpy()
    q=d.loc[valid].sort_values("photo_id",kind="stable").drop_duplicates(
        ["inat_taxon_id","observer_id"]).copy()
    check(len(q)==218709,"prior observer-thinned frame drift")
    q=sample_climate(q,a.bio);q=add_srad(q,a.srad)
    q=q.loc[np.isfinite(q[list(AXES)].to_numpy(float)).all(1)].copy()
    check(len(q)==218682,"prior climate frame drift")
    raw=t.groupby("inat_taxon_id",sort=True)
    rows=[]
    for tid,g0 in q.groupby("inat_taxon_id",sort=True):
        g=g0.reset_index(drop=True)
        f=g.observer_id.map(observer_folds(raw.get_group(tid))).to_numpy()
        check(np.isin(f,[0,1]).all(),"unassigned original observer")
        check(not(set(g.loc[f==0,"observer_id"])&set(g.loc[f==1,"observer_id"])),
              "observer leakage")
        M=common_domains(g)
        row={"inat_taxon_id":int(tid),"n_full_climate_photos":len(g)}
        for i,letter in enumerate(("A","B")):
            mask=f==i
            z=fold_fit(g.loc[mask].reset_index(drop=True),M[mask])
            row.update({k+"_"+letter:v for k,v in z.items()})
        rows.append(row)
    ledger=reg.merge(pd.DataFrame(rows),on="inat_taxon_id",how="left",validate="one_to_one")
    check(len(ledger)==4475 and ledger.inat_taxon_id.is_unique,"registry coverage")
    solar=ledger.status_solar_A.eq("estimable")&ledger.status_solar_B.eq("estimable")
    therm=ledger.status_thermal_A.eq("estimable")&ledger.status_thermal_B.eq("estimable")
    good=ledger.loc[solar&therm].copy()
    check(len(good)>30,"too few paired solar and thermal estimates")
    print("Source frame 447500/227107/218682; solar both, thermal both, common",
          int(solar.sum()),int(therm.sum()),len(good),flush=True)
    result={"schema":"fcp_solar_vs_thermal_signed_species_response_prediction_v1",
       "status":"ORIGINAL_PHOTO_CLIMATE_HASH_VERIFIED",
       "source":{"photos":len(t),"taxa":len(reg),"classified":len(d),
                 "climate_complete_thinned":len(q),"photo_sha256":SOURCE_SHA,
                 "bio_sha256":BIO_ARCHIVE_SHA,"srad_sha256":SRAD_ZIP_SHA},
       "design":{"question":"Are opposite species-specific solar/white responses reproducible on other photographers?",
           "observer_split":"existing deterministic species-specific observer_folds from all 100 original attempts, one valid classified photo per observer",
           "response":"original normalized nine-palette white fraction",
           "climate_units":"shared pooled species SD for each fixed five-domain climate score, computed using positions only, both folds",
           "model":"signed solar and thermal white-fraction slopes conditional on calendar month, XYZ, other four climate domains within each fold",
           "slope_units":"white photo-fraction percentage points per one pooled within-species climate-domain SD",
           "min_fold_photos":MIN_PHOTOS,"min_fold_cells":MIN_CELLS,"min_fold_residual_df":MIN_RESID_DF,
           "min_unique_fraction":MIN_RETAIN,"min_unique_sd_pooled_units":MIN_UNIQUE_SD,
           "primary":"signed solar slope Spearman crossfold correlation, on common solar/thermal support; matched thermal and direct difference form a simultaneous three-contrast family",
           "sign_prediction":"centered out-of-fold sign times opposite-fold signed slope, symmetrized; zero under independent directions after shared mean removal",
           "inference":"5000 whole-genus bootstrap with reranking within each draw; fixed observed climate/photo frames"},
       "support":{"solar_both_folds":int(solar.sum()),"thermal_both_folds":int(therm.sum()),
                  "common_both_domains_both_folds":len(good)},
       "common_species":summary(good),
       "limitations":["Repeatability concerns observed photo white fractions, not genetically verified petal morphs.",
         "Observer folds do not separate geographic populations, individual plants, platforms, or all measurement errors.",
         "Half-sample slopes may be noisy. A null repeatability estimate does not prove absence of species-specific biological responses.",
         "A positive repeatability could reflect repeated environmental or spatial confounding; not direct solar causation, UV protection or selection.",
         "Common-support species have enough residual climate variation in both folds and are not all 4475 taxa.",
         "The signed slope estimand differs from previous unsigned pair-distance divergence; neither algebraically implies the other."]}
    ledger.to_csv(a.out/"fcp_solar_thermal_observer_fold_species_ledger.csv",index=False)
    (a.out/"fcp_solar_thermal_species_response_result.json").write_text(
        json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("FCP_SOLAR_RESPONSE_RESULT_BEGIN",flush=True)
    print(json.dumps(result,indent=2,allow_nan=False),flush=True)
    print("FCP_SOLAR_RESPONSE_RESULT_END",flush=True)

if __name__=="__main__":main()
