#!/usr/bin/env python3
"""Same-species, calendar-month-matched regional vs nearby photographic colour.

One primary contrast: normalized regional replacement of white-vs-chromatic
labels versus between distinct chromatic labels. Original photo counts/positions,
no invented morphs or fitness. Distances are geographic, not populations.
"""
from __future__ import annotations
import argparse, hashlib, json, math
from pathlib import Path
import numpy as np
import pandas as pd

SHA="8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6"
SEED=20261127
BOOT=5000
MIN_COMMON_PAIRS=30
MIN_EXPECTED_AXIS_PAIRS=5
LABELS=["white","yellow_orange","red_pink","blue_purple"]

def digest(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(2**20),b""):h.update(b)
    return h.hexdigest()

def xyz(la,lo):
    a=np.deg2rad(la);b=np.deg2rad(lo);c=np.cos(a)
    return np.column_stack([c*np.cos(b),c*np.sin(b),np.sin(a)])

def month_stats(g):
    # All pair opportunities derive from original, independent observer/observation IDs.
    n=len(g)
    if n<3:return None
    i,j=np.triu_indices(n,1)
    good=(g.observer_id.to_numpy()[i]!=g.observer_id.to_numpy()[j])&(
        g.observation_id.to_numpy()[i]!=g.observation_id.to_numpy()[j])
    i,j=i[good],j[good]
    if not len(i):return None
    co=xyz(g.latitude.to_numpy(float),g.longitude.to_numpy(float))
    dist=2*6371.0088*np.arcsin(np.clip(np.linalg.norm(co[i]-co[j],axis=1)/2,0,1))
    a,b=g.label.to_numpy(int)[i],g.label.to_numpy(int)[j]
    near=dist<=50.;far=dist>=250.
    nn=int(near.sum());nf=int(far.sum())
    k=min(nn,nf)
    if not k:return None
    wc=(a==0)!=(b==0)
    cc=(a!=0)&(b!=0)&(a!=b)
    c=np.bincount(g.label.to_numpy(int),minlength=4)
    # Expected contrast probabilities from exact uniform vertex relabeling,
    # conditional on the original month-specific category counts.
    ew=2*c[0]*(n-c[0])/(n*(n-1))
    ec=(int(c[1:].sum())**2-int(np.sum(c[1:]**2)))/(n*(n-1))
    L=np.array([wc[near].mean(),cc[near].mean()],float)
    F=np.array([wc[far].mean(),cc[far].mean()],float)
    E=np.array([ew,ec],float)
    assert np.isclose(E.sum(),1-np.sum(c*(c-1))/(n*(n-1)))
    assert np.isclose(L.sum(),(a[near]!=b[near]).mean())
    assert np.isclose(F.sum(),(a[far]!=b[far]).mean())
    return k,nn,nf,L,F,E

def one_taxon(g):
    a=[r for _,m in g.groupby("month",sort=True) if (r:=month_stats(m)) is not None]
    row={"n_valid_thinned_photos":len(g),"n_matched_months":len(a),
         "status":"no_near_and_far_in_same_month"}
    if not a:return row
    K=sum(x[0] for x in a);N=sum(x[1] for x in a);F=sum(x[2] for x in a)
    L=sum(x[0]*x[3] for x in a)/K
    R=sum(x[0]*x[4] for x in a)/K
    E=sum(x[0]*x[5] for x in a)/K
    row.update(common_pair_opportunity=K,near_pairs=N,far_pairs=F,
        near_wc=float(L[0]),far_wc=float(R[0]),expected_wc=float(E[0]),
        near_cc=float(L[1]),far_cc=float(R[1]),expected_cc=float(E[1]),
        near_all=float(L.sum()),far_all=float(R.sum()),expected_all=float(E.sum()),
        raw_wc_far_minus_near=float(R[0]-L[0]),
        raw_cc_far_minus_near=float(R[1]-L[1]),
        raw_all_far_minus_near=float((R-L).sum()))
    row["status"]="too_few_common_pair_opportunities"
    if K<MIN_COMMON_PAIRS:return row
    row["status"]="total_contrast_estimable_axes_insufficient"
    # Use the label-conditioned expected opportunity, not observed regional
    # depletion or observed contrast sign, for the two-axis support gate.
    if np.any(K*E<MIN_EXPECTED_AXIS_PAIRS):return row
    D=(R-L)/E
    row.update(status="both_axes_estimable",wc_normalized_replacement=float(D[0]),
               cc_normalized_replacement=float(D[1]),
               wc_minus_cc_normalized_replacement=float(D[0]-D[1]))
    return row

def group_boot(frame,cols):
    if len(frame)==0:return None
    x=frame[cols].to_numpy(float)
    assert np.isfinite(x).all()
    gr=frame.groupby("genus",sort=True)
    sums=gr[cols].sum().to_numpy(float);ns=gr.size().to_numpy(int)
    rng=np.random.default_rng(SEED)
    sims=np.empty((BOOT,len(cols)))
    for b in range(BOOT):
        i=rng.integers(len(sums),size=len(sums))
        sims[b]=sums[i].sum(axis=0)/ns[i].sum()
    ci=np.quantile(sims,[.025,.975],axis=0)
    return {"n_species":int(len(frame)),"n_genera":len(sums),
            "means":{c:float(x[:,i].mean()) for i,c in enumerate(cols)},
            "genus_bootstrap95":{c:[float(ci[0,i]),float(ci[1,i])] for i,c in enumerate(cols)},
            "median":{c:float(np.median(x[:,i])) for i,c in enumerate(cols)},
            "fraction_positive":{c:float((x[:,i]>0).mean()) for i,c in enumerate(cols)}}

def self_test():
    # Two distance clusters, white in one and chromatic in the other;
    # expected white/chromatic difference is stronger across regions.
    g=pd.DataFrame({"latitude":[0.]*6,"longitude":[0,.01,.02,4,4.01,4.02],
        "label":[0,0,0,1,1,2],"observer_id":list(range(6)),
        "observation_id":list(range(6)),"month":[5]*6})
    a=month_stats(g)
    assert a is not None and a[0]==6 and a[1]==6 and a[2]==9
    assert np.allclose(a[3],[0,1/3]) and np.allclose(a[4],[1,0])
    assert np.allclose(a[5],[.6,2/15])

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    self_test()
    assert digest(a.source)==SHA
    fields=["inat_taxon_id","species","fcp_source_cohort","photo_id","observation_id",
       "observer_id","latitude","longitude","positional_accuracy_m","observed_on",
       "global_classifiable","morph"]
    d=pd.read_csv(a.source,usecols=fields,low_memory=False)
    assert len(d)==447500 and d.photo_id.nunique()==447500
    assert d.inat_taxon_id.nunique()==4475 and d.groupby("inat_taxon_id").size().eq(100).all()
    flag=d.global_classifiable.astype(str).str.strip().str.lower().isin(["true","yes","1"])
    assert int(flag.sum())==227107 and d.loc[flag,"morph"].isin(LABELS).all()
    reg=d[["inat_taxon_id","species","fcp_source_cohort"]].drop_duplicates("inat_taxon_id").copy()
    reg["genus"]=reg.species.str.split().str[0]
    d["date"]=pd.to_datetime(d.observed_on,format="mixed",errors="coerce",utc=True)
    valid=(flag&d.latitude.between(-90,90)&d.longitude.between(-180,180)&
           d.positional_accuracy_m.between(0,5000)&d.date.notna()&
           d.observer_id.notna()&d.observation_id.notna())
    q=d.loc[valid].sort_values("photo_id",kind="stable").drop_duplicates(
        ["inat_taxon_id","observer_id"]).copy()
    q["month"]=q.date.dt.month
    q["label"]=pd.Categorical(q.morph,categories=LABELS).codes
    assert (q.label>=0).all()
    rows=[]
    for tid,g in q.groupby("inat_taxon_id",sort=True):
        rows.append({"inat_taxon_id":int(tid),**one_taxon(g)})
    ledger=reg.merge(pd.DataFrame(rows),on="inat_taxon_id",how="left",validate="one_to_one")
    ledger["status"]=ledger.status.fillna("no_valid_thinned_photo")
    assert len(ledger)==4475 and ledger.inat_taxon_id.is_unique
    total=ledger.loc[ledger.status.isin(["both_axes_estimable","total_contrast_estimable_axes_insufficient"])].copy()
    both=ledger.loc[ledger.status.eq("both_axes_estimable")].copy()
    rawcols=["near_all","far_all","expected_all","raw_all_far_minus_near"]
    axescols=["near_wc","far_wc","expected_wc","near_cc","far_cc","expected_cc",
              "raw_wc_far_minus_near","raw_cc_far_minus_near",
              "wc_normalized_replacement","cc_normalized_replacement",
              "wc_minus_cc_normalized_replacement"]
    cols=rawcols+axescols
    result={"schema":"fcp_own_photo_same_species_axis_replacement_v1",
        "status":"ORIGINAL_447500_PHOTO_HASH_VERIFIED",
        "source_sha256":SHA,"source_n_photos":len(d),"source_species":len(reg),
        "source_classified":int(flag.sum()),"dated_geo_observer_thinned":len(q),
        "design":{"near_km_le":50,"regional_km_ge":250,"same_calendar_month_not_same_year":True,
          "month_weight":"min(near pair count, far pair count) for identical month weighting of both bands",
          "minimum_common_pair_opportunity":MIN_COMMON_PAIRS,
          "minimum_expected_pairs_each_axis":MIN_EXPECTED_AXIS_PAIRS,
          "white_axis":"white versus any chromatic category",
          "chromatic_axis":"two different nonwhite coarse categories",
          "expected_axis_frequency":"exact vertex-label shuffle expectation within original species-month",
          "normalized_replacement":"(regional minus nearby mismatch probability for axis) / expected month-matched mismatch probability",
          "primary":"white normalized replacement minus chromatic normalized replacement on identical eligible taxa",
          "uncertainty":"5000 whole-genus bootstrap intervals conditional on original photo pairs; no independent plant/population assumption"},
        "status_counts":{str(k):int(v) for k,v in ledger.status.value_counts().items()},
        "all_total_supported":group_boot(total,rawcols),
        "matched_axis_supported":group_boot(both,cols),
        "cohorts":{str(k):group_boot(g,cols) for k,g in both.groupby("fcp_source_cohort",sort=True)},
        "limitations":["Existing FCP photograph sample, not a census of individuals or known mating populations.",
          "50-km and 250-km bands are chosen once; same calendar month is not synchronous flowering.",
          "Different observer IDs and observation IDs do not guarantee different individual plants.",
          "The two-axis comparison requires photographed white and multiple chromatic categories; not all taxa can support it.",
          "The normalized differences are scale contrasts, not Fst or fitness/selection coefficients; sparse opportunities and photographic organ/colour mismeasurement remain.",
          "Genus bootstrap is conditional on original sampled geography and colours; not a neutral population-genetic or frequency-dependent selection null."]}
    ledger.to_csv(a.out/"fcp_spatial_axis_species_ledger.csv",index=False)
    (a.out/"fcp_spatial_axis_result.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("FCP_SPATIAL_AXIS_RESULT_BEGIN",flush=True)
    print(json.dumps(result,indent=2,allow_nan=False),flush=True)
    print("FCP_SPATIAL_AXIS_RESULT_END",flush=True)

if __name__=="__main__":main()
