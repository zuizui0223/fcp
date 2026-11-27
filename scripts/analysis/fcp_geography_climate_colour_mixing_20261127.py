#!/usr/bin/env python3
"""White/nonwhite photographic neighbourhood mixing in geographic vs climatic space.

One fixed paired-space comparison on the original FCP census. These are
rank-neighbourhood mixing indices, not estimates of fundamental niche
hypervolumes, plant abundances, mating populations or selection.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from fcp_full4475_hypothesis_census_20261021 import (
    load_original, sample_climate, all_coords_ok, SOURCE_SHA, BIO_ARCHIVE_SHA,
)
from fcp_all4475_seven_abiotic_axes_20261021 import add_srad, AXES, SRAD_ZIP_SHA
from fcp_spatial_environment_architecture_20261127 import climate_domains, DOMAINS
from fcp_white_chromatic_spatial_partition_20261127 import xyz

SEED=20261127
B=5000
CLOSE_FRACTION=.25
MIN_NEIGHBOUR_PAIR_MASS=30
MIN_EXPECTED_CROSS_PAIRS=5
MIN_PER_COLOUR=5
MIN_RANK_RESIDUAL_FRACTION=.05
EXPECTED_THINNED=218709
EXPECTED_CLIMATE=218682

def check(ok,msg):
    if not bool(ok):
        raise RuntimeError("FCP_PAIRED_SPACE_OVERLAP: "+msg)

def close_weights(distance):
    """Select exactly the nearest quarter of pair opportunities, fractionally
    weighting any ties at the boundary, without arbitrary photo-ID tie breaks.
    """
    n=len(distance)
    target=CLOSE_FRACTION*n
    boundary=np.sort(distance)[int(np.ceil(target))-1]
    low=distance<boundary
    equal=distance==boundary
    w=low.astype(float)
    w[equal]=(target-low.sum())/equal.sum()
    check(np.all(w>=-1e-12) and np.all(w<=1+1e-12)
          and abs(w.sum()-target)<1e-8,"nearest-quarter mass not conserved")
    return w

def month_stats(g):
    n=len(g)
    if n<3 or g.white.nunique()<2:return None
    i,j=np.triu_indices(n,1)
    good=(g.observer_id.to_numpy()[i]!=g.observer_id.to_numpy()[j])&(
          g.observation_id.to_numpy()[i]!=g.observation_id.to_numpy()[j])
    i,j=i[good],j[good]
    m=len(i)
    if m<3:return None
    pos=xyz(g.latitude,g.longitude)
    geo=2*6371.0088*np.arcsin(np.clip(np.linalg.norm(pos[i]-pos[j],axis=1)/2,0,1))
    env=g[DOMAINS].to_numpy(float)
    ed=np.mean((env[i]-env[j])**2,axis=1)
    # An all-identical distance axis is not an identifiable neighbourhood.
    if np.ptp(geo)<=1e-8 or np.ptp(ed)<=1e-12:return None
    w=g.white.to_numpy(int)
    y=(w[i]!=w[j]).astype(float)
    p=int(w.sum())
    expect=2*p*(n-p)/(n*(n-1))
    mass=CLOSE_FRACTION*m
    wg=close_weights(geo);we=close_weights(ed)
    # Exact uniform vertex-label randomization expectation on the same graph.
    # Both neighbourhoods have identical opportunity and expected mismatch.
    expected=mass*expect
    obsgeo=float(wg@y);obsenv=float(we@y)
    # Separate geographic and climatic distance ranks on the same pairs.
    # Average ranks handle ties. Within-month centering removes differences
    # in seasonal colour composition and mean distance ranks.
    x=np.column_stack((rankdata(geo,method="average"),
                       rankdata(ed,method="average")))/m
    x=x-x.mean(0)
    yc=y-y.mean()
    # Retain the same pair weights in both rank dimensions.
    return dict(n=n,n_white=p,n_nonwhite=n-p,m=m,mass=mass,
                expected=expected,obsgeo=obsgeo,obsenv=obsenv,
                xx=x.T@x,xy=x.T@yc,yy=float(yc@yc))

def summarize_taxon(g):
    row={"n_climate_complete":len(g),"status":"no_comparable_mixed_month",
         "n_eligible_months":0,"n_used_photos":0}
    parts=[r for _,sub in g.groupby("month",sort=True)
           if (r:=month_stats(sub)) is not None]
    if not parts:return row
    def sumkey(k):return sum(x[k] for x in parts)
    mass=sumkey("mass");expect=sumkey("expected")
    ng=sumkey("obsgeo");ne=sumkey("obsenv")
    nw=sumkey("n_white");nc=sumkey("n_nonwhite")
    row.update(n_eligible_months=len(parts),n_used_photos=sumkey("n"),
               n_white_used=nw,n_nonwhite_used=nc,n_all_pairs=sumkey("m"),
               nearest_pair_mass=mass,expected_cross_nearest_pairs=expect,
               observed_cross_geo=ng,observed_cross_climate=ne)
    row["status"]="insufficient_two_colour_neighbour_opportunity"
    if (mass<MIN_NEIGHBOUR_PAIR_MASS or expect<MIN_EXPECTED_CROSS_PAIRS or
        min(nw,nc)<MIN_PER_COLOUR):return row
    # 1 = random colour mixing under the month-specific composition;
    # >1 = excess mixing. This is not a bounded geometric overlap metric.
    row.update(status="neighbour_overlap_estimable",
               geographic_mixing=ng/expect,climate_mixing=ne/expect,
               climate_minus_geographic_mixing=(ne-ng)/expect,
               geographic_cross_share=ng/mass,climate_cross_share=ne/mass,
               expected_cross_share=expect/mass,
               raw_climate_minus_geographic_cross_pp=100*(ne-ng)/mass)
    xx=sumkey("xx");xy=sumkey("xy")
    check(xx.shape==(2,2) and np.all(np.diag(xx)>0),"degenerate rank design")
    rho=xx[0,1]/np.sqrt(xx[0,0]*xx[1,1])
    residfrac=max(0.,1-float(rho)**2)
    row["rank_distance_residual_fraction"]=residfrac
    row["conditional_rank_status"]="rank_axes_too_collinear"
    if residfrac<MIN_RANK_RESIDUAL_FRACTION:return row
    beta=np.linalg.solve(xx,xy)
    check(np.isfinite(beta).all(),"nonfinite rank fit")
    # A 10-percentile increase in each pair-distance rank, holding other rank
    # fixed, expressed as probability percentage points.
    row.update(conditional_rank_status="estimable",
               geo_rank_beta_pp_per10pct=10*float(beta[0]),
               climate_rank_beta_pp_per10pct=10*float(beta[1]),
               climate_minus_geo_rank_beta_pp_per10pct=10*float(beta[1]-beta[0]))
    return row

def summarize(frame,cols):
    if not len(frame):return None
    x=frame[cols].to_numpy(float)
    check(np.isfinite(x).all(),"summary has missing estimates")
    groups=frame.groupby("genus",sort=True)
    sums=groups[cols].sum().to_numpy(float)
    ns=groups.size().to_numpy(int)
    rng=np.random.default_rng(SEED)
    draws=np.empty((B,len(cols)))
    for b in range(B):
        idx=rng.integers(len(ns),size=len(ns))
        draws[b]=sums[idx].sum(0)/ns[idx].sum()
    ci=np.quantile(draws,[.025,.975],axis=0)
    return {"n_species":len(frame),"n_genera":len(ns),
            "n_photos_in_eligible_months":int(frame.n_used_photos.sum()),
            "mean":{k:float(x[:,j].mean()) for j,k in enumerate(cols)},
            "conditional_genus_bootstrap95":{
                k:[float(ci[0,j]),float(ci[1,j])] for j,k in enumerate(cols)},
            "median":{k:float(np.median(x[:,j])) for j,k in enumerate(cols)},
            "fraction_positive":{k:float((x[:,j]>0).mean()) for j,k in enumerate(cols)}}

def self_test():
    check(np.allclose(close_weights(np.array([0.,0.,1.,2.])),[.5,.5,0,0]),
          "ties not split uniformly")
    check(np.allclose(close_weights(np.array([0.,1.,2.,3.,4.])),[1,.25,0,0,0]),
          "fractional nearest-quarter mass")
    # Climate perfectly separates colour while geography mixes it.
    g=pd.DataFrame({"latitude":[0.]*8,"longitude":np.arange(8)*.01,
        "white":[0,1,0,1,0,1,0,1],"observer_id":np.arange(8),
        "observation_id":np.arange(8),"month":[5]*8})
    for k in DOMAINS:g[k]=g.white
    a=month_stats(g)
    check(a is not None and a["obsenv"]==0 and a["obsgeo"]>0,
          "synthetic climatic sorting not recovered")
    check(abs(a["expected"]-a["mass"]*4/7)<1e-12,
          "exact vertex mismatch expectation")
    # Null expectation must not differ across the two neighbourhood graphs.
    check(abs(a["mass"]-7)<1e-12,"synthetic pair mass")

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--photos",type=Path,required=True)
    p.add_argument("--bio",type=Path,required=True)
    p.add_argument("--srad",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    self_test()
    t,d=load_original(a.photos)
    check(len(t)==447500 and t.inat_taxon_id.nunique()==4475 and len(d)==227107,
          "original source denominators")
    reg=t[["inat_taxon_id","species","fcp_source_cohort"]].drop_duplicates("inat_taxon_id").copy()
    reg["genus"]=reg.species.str.split().str[0]
    d["date"]=pd.to_datetime(d.observed_on,format="mixed",errors="coerce",utc=True)
    valid=all_coords_ok(d)&d.date.notna().to_numpy()&d.observer_id.notna().to_numpy()&d.observation_id.notna().to_numpy()
    q=d.loc[valid].sort_values("photo_id",kind="stable").drop_duplicates(
        ["inat_taxon_id","observer_id"]).copy()
    check(len(q)==EXPECTED_THINNED,"prior source geo-date observer filter drift")
    q["month"]=q.date.dt.month.astype(int)
    q["white"]=q.morph.eq("white").astype(int)
    q=sample_climate(q,a.bio)
    q=add_srad(q,a.srad)
    q=q.loc[np.isfinite(q[list(AXES)].to_numpy(float)).all(1)].copy()
    check(len(q)==EXPECTED_CLIMATE,"prior complete-climate frame drift")
    scaling=climate_domains(q)
    print("Exact source frame reproduced",len(t),len(d),len(q),flush=True)
    records=[{"inat_taxon_id":int(tid),**summarize_taxon(g)}
             for tid,g in q.groupby("inat_taxon_id",sort=True)]
    ledger=reg.merge(pd.DataFrame(records),on="inat_taxon_id",how="left",
                     validate="one_to_one")
    ledger["status"]=ledger.status.fillna("no_climate_complete_photo")
    check(len(ledger)==4475 and ledger.inat_taxon_id.is_unique,"registry coverage")
    good=ledger.loc[ledger.status.eq("neighbour_overlap_estimable")].copy()
    matched=good.loc[good.conditional_rank_status.eq("estimable")].copy()
    cols=["geographic_mixing","climate_mixing","climate_minus_geographic_mixing",
          "geographic_cross_share","climate_cross_share","expected_cross_share",
          "raw_climate_minus_geographic_cross_pp"]
    rankcols=["geo_rank_beta_pp_per10pct","climate_rank_beta_pp_per10pct",
              "climate_minus_geo_rank_beta_pp_per10pct"]
    check(len(good)>30,"inadequate macro-neighbourhood opportunities")
    check(np.allclose(good.climate_mixing-good.geographic_mixing,
                      good.climate_minus_geographic_mixing),"overlap contrast identity")
    result={"schema":"fcp_original_photo_geographic_climatic_colour_mixing_v1",
            "status":"ORIGINAL_PHOTO_CLIMATE_HASH_VERIFIED",
            "source_sha256":SOURCE_SHA,"bio_sha256":BIO_ARCHIVE_SHA,"srad_sha256":SRAD_ZIP_SHA,
            "source_photos":len(t),"source_species":len(reg),"classified_photos":len(d),
            "climate_complete_thinned_photos":len(q),
            "design":{"response":"white vs any nonwhite frozen categorical photo labels",
                      "same_calendar_month_across_years":True,
                      "one_photo_per_species_observer":True,
                      "geographic_distance":"great-circle distance",
                      "climate_distance":"mean squared difference in five globally standardized climate domains",
                      "common_nearest_pair_fraction":CLOSE_FRACTION,
                      "ties":"fractional weights at boundary preserve exact equal pair mass",
                      "minimum_neighbour_pair_mass":MIN_NEIGHBOUR_PAIR_MASS,
                      "minimum_expected_cross_pairs":MIN_EXPECTED_CROSS_PAIRS,
                      "minimum_photos_each_colour":MIN_PER_COLOUR,
                      "mixing_index":"observed nearest-quarter opposite-colour pair mass / exact same-month random-label expected mass",
                      "primary":"climate minus geographic mixing index on identical eligible species and month-specific pair opportunities",
                      "joint_model":"two within-month distance ranks simultaneously predict opposite-colour pair probability, slopes per 10-percentile rank increment",
                      "joint_rank_min_residual_fraction":MIN_RANK_RESIDUAL_FRACTION,
                      "uncertainty":"whole-genus 5000 bootstrap conditional on fixed species-level photo/climate graphs"},
            "global_climate_scaling":scaling,
            "status_counts":{str(k):int(v) for k,v in ledger.status.value_counts().items()},
            "primary_neighbourhood_comparison":summarize(good,cols),
            "joint_rank_support":{str(k):int(v) for k,v in good.conditional_rank_status.value_counts().items()},
            "same_pairs_joint_rank_decomposition":summarize(matched,cols+rankcols),
            "limitations":["Photographed white/nonwhite labels do not establish hereditary morphs or homologous petal identity.",
              "The mixing index is relative to sampled same-month composition and nearest-quarter pair opportunity; it is not Schoener D, bounded niche overlap or an estimated fundamental niche.",
              "Photo neighbourhoods are relative pair-distance ranks, not fixed pollinator neighbourhoods, breeding populations or exact habitat overlap.",
              "Only species-months containing both recorded labels and variable distances contribute; this is a conditional comparison among supported species, not all 4475.",
              "WorldClim normals at sampled locations do not capture temporal anomalies or all microhabitat/environmental differences.",
              "Joint rank slopes separate descriptive distance associations, not causal effects or migration, and can retain nonlinear geographic confounding.",
              "Pairs share photos and genus intervals condition on the fixed sampled graphs; neither raw pair count nor observer count is independent biological replication.",
              "No fitness, allele frequency, adaptive selection or frequency-dependent maintenance inferred."]}
    ledger.to_csv(a.out/"fcp_geography_climate_mixing_species_ledger.csv",index=False)
    (a.out/"fcp_geography_climate_mixing_result.json").write_text(
        json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("FCP_SPACE_MIXING_RESULT_BEGIN",flush=True)
    print(json.dumps(result,indent=2,allow_nan=False),flush=True)
    print("FCP_SPACE_MIXING_RESULT_END",flush=True)

if __name__=="__main__":main()
