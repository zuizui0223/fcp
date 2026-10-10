#!/usr/bin/env python3
"""Post hoc 50-km local-pair support ladder for photographed FCP geographic allocation.

Replays the *same* source-locked 4-state statistic and per-species vertex
label permutations as the current paper. No new photos, colour classes,
geographic radii or biological hypotheses are introduced.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PREVIOUS = ROOT / "scripts/analysis/run_polymorphism_distributed_polymorphism_posthoc_20261007.py"
spec=importlib.util.spec_from_file_location("fcp_frozen_dist",PREVIOUS)
base=importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)

THRESHOLDS=(30,50,100,200,500,1000)
PRIMARY_STRESS_THRESHOLD=100
N_BOOT=1999
SEED=2026101067
EXPECTED={
    "discovery":(166,0.020529254583812922),
    "validation":(181,0.018672971642749295),
    "third":(204,0.01468491968437185),
}

def deterministic_seed(*parts):
    payload="|".join(map(str,(SEED,*parts)))
    return int.from_bytes(hashlib.sha256(payload.encode()).digest()[:8],"little")

def maximum_single_photo_edge_share(g, radius_km=50.):
    """Fraction of *local pair edges* incident to the most-connected photo."""
    loc=g.dropna(subset=["latitude","longitude"])
    d=base.pairwise_geo_km(loc.latitude.to_numpy(float),loc.longitude.to_numpy(float))
    u,v=np.triu_indices(len(loc),k=1)
    keep=d[u,v]<=radius_km
    m=int(keep.sum())
    if not m:return np.nan
    degree=np.bincount(np.r_[u[keep],v[keep]],minlength=len(loc))
    return float(degree.max()/m)

def summarize_threshold(rows:pd.DataFrame, nulls:list[np.ndarray], threshold:int,cohort:str)->dict:
    select=(rows.local_pairs_50km.to_numpy(int)>=threshold) & np.isfinite(rows.depletion_50km.to_numpy(float))
    n=int(select.sum())
    result={"min_original_local_pairs_per_species":int(threshold),"n_species":n,
            "share_of_primary_30_pair_eligible_species":None,
            "status":"HOLD_INSUFFICIENT_SPECIES" if n<30 else "POSTHOC_ESTIMABLE"}
    if n==0:return result
    data=rows.loc[select]
    dep=data.depletion_50km.to_numpy(float)
    values=np.vstack([x for x,t in zip(nulls,select) if t])
    if values.shape!=(n,base.PERMUTATIONS):
        raise RuntimeError("per-species null vector/source group mismatch")
    means=values.mean(axis=0)
    result.update({
        "mean_depletion":float(dep.mean()),
        "median_depletion":float(np.median(dep)),
        "positive_species_fraction":float(np.mean(dep>0)),
        "mean_specieswide_pair_discordance":float(data.D_pair.mean()),
        "mean_local_pair_discordance":float(data.D_local_50km.mean()),
        "relative_depletion_fraction":float(dep.mean()/data.D_pair.mean()),
        "median_original_50km_pairs_per_species":float(data.local_pairs_50km.median()),
        "max_single_photo_local_edge_share_median":float(data.max_single_photo_edge_share.median()),
        "max_single_photo_local_edge_share_q90":float(data.max_single_photo_edge_share.quantile(.90)),
        "null_mean_depletion":float(means.mean()),
        "permutation_p_upper":float((1+(means>=dep.mean()).sum())/(base.PERMUTATIONS+1)),
        "null_95pct":[float(v) for v in np.quantile(means,[.025,.975])],
    })
    if n>=30:
        rng=np.random.default_rng(deterministic_seed(cohort,threshold,"bootstrap"))
        b=np.mean(dep[rng.integers(n,size=(N_BOOT,n))],axis=1)
        result["species_bootstrap_ci95"]=[float(x) for x in np.quantile(b,[.025,.975])]
    else:
        result["species_bootstrap_ci95"]=None
    result["supported_descriptive"] = bool(
        n>=30 and dep.mean()>0 and result["permutation_p_upper"]<.05)
    return result

def run_cohort(path:Path, cohort:str)->tuple[dict,pd.DataFrame]:
    # SHA and pre-existing n>=40-classifiable rule are verified by prior loader.
    source=base.load(path,cohort)
    rows=[];nulls=[]
    for _,g in source.groupby("inat_taxon_id",sort=True):
        row,nul=base.analyse_species(g,cohort)
        row["max_single_photo_edge_share"]=maximum_single_photo_edge_share(g)
        rows.append(row)
        nulls.append(nul[50.])
    data=pd.DataFrame(rows).sort_values("inat_taxon_id",kind="stable").reset_index(drop=True)
    primary=summarize_threshold(data,nulls,30,cohort)
    exp_n,exp_dep=EXPECTED[cohort]
    if primary["n_species"]!=exp_n or abs(primary["mean_depletion"]-exp_dep)>1e-10:
        raise RuntimeError(f"Source frozen 50-km baseline mismatch in {cohort}: {primary}")
    ladder={str(th):summarize_threshold(data,nulls,th,cohort) for th in THRESHOLDS}
    for th in ladder:
        ladder[th]["share_of_primary_30_pair_eligible_species"] = (
            ladder[th]["n_species"]/exp_n)
    return {
        "n_eligible_colour_species":len(data),
        "n_original_50km_baseline_species":exp_n,
        "original_50km_baseline_depletion":exp_dep,
        "primary_stress_threshold_pairs":PRIMARY_STRESS_THRESHOLD,
        "ladder":ladder,
    },data

def main():
    ap=argparse.ArgumentParser()
    for c in EXPECTED:
        ap.add_argument("--"+c,type=Path,required=True)
    ap.add_argument("--outdir",type=Path,required=True)
    a=ap.parse_args()
    d={"schema":"fcp_50km_original_pair_support_sensitivity_20261010_v1",
       "role":"post_outcome_observation_support_falsification",
       "status":"COMPLETE_SOURCE_FROZEN_POSTHOC",
       "confirmatory_decisions_changed":False,
       "source_sha256":base.SHA,"source_null_permutations":base.PERMUTATIONS,
       "radius_km":50,"minimum_local_pair_thresholds":list(THRESHOLDS),
       "independent_photo_species_biological_population_equivalence":False,
       "hard_nonclaims":[
           "The 50-km photo-pair count is sample geometry/opportunity, not population census.",
           "This post hoc ladder does not upgrade original post hoc inference to pre-registered confirmation.",
           "More pairs need not imply independent photo observations; network edges share endpoints.",
           "Spatially observed flower colours are not validated genotypes, fitness, or selection effects.",
           "Threshold-selected cohorts differ; effect changes cannot be interpreted as causal sample-size effects.",
       ]}
    species=[]
    for cohort in EXPECTED:
        d[cohort],details=run_cohort(getattr(a,cohort),cohort)
        details.insert(0,"cohort",cohort)
        species.append(details)
    d["replication"]={
       str(th):{
           "positive_statistic_all_three":all(d[c]["ladder"][str(th)].get("mean_depletion",0)>0 for c in EXPECTED),
           "meets_30_species_each":all(d[c]["ladder"][str(th)]["n_species"]>=30 for c in EXPECTED),
           "supported_permutation_all_three_given_30species":all(
               d[c]["ladder"][str(th)].get("supported_descriptive",False) for c in EXPECTED)
       } for th in THRESHOLDS}
    out=a.outdir;out.mkdir(parents=True,exist_ok=True)
    pd.concat(species,ignore_index=True).to_csv(out/"all_original_species_50km_pair_support.csv",index=False)
    (out/"result.json").write_text(json.dumps(d,indent=2)+"\n",encoding="utf8")
    print(json.dumps({c:d[c]["ladder"] for c in EXPECTED},indent=2))
    return d
if __name__=="__main__": main()
