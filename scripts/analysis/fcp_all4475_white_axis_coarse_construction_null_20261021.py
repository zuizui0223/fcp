#!/usr/bin/env python3
"""Post hoc coarse-state- and source-cohort-preserving W null on all 4,475 taxa.

Identical frozen two-mode Hellinger geometry, nine-palette white-vs-rest axis,
full source photo composition. Unlike the 158-species prospectively frozen test,
this is an EXPLORATORY test on an already-open global photo source. Every
4,475 species has a W measurement; the original >=40, 10%/10% sensitivity
set remains a declared side analysis, not a replacement denominator.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from fcp_full4475_hypothesis_census_20261021 import (
    SOURCE_SHA,TOTAL_TAXA,CLASSIFIABLE_ROWS,PERMS,Q,
    PHOTO_LABELS,load_original,check,
)
from disttrait.geometry import two_mode_axis

SEED=20261021
COHORTS=("discovery","validation","prospective","wave4","new_U100")

def fit_axis(p):
    try:
        a=two_mode_axis(p).unit_axis
    except ValueError:
        return None
    return float((a @ Q)**2)

def prepare(d):
    d=d.sort_values(["inat_taxon_id","photo_id"],kind="mergesort").reset_index(drop=True)
    palette=np.stack(d.palette_norm9.to_numpy())
    species=d.inat_taxon_id.to_numpy(dtype=int)
    source=d.fcp_source_cohort.to_numpy(dtype=str)
    coarse=d.morph.to_numpy(dtype=str)
    tids,starts,counts=np.unique(species,return_index=True,return_counts=True)
    check(len(tids)==TOTAL_TAXA and counts.min()>=2,"less than 4475 source taxa")
    sl=[slice(int(start),int(start+count)) for start,count in zip(starts,counts)]
    idx={}
    for cohort in COHORTS:
        for state in PHOTO_LABELS:
            group=np.flatnonzero((source==cohort)&(coarse==state))
            if len(group):
                idx[(cohort,state)]=group
    check(sum(len(z) for z in idx.values())==CLASSIFIABLE_ROWS,
          "not all original classifiable rows in coarse-cohort preserved null")
    return tids,sl,palette,idx

def summarize_w(vals,groups,primary):
    states={"all_4475":np.ones(len(vals),bool),
            "frozen_primary_gate_on_all_4475_frame":primary}
    for c in COHORTS:states["cohort_"+c]=(groups==c)
    out={}
    for key,mask in states.items():
        v=vals[mask]
        out[key]={
          "n_species":int(len(v)),
          "observed_mean_W":float(v.mean()),
          "null_median_W":None,
          "null_q025_q975":None,
          "null_upper_p":None,
          "n_null_worlds":0,
        }
    return states,out

def synthetic():
    a=np.array([[.85,.15]+[0.]*7,[.9,.1]+[0.]*7,
                [.15,.85]+[0.]*7,[.1,.9]+[0.]*7])
    x=fit_axis(a)
    check(x is not None and 0<=x<=1,"synthetic W must be bounded")
    print("FCP_ALL4475_W_SOURCE_COMPOSITION_NULL_SYNTHETIC_PASS")

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--photos",type=Path)
    p.add_argument("--species-master",type=Path)
    p.add_argument("--out",type=Path)
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    synthetic()
    if a.self_test:return
    check(a.photos and a.species_master and a.out,"missing source")
    t,d=load_original(a.photos)
    tids,sl,pal,groups=prepare(d)
    master=pd.read_csv(a.species_master,low_memory=False).set_index("inat_taxon_id",verify_integrity=True)
    check(len(master)==TOTAL_TAXA and set(master.index)==set(tids.tolist()),"master species missing")
    masters=master.loc[tids]
    primary=masters.H1_original_primary_gate_pass.astype(bool).to_numpy()
    cohorts=masters.source_cohort.astype(str).to_numpy()
    observed=np.empty(len(tids))
    for i,s in enumerate(sl):
        z=fit_axis(pal[s])
        check(z is not None,"degenerate source nonzero displacement for "+str(tids[i]))
        observed[i]=z
    check(abs(float(observed.mean())-0.4788719110425215)<1e-9,
          "source W differs from unmodified authentic core full-population result")
    state,output=summarize_w(observed,cohorts,primary)
    nulls={k:np.empty(PERMS,float) for k in output}
    rng=np.random.default_rng(SEED)
    any_degeneracy=0
    for replicate in range(PERMS):
        perm=pal.copy()
        # A SINGLE WORLD reshuffles complete 9-colour vectors only among
        # observations of the SAME source cohort and SAME four-colour label.
        # Species, individual colour labels, and observer/photo membership are
        # unchanged. No species dropped from any world.
        for ix in groups.values():
            if len(ix)>1:
                perm[ix]=pal[ix[rng.permutation(len(ix))]]
        w=np.empty(len(sl),dtype=float)
        for i,s in enumerate(sl):
            v=fit_axis(perm[s])
            if v is None:
                any_degeneracy+=1
                raise RuntimeError("FCP_NULL_W_DEGENERATE_PER_SPECIES: no matched full-source W; "+
                                   str(tids[i])+" world "+str(replicate))
            w[i]=v
        for key,mask in state.items():
            nulls[key][replicate]=float(w[mask].mean())
        if (replicate+1)%25==0:
            print("FCP_ALL4475_WORLD_WHITE_AXIS_NULLS_COMPLETED",replicate+1,flush=True)
    for name,values in nulls.items():
        obs=output[name]["observed_mean_W"]
        output[name].update({
          "null_median_W":float(np.median(values)),
          "null_mean_W":float(values.mean()),
          "null_q025_q975":[float(x) for x in np.quantile(values,[.025,.975])],
          "null_upper_p":float((1+int(np.count_nonzero(values>=obs)))/(PERMS+1)),
          "n_null_worlds":PERMS,
          "fixed_species_in_every_world":True,
          "frozen_source_cohort_and_coarse_colour_count_in_every_world":True,
          "source_descendants_with_no_assigned_morph":0,
        })
    result={
      "schema":"fcp_all4475_source_cohort_coarse_preserving_white_axis_posthoc_null_v1",
      "source_original_photo_csv_sha256":SOURCE_SHA,
      "taxon_count":TOTAL_TAXA,
      "genuine_original_photos":447500,
      "original_classifiable_four_state_photos":CLASSIFIABLE_ROWS,
      "W_axis":"squared dot of each species' unlabeled Hellinger 2-means composition displacement with q=[white positive,8 nonwhite equal negative], normalized",
      "null":"within the five frozen source cohorts and four frozen coarse label classes, shuffle entire nine-palette photo vectors among species, refit same Hellinger 2-mode for all species in each world",
      "N_PERMUTATIONS":PERMS,
      "observed_photo_source_already_opened":True,
      "prospective_158_species_original_confirmatory_W_result_unchanged":True,
      "not_independent_technical_or_genetic_validation":True,
      "does_not_estimate_ancestral_white_or_molecular_loss_rate":True,
      "no_taxa_excluded_from_full_primary_W":True,
      "all_coarse_labels_preserved_exactly_per_taxon":True,
      "H1_W":output,
    }
    a.out.mkdir(parents=True,exist_ok=True)
    (a.out/"all4475_white_axis_coarse_preserving_null_result.json").write_text(
       json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("FCP_ALL4475_FULL_CONSTRUCTION_NULL_COMPLETE",json.dumps(result,sort_keys=True,allow_nan=False))

if __name__=="__main__":
    main()
