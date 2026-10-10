#!/usr/bin/env python3
"""Constructive, composition-preserving adversarial colour-label sensitivity for FCP.

NOT a measured misclassification rate or proof that measurement errors exist.
"""
from __future__ import annotations

import argparse
import hashlib
import heapq
import json
from pathlib import Path
import numpy as np
import pandas as pd

SHA = {
    "discovery": "ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4",
    "validation": "0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6",
    "third": "57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186",
}
STATES = ["white", "yellow_orange", "red_pink", "blue_purple"]
EXPECTED_N = {"discovery": 166, "validation": 181, "third": 204}
EXPECTED_DEPLETION = {
    "discovery": 0.020529254583812922,
    "validation": 0.018672971642749295,
    "third": 0.01468491968437185,
}
EARTH_RADIUS_KM = 6371.0088
RADIUS_KM = 50.
MIN_LOCAL_PAIRS = 30
BUDGET_PCT = [0., 0.1, 0.25, 0.5, 1., 2., 5., 10.]
MAX_BUDGET_PCT = 10.


def checksum(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for part in iter(lambda: fh.read(1 << 20), b""):
            h.update(part)
    return h.hexdigest()


def as_bool(series):
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).astype(bool)
    return series.fillna("").astype(str).str.strip().str.casefold().isin(
        ["true", "1", "yes", "y"])


def pair_km(lat, lon):
    la = np.deg2rad(np.asarray(lat, float))
    lo = np.deg2rad(np.asarray(lon, float))
    xyz = np.column_stack([np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)])
    return EARTH_RADIUS_KM*np.arccos(np.clip(xyz@xyz.T, -1, 1))


def load(path, cohort):
    if checksum(path) != SHA[cohort]:
        raise RuntimeError(f"{cohort}: frozen source SHA mismatch")
    frame = pd.read_csv(path, low_memory=False)
    required = {"inat_taxon_id", "species", "photo_id", "latitude",
                "longitude", "morph", "global_classifiable"}
    if required-set(frame):
        raise RuntimeError(f"missing source columns {sorted(required-set(frame))}")
    keep = as_bool(frame.global_classifiable)&frame.morph.astype(str).isin(STATES)
    frame = frame.loc[keep].copy()
    counts = frame.groupby("inat_taxon_id").size()
    frame = frame.loc[frame.inat_taxon_id.isin(counts[counts>=40].index)]
    if frame.photo_id.duplicated().any():
        raise RuntimeError("duplicate immutable source photo identifiers")
    return frame.sort_values(["inat_taxon_id", "photo_id"], kind="stable")


class Species:
    def __init__(self, frame):
        g = frame.dropna(subset=["latitude", "longitude"]).reset_index(drop=True)
        if len(g) < 40:
            raise RuntimeError("existing eligible species lost coordinates")
        self.taxon = int(g.inat_taxon_id.iloc[0])
        self.name = str(g.species.iloc[0])
        self.photo = g.photo_id.astype(str).to_numpy()
        self.orig = pd.Categorical(g.morph.astype(str), categories=STATES).codes.astype(np.int8)
        if (self.orig < 0).any():
            raise RuntimeError("nonbiological colour label")
        self.labels = self.orig.copy()
        self.n = len(g)
        a,b = np.triu_indices(self.n,1)
        self.all_a, self.all_b = a,b
        near = pair_km(g.latitude.to_numpy(), g.longitude.to_numpy())[a,b]<=RADIUS_KM
        self.edge_a, self.edge_b = a[near],b[near]
        self.m = len(self.edge_a)
        self.touched = np.zeros(self.n,dtype=bool)
        self.adj = np.zeros((self.n,self.n), dtype=np.int8)
        self.adj[self.edge_a,self.edge_b]=1
        self.adj[self.edge_b,self.edge_a]=1
        self.fixed_global = float(np.mean(self.orig[a]!=self.orig[b]))
        self.local = float(np.mean(self.orig[self.edge_a]!=self.orig[self.edge_b])) if self.m else 0.
        self.num_swaps = 0

    @property
    def depletion(self):
        return self.fixed_global-self.local

    def candidate(self):
        """Exact change in LOCAL mismatch count for a swap, not an approximation.

        For i with colour a and j colour b, swapping labels changes the local
        mismatch count by C_i(a)-C_i(b)+C_j(b)-C_j(a). i-j edge is unchanged.
        """
        if self.m < MIN_LOCAL_PAIRS:
            return None
        valid = (~self.touched)
        if valid.sum()<2:
            return None
        labels=self.labels.astype(int)
        count = self.adj.astype(np.int32)@np.eye(4,dtype=np.int32)[labels]
        indices=np.arange(self.n)
        ci_self=count[indices, labels]
        score = (ci_self[:,None]-count[:,labels]) + (
            ci_self[None,:]-count[:,labels].T)
        allowed = np.triu(np.ones((self.n,self.n), dtype=bool),1)
        allowed &= valid[:,None]&valid[None,:]&(labels[:,None]!=labels[None,:])
        score=np.where(allowed,score,-10**9)
        flat=np.argmax(score)
        delta=int(score.ravel()[flat])
        if delta<=0:
            return None
        i,j=map(int,np.unravel_index(flat,score.shape))
        return (delta/self.m, delta, i,j)

    def swap(self,i,j,expected_delta):
        assert not self.touched[i] and not self.touched[j] and self.labels[i]!=self.labels[j]
        before=float(self.local)
        aa,bb=int(self.labels[i]),int(self.labels[j])
        self.labels[i],self.labels[j]=self.labels[j],self.labels[i]
        self.touched[[i,j]]=True
        self.num_swaps+=1
        self.local=float(np.mean(self.labels[self.edge_a]!=self.labels[self.edge_b]))
        if not np.isclose(self.local-before,expected_delta/self.m,atol=1e-12):
            raise RuntimeError("exact colour-swap delta formula failed")
        if not np.array_equal(np.bincount(self.orig,minlength=4),
                              np.bincount(self.labels,minlength=4)):
            raise RuntimeError("within-species original colour-state counts changed")
        return {"inat_taxon_id":self.taxon,"species":self.name,
                "photo_id_i":self.photo[i], "photo_id_j":self.photo[j],
                "original_label_i":STATES[aa],"original_label_j":STATES[bb],
                "new_label_i":STATES[bb],"new_label_j":STATES[aa],
                "local_discordant_pair_increment":expected_delta,
                "n_local_pairs":self.m}


def build_species(frame):
    s=[]
    for _,g in frame.groupby("inat_taxon_id",sort=True):
        x=Species(g)
        if x.m >= MIN_LOCAL_PAIRS:
            s.append(x)
    return s


def analyse(cohort, frame):
    species=build_species(frame)
    if len(species)!=EXPECTED_N[cohort]:
        raise RuntimeError(f"{cohort}: eligible 50km colour species count mismatch {len(species)}")
    original_mean=float(np.mean([s.depletion for s in species]))
    if abs(original_mean-EXPECTED_DEPLETION[cohort])>1e-10:
        raise RuntimeError(f"{cohort}: original 50km source local depletion changed: {original_mean}")
    denominator=sum(s.n for s in species)
    max_swaps=int(np.floor(MAX_BUDGET_PCT/100*denominator/2))
    heap=[]
    for ix,s in enumerate(species):
        best=s.candidate()
        if best is not None:
            gain, delta,i,j=best
            heapq.heappush(heap,(-gain,ix,delta,i,j))
    trajectory=[{"edited_photos":0,"fraction_edited_pct":0.,"mean_depletion":original_mean}]
    changes=[]
    total_depletion=sum(s.depletion for s in species)
    first_nonpositive=None
    while heap and len(changes)<max_swaps and total_depletion>0:
        neg_gain,idx,delta,i,j=heapq.heappop(heap)
        part=species[idx]
        change=part.swap(i,j,delta)
        total_depletion-=(-neg_gain)
        mean=float(np.mean([s.depletion for s in species]))
        if not np.isclose(mean,total_depletion/len(species),atol=1e-11):
            raise RuntimeError("global monotone sum inconsistency")
        rec={"edited_photos":2*(len(changes)+1),
             "fraction_edited_pct":200*(len(changes)+1)/denominator,
             "mean_depletion":mean}
        trajectory.append(rec)
        change["cohort"]=cohort
        change["mean_depletion_after_swap"]=mean
        changes.append(change)
        if mean<=0:
            first_nonpositive=rec
            break
        candidate=part.candidate()
        if candidate is not None:
            gain,new_delta,new_i,new_j=candidate
            heapq.heappush(heap,(-gain,idx,new_delta,new_i,new_j))
    budgets={}
    for target in BUDGET_PCT:
        eligible=[r for r in trajectory if r["fraction_edited_pct"]<=target+1e-12]
        budgets[str(target)] = {"budget_pct":target,"attained":eligible[-1],
                                "n_swaps":eligible[-1]["edited_photos"]//2}
    result={
        "n_species":len(species),"n_classifiable_photos_in_50km_evaluable_species":denominator,
        "original_depletion":original_mean,
        "first_constructive_nonpositive":first_nonpositive,
        "no_nonpositive_within_maximum_budget":first_nonpositive is None,
        "last_evaluated":trajectory[-1],
        "fraction_edited_denominator":"original classifiable photos belonging to species with >=30 local 50km pairs",
        "all_swaps_are_within_species_and_involve_two_distinct_previously_unedited_photos":True,
        "species_wide_four_colour_frequencies_preserved_exactly":True,
        "budget_curve":budgets,
        "number_of_accepted_swaps":len(changes),
        "optimality":"Greedy constructive feasible attack, NOT mathematically minimal error rate.",
        "inference":"This is an adversarial label-swap fragility bound, NOT calibrated image misclassification."
    }
    return result, pd.DataFrame(changes),pd.DataFrame(trajectory)


def run(paths,out):
    final={"schema":"fcp_adversarial_morph_labels_20261010_v1",
          "date_jst":"2026-10-10","role":"posthoc_adversarial_colour_specific_error_sensitivity",
          "primary_km":RADIUS_KM,"budget_percentages":BUDGET_PCT,
          "source_sha256":SHA,"confirmatory_decisions_changed":False,
          "hard_nonclaims":[
              "Constructed swaps are not empirically observed photo-label errors.",
              "Attack fraction is not an estimated error prevalence, confidence interval or optimal minimum.",
              "A greedy worst-case attack does not estimate random or image-realistic error susceptibility.",
              "No genetic FCP, fitness, selection, plasticity, or independently validated morphology."]}
    edits=[];curves=[]
    for cohort,p in paths.items():
        r,rows,curve=analyse(cohort,load(p,cohort))
        final[cohort]=r
        edits.append(rows);curves.append(curve.assign(cohort=cohort))
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    pd.concat(edits,ignore_index=True).to_csv(out/"constructed_label_swaps.csv",index=False)
    pd.concat(curves,ignore_index=True).to_csv(out/"greedy_sensitivity_trajectory.csv",index=False)
    (out/"result.json").write_text(json.dumps(final,indent=2)+"\n")
    return final


def main():
    parser=argparse.ArgumentParser()
    for c in SHA: parser.add_argument("--"+c,required=True,type=Path)
    parser.add_argument("--outdir",type=Path,required=True)
    a=parser.parse_args()
    result=run({c:getattr(a,c) for c in SHA},a.outdir)
    print(json.dumps({c:result[c] for c in SHA},indent=2))


if __name__=="__main__":
    main()
