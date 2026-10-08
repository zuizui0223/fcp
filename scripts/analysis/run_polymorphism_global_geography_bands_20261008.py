#!/usr/bin/env python3
"""Species-balanced FCP geographic atlas: visible colour across latitude/elevation bands.

All inferences are WITHIN-SPECIES comparisons among photographic regions,
not prevalence estimates for all flowering plants. Each species supplies
the same weight in every supported band. Geography and exact species-wide
morph counts stay fixed in the within-species label permutation null.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

import run_polymorphism_white_chromatic_clines_20261008 as base

MORPHS = base.MORPHS
PHOTO_SHA256 = base.INPUT_SHA256
EDGES = {
    "absolute_latitude": [0, 15, 30, 45, 60, np.inf],
    "elevation": [-np.inf, 250, 1000, 2000, np.inf],
}
NAMES = {
    "absolute_latitude": ["0–15°", "15–30°", "30–45°", "45–60°", "60–90°"],
    "elevation": ["<250 m", "250–1000 m", "1000–2000 m", "≥2000 m"],
}
MIN_SPECIES_PHOTOS = 40
MIN_PHOTOS_IN_BAND = 8
MIN_PHOTOS_OUTSIDE_BAND = 8
MIN_INFORMATIVE_SPECIES_PER_BAND = 25
PERMUTATIONS = 199
BOOTSTRAPS = 999
MASTER_SEED = 2026100851


def rng_seed(*args: object) -> int:
    val = "|".join(map(str, (MASTER_SEED, *args)))
    return int.from_bytes(hashlib.sha256(val.encode()).digest()[:8], "little")


def pair_diversity(counts: np.ndarray) -> float:
    """Unbiased pairwise probability of distinct photo colour labels."""
    n = int(np.sum(counts))
    if n < 2:
        return float("nan")
    return float((n*n-np.sum(np.asarray(counts, float)**2))/(n*(n-1)))


def region_assignments(x: np.ndarray, axis: str) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    bins = np.full(len(x), -1, dtype=np.int16)
    good = np.isfinite(x)
    if axis == "absolute_latitude":
        good &= (x>=0) & (x<=90)
    else:
        good &= (x>=-500) & (x<=9000)
    if good.any():
        b = np.searchsorted(np.asarray(EDGES[axis]), x[good], side="right")-1
        if np.any((b<0)|(b>=len(NAMES[axis]))):
            raise ValueError("out-of-range geographic band assignment")
        bins[good] = b
    return bins


def species_band_metrics(
    values: np.ndarray, morph: np.ndarray, axis: str, *,
    nperm: int, seed: int
) -> tuple[list[dict], np.ndarray, np.ndarray]:
    values = np.asarray(values, dtype=float)
    morph = np.asarray(morph, dtype=str)
    if len(values)!=len(morph):
        raise ValueError("axis and flower-label lengths differ")
    b = region_assignments(values,axis)
    usable = b>=0
    k=len(NAMES[axis])
    null_sum = np.zeros((k,nperm,2),dtype=float)
    null_mixed = np.zeros((k,nperm),dtype=float)
    if int(usable.sum())<MIN_SPECIES_PHOTOS:
        return [],null_sum,null_mixed
    b=b[usable]
    lab=pd.Categorical(morph[usable],categories=MORPHS).codes
    if np.any(lab<0):
        raise ValueError("unknown flower morph code")
    n=len(b)
    all_counts=np.bincount(lab,minlength=4)
    global_white=float(all_counts[0]/n)
    global_D=pair_diversity(all_counts)
    counts=np.bincount(b,minlength=k)
    qualified = (counts>=MIN_PHOTOS_IN_BAND)&(n-counts>=MIN_PHOTOS_OUTSIDE_BAND)
    if not qualified.any():
        return [],null_sum,null_mixed
    tab=np.bincount(4*b+lab,minlength=4*k).reshape(k,4)
    result=[]
    for i in np.flatnonzero(qualified):
        c=tab[i]
        white_fraction=float(c[0]/counts[i])
        observed_D=pair_diversity(c)
        result.append({
            "band_index":int(i),"n_photos_in_band":int(counts[i]),
            "n_photos_outside":int(n-counts[i]),
            "species_global_visible_white_fraction":global_white,
            "band_visible_white_fraction":white_fraction,
            "delta_white_vs_same_species":white_fraction-global_white,
            "band_four_state_pair_diversity":observed_D,
            "species_global_pair_diversity":global_D,
            "delta_diversity_vs_same_species":observed_D-global_D,
            "white_chromatic_photographic_mixture":bool(c[0]>=3 and int(c[1:].sum())>=3),
            "nonwhite_hue_mixture":bool(int((c[1:]>=3).sum())>=2),
        })
    rng=np.random.default_rng(seed)
    for j in range(nperm):
        perm=rng.permutation(lab)
        shuffled=np.bincount(4*b+perm,minlength=4*k).reshape(k,4)
        for i in np.flatnonzero(qualified):
            c=shuffled[i]
            null_sum[i,j,0]=(c[0]/counts[i])-global_white
            null_sum[i,j,1]=pair_diversity(c)-global_D
            null_mixed[i,j]=int(c[0]>=3 and int(c[1:].sum())>=3)
    return result,null_sum,null_mixed


def summarize_band(rows: list[dict], null_values: np.ndarray, mixed_null: np.ndarray,
                   cohort: str, axis: str, idx: int) -> dict:
    n=len(rows)
    if n==0:
        return {"n_informative_species":0,"status":"NO_WITHIN_SPECIES_BAND_OPPORTUNITY"}
    df=pd.DataFrame(rows)
    mat=np.asarray(null_values,float)
    if mat.shape!=(PERMUTATIONS,2):
        raise ValueError("permutation result shape mismatch")
    rng=np.random.default_rng(rng_seed(cohort,axis,idx,"bootstrap"))
    effects=df[["delta_white_vs_same_species","delta_diversity_vs_same_species"]].to_numpy(float)
    obs=effects.mean(axis=0)
    boot=effects[rng.integers(0,n,(BOOTSTRAPS,n))].mean(axis=1)
    null_sd=mat.std(axis=0,ddof=1)
    groups=df.groupby("genus")[["delta_white_vs_same_species","delta_diversity_vs_same_species"]].mean()
    return {
        "n_informative_species":int(n),
        "n_genera":int(len(groups)),
        "n_classifiable_photo_records":int(df.n_photos_in_band.sum()),
        "n_observed_photo_white_colour_mixed_species":int(df.white_chromatic_photographic_mixture.sum()),
        "fraction_observed_white_colour_mixed_species":float(df.white_chromatic_photographic_mixture.mean()),
        "null_expected_mixed_species":float(np.asarray(mixed_null).mean()),
        "mean_species_equal_white_fraction_in_band":float(df.band_visible_white_fraction.mean()),
        "mean_species_equal_fourstate_D_in_band":float(df.band_four_state_pair_diversity.mean()),
        "mean_white_deviation_vs_specieswide":float(obs[0]),
        "mean_D_deviation_vs_specieswide":float(obs[1]),
        "mean_genus_balanced_white_deviation":float(groups.iloc[:,0].mean()),
        "mean_genus_balanced_D_deviation":float(groups.iloc[:,1].mean()),
        "bootstrap_95CI_white_deviation":[float(t) for t in np.quantile(boot[:,0],[.025,.975])],
        "bootstrap_95CI_D_deviation":[float(t) for t in np.quantile(boot[:,1],[.025,.975])],
        "null_white_sd":float(null_sd[0]),
        "null_D_sd":float(null_sd[1]),
        "permutation_p_two_sided_unadjusted":[float(
            (1+np.count_nonzero(np.abs(mat[:,j]-mat[:,j].mean()) >= abs(obs[j]-mat[:,j].mean())))
            /(PERMUTATIONS+1)) for j in range(2)],
        "minimum_25_species_pass":bool(n>=MIN_INFORMATIVE_SPECIES_PER_BAND),
        "status":("COVERAGE_ELIGIBLE_EXPLORATORY" if n>=MIN_INFORMATIVE_SPECIES_PER_BAND
                  else "HOLD_SPARSE_WITHIN_SPECIES_BAND_OPPORTUNITY"),
    }


def hemisphere_coverage(g: pd.DataFrame) -> list[dict]:
    """Descriptive coverage only; the hemisphere species pools are different."""
    g=g.copy()
    g["bin"]=region_assignments(np.abs(g.latitude.to_numpy(float)),"absolute_latitude")
    g["hemisphere"]=np.where(g.latitude>=0,"N","S")
    records=[]
    for (hemi,band), grp in g.groupby(["hemisphere","bin"],sort=True):
        if int(band)<0: continue
        selected=grp.groupby("inat_taxon_id").filter(lambda x:len(x)>=MIN_PHOTOS_IN_BAND)
        if selected.empty: continue
        per=selected.groupby("inat_taxon_id").morph.apply(
            lambda a: float((a=="white").mean()))
        records.append({"hemisphere":str(hemi),"band":NAMES["absolute_latitude"][int(band)],
                        "n_species_with_at_least_8_photos":int(len(per)),
                        "n_photos_in_supported_species":int(len(selected)),
                        "mean_species_equal_photo_white_fraction":float(per.mean()),
                        "only_descriptive_not_a_between_species_adaptation_test":True})
    return records


def analyse(d: pd.DataFrame,cohort: str,axis: str) -> tuple[dict,pd.DataFrame]:
    k=len(NAMES[axis])
    all_rows=[]
    null_sum=np.zeros((k,PERMUTATIONS,2),float)
    mixed_null=np.zeros((k,PERMUTATIONS),float)
    species_supported=np.zeros(k,int)
    photos_in_bin=np.zeros(k,int)
    for taxon, group in d.groupby("inat_taxon_id",sort=True):
        values=(np.abs(group.latitude.to_numpy(float)) if axis=="absolute_latitude"
                else group.elevation_m.to_numpy(float))
        photo_bins=region_assignments(values,axis)
        for i in range(k):
            if int((photo_bins==i).sum())>=MIN_PHOTOS_IN_BAND:
                species_supported[i]+=1
                photos_in_bin[i]+=int((photo_bins==i).sum())
        values_rows,nul,n_mix=species_band_metrics(
            values,group.morph.astype(str).to_numpy(),axis,
            nperm=PERMUTATIONS,seed=rng_seed(cohort,axis,int(taxon)))
        null_sum+=nul
        mixed_null+=n_mix
        for row in values_rows:
            row.update({"cohort":cohort,"species":str(group.species.iloc[0]),
                        "genus":str(group.species.iloc[0]).split()[0],
                        "inat_taxon_id":int(taxon),"axis":axis,
                        "band":NAMES[axis][row["band_index"]]})
            all_rows.append(row)
    records=pd.DataFrame(all_rows)
    by_band={}
    for i,name in enumerate(NAMES[axis]):
        sel=records.loc[records.band_index==i] if len(records) else []
        payload=sel.to_dict("records") if len(records) else []
        n=len(payload)
        nul=null_sum[i]/n if n else np.zeros((PERMUTATIONS,2))
        mixed=mixed_null[i]/n if n else np.zeros(PERMUTATIONS)
        val=summarize_band(payload,nul,mixed,cohort,axis,i)
        val["n_species_with_min8_photos_in_bin_even_without_outside_opportunity"]=int(species_supported[i])
        val["n_photos_in_species_with_min8_in_bin"]=int(photos_in_bin[i])
        by_band[name]=val
    # Multiplicity over all coverage-eligible latitude or elevation bins,
    # separately for white and full four-state diversity.
    for metric,col in (("white",0),("D",1)):
        ids=[i for i,n in enumerate(NAMES[axis]) if
             by_band[n]["minimum_25_species_pass"] and
             by_band[n][f"null_{metric}_sd"]>1e-12]
        if not ids: continue
        nulls=np.vstack([null_sum[i,:,col]/by_band[NAMES[axis][i]]["n_informative_species"]
                          for i in ids])
        mean=nulls.mean(axis=1)
        sd=nulls.std(axis=1,ddof=1)
        max_null=np.max(np.abs((nulls-mean[:,None])/sd[:,None]),axis=0)
        for j,i in enumerate(ids):
            name=NAMES[axis][i]
            obs=by_band[name][("mean_white_deviation_vs_specieswide" if col==0
                               else "mean_D_deviation_vs_specieswide")]
            z=abs((obs-mean[j])/sd[j])
            by_band[name][f"maxT_across_{axis}_bands_p_{metric}"]=float(
                (1+np.count_nonzero(max_null>=z))/(PERMUTATIONS+1))
    return {"n_source_40_classifiable_species":int(d.inat_taxon_id.nunique()),
            "axis":axis,"bands":by_band,
            "n_bands_eligible_with_25_species":sum(x["minimum_25_species_pass"] for x in by_band.values())},records


def main() -> None:
    p=argparse.ArgumentParser()
    for cohort in PHOTO_SHA256:
        p.add_argument("--"+cohort,required=True,type=Path)
    p.add_argument("--elevation",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    args=p.parse_args()
    args.outdir.mkdir(parents=True,exist_ok=True)
    out={
        "schema":"fcp_geographic_latitude_elevation_band_atlas_v1",
        "date_jst":"2026-10-08",
        "status":"complete",
        "role":"retrospective_species_equal_geographical_description_and_within_species_permutation_diagnostic",
        "confirmatory_decisions_changed":False,
        "source_measured_sha256":PHOTO_SHA256,
        "elevation_raster_sha256":base.file_sha(args.elevation),
        "latitude_band_edges_degrees":[0,15,30,45,60,90],
        "elevation_bands_m":["below 250","250–1000","1000–2000","2000+"],
        "minimum_photo_counts":{"species_global":MIN_SPECIES_PHOTOS,"in_band":MIN_PHOTOS_IN_BAND,
             "outside_band":MIN_PHOTOS_OUTSIDE_BAND,"band_informative_species":MIN_INFORMATIVE_SPECIES_PER_BAND},
        "sampling_design":"species-equal in each bin; requires repeated locations of SAME species, not species turnover across bins",
        "permutations":PERMUTATIONS,"source_system":"iNaturalist same classifier, three species-disjoint cohorts",
        "hard_nonclaims":[
            "the photographic species frame is not an unbiased world flora census",
            "the number of species in a latitude/elevation band is NOT the incidence of natural floral polymorphism",
            "visible white has demonstrable exposure coupling; cannot be read as biochemical anthocyanin absence",
            "WorldClim 10-minute elevation is not the individual plant's field elevation",
            "the hemisphere-specific display changes species composition and is descriptive only",
            "band-level deviations can reflect observer/site/phenology confounding or historical spatial structure",
            "genus balance is not a species phylogenetic covariance correction",
            "a band with fewer than 25 informative species is coverage limited, NOT evidence of no effect",
            "no genetic polymorphism, adaptive selection, pollinator preference, or causal high-altitude pigmentation claimed",
        ],
        "cohorts":{},
    }
    details=[]
    for cohort in PHOTO_SHA256:
        d=base.read_cohort(getattr(args,cohort),cohort)
        d["elevation_m"]=base.sample_elevation(args.elevation,d)
        a,dr=analyse(d,cohort,"absolute_latitude")
        e,er=analyse(d,cohort,"elevation")
        out["cohorts"][cohort]={"absolute_latitude":a,"elevation":e,
                                "hemisphere_coverage_descriptive":hemisphere_coverage(d)}
        details.extend([dr,er])
    pd.concat(details,ignore_index=True).to_csv(args.outdir/"species_geographic_band_estimates.csv",index=False)
    (args.outdir/"result.json").write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps(out,indent=2,ensure_ascii=False))


if __name__=="__main__":
    main()
