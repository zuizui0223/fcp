#!/usr/bin/env python3
"""Post hoc evidence screen: do photographed floral morphs recur locally across years?

This counts geographically proximate multi-year observations, NOT biologically
demonstrated mating-population polymorphism, balancing selection, or fitness.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest

MORPHS = ("white", "yellow_orange", "red_pink", "blue_purple")
NONWHITE = MORPHS[1:]
COHORT_SHA256 = {
    "discovery": "ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4",
    "validation": "0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6",
    "third": "57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186",
}
DIAMETERS_KM = (10.0, 25.0, 50.0)
# Each anchor neighbourhood has radius = half of the reported maximum
# pairwise diameter; avoids calling two points up to 100 km apart "50 km local".
MIN_CLASSIFIABLE = 40
MIN_GLOBAL_MORPH = 5
MIN_ANCHOR_MORPH = 3
MIN_DISTINCT_YEARS = 2
MIN_DISTINCT_OBSERVERS = 2
MIN_TEMPORAL_PAIRS_EACH = 10
DISTANCE_BINS = (0.0, 10.0, 25.0, 50.0)
EARTH_RADIUS_KM = 6371.0088
SEED = 2026100813
N_BOOT = 1999


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for part in iter(lambda: f.read(1 << 20), b""):
            h.update(part)
    return h.hexdigest()


def as_bool(v: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(v):
        return v.fillna(False).astype(bool)
    return v.fillna("").astype(str).str.casefold().str.strip().isin({"true", "1", "yes", "y"})


def load(path: Path, cohort: str, *, enforce_hash: bool = True) -> pd.DataFrame:
    if enforce_hash and sha256(path) != COHORT_SHA256[cohort]:
        raise RuntimeError(f"{cohort}: measured table checksum mismatch")
    d = pd.read_csv(path, low_memory=False)
    cols = ["inat_taxon_id", "species", "photo_id", "latitude", "longitude",
            "observer_id", "observed_on", "global_classifiable", "morph"]
    missing = sorted(set(cols)-set(d.columns))
    if missing:
        raise RuntimeError(f"{cohort}: missing columns {missing}")
    d = d.loc[as_bool(d["global_classifiable"]) &
              d["morph"].astype(str).isin(MORPHS), cols].copy()
    d["latitude"] = pd.to_numeric(d.latitude, errors="coerce")
    d["longitude"] = pd.to_numeric(d.longitude, errors="coerce")
    d = d.loc[d.latitude.between(-90,90) & d.longitude.between(-180,180)].copy()
    count = d.groupby("inat_taxon_id").size()
    d = d.loc[d.inat_taxon_id.isin(count[count>=MIN_CLASSIFIABLE].index)].copy()
    d["year"] = pd.to_datetime(d.observed_on,errors="coerce").dt.year
    d.loc[~d.year.between(1990,2026), "year"] = np.nan
    d["observer"] = d.observer_id.fillna("").astype(str).str.strip()
    d.loc[d.observer.isin(("nan", "None", "<NA>")), "observer"] = ""
    if d.photo_id.duplicated().any():
        raise RuntimeError(f"{cohort}: duplicated selected photo IDs")
    d["cohort"] = cohort
    return d.sort_values(["inat_taxon_id","photo_id"],kind="stable").reset_index(drop=True)


def pair_distance_km(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    lat=np.deg2rad(np.asarray(lat,float))
    lon=np.deg2rad(np.asarray(lon,float))
    c=np.cos(lat)
    xyz=np.stack((c*np.cos(lon),c*np.sin(lon),np.sin(lat)),axis=1)
    return np.arccos(np.clip(xyz@xyz.T,-1,1))*EARTH_RADIUS_KM


def supports_repeated_morph(mask: np.ndarray, morph: np.ndarray,
                            year: np.ndarray, obs: np.ndarray,
                            category: str) -> bool:
    chosen = mask & (morph == category) & np.isfinite(year)
    if int(chosen.sum()) < MIN_ANCHOR_MORPH:
        return False
    if len(np.unique(year[chosen])) < MIN_DISTINCT_YEARS:
        return False
    valid_obs = obs[chosen]
    valid_obs = valid_obs[valid_obs != ""]
    return len(np.unique(valid_obs)) >= MIN_DISTINCT_OBSERVERS


def strict_anchor_neighborhoods(g: pd.DataFrame, diameter: float) -> dict:
    """Find a bounded-diameter anchor neighborhood (not a field population)."""
    dist = pair_distance_km(g.latitude.to_numpy(float),g.longitude.to_numpy(float))
    groups = g.morph.astype(str).to_numpy()
    years = g.year.to_numpy(float)
    observers = g.observer.astype(str).to_numpy()
    evidence_white_nonwhite = False
    evidence_two_hues = False
    any_white_nonwhite = False
    any_two_hues = False
    any_temporal_white_nonwhite = False
    n_anchors_white_nonwhite = 0
    n_anchors_two_hues = 0
    best_wc_observers=0
    for a in range(len(g)):
        # Restriction is conservative: all included photo coordinates within
        # diameter/2 of this anchor, hence mutual separation <= diameter.
        near = dist[a] <= diameter/2
        if int(near.sum()) < 2*MIN_ANCHOR_MORPH:
            continue
        white = near & (groups == "white")
        coloured = near & (groups != "white")
        count_white, count_coloured = int(white.sum()), int(coloured.sum())
        if count_white >= MIN_ANCHOR_MORPH and count_coloured >= MIN_ANCHOR_MORPH:
            any_white_nonwhite = True
            # Evidence is observed on >=2 different years and by >=2 different
            # observers PER colour side, not merely in aggregate.
            good_w = supports_repeated_morph(near,groups,years,observers,"white")
            colour_years = np.unique(years[coloured & np.isfinite(years)])
            colour_observers = np.unique(observers[coloured & (observers != "") &
                                                   np.isfinite(years)])
            good_c = (len(colour_years) >= MIN_DISTINCT_YEARS and
                      len(colour_observers) >= MIN_DISTINCT_OBSERVERS and
                      int(np.sum(coloured & np.isfinite(years))) >= MIN_ANCHOR_MORPH)
            evidence_white_nonwhite |= good_w and good_c
            n_anchors_white_nonwhite += int(good_w and good_c)
            best_wc_observers = max(best_wc_observers,
                                    min(len(np.unique(observers[white & (observers!="")])),
                                        len(colour_observers)))
            if len(colour_years) >= 2 and len(np.unique(years[white & np.isfinite(years)]))>=2:
                any_temporal_white_nonwhite = True
        good_hues = [hue for hue in NONWHITE
                     if supports_repeated_morph(near,groups,years,observers,hue)]
        evidence_two_hues |= len(good_hues) >= 2
        n_anchors_two_hues += int(len(good_hues) >= 2)
        any_two_hues |= (sum(np.sum(near & (groups==hue))>=MIN_ANCHOR_MORPH
                            for hue in NONWHITE) >= 2)
    return {
        "n_white_nonwhite_anchors": n_anchors_white_nonwhite,
        "n_nonwhite_hue_anchors": n_anchors_two_hues,
        "white_nonwhite_same_neighborhood": bool(any_white_nonwhite),
        "two_nonwhite_hues_same_neighborhood": bool(any_two_hues),
        "white_nonwhite_multi_year_two_observers": bool(evidence_white_nonwhite),
        "two_nonwhite_hues_multi_year_two_observers": bool(evidence_two_hues),
        "white_nonwhite_two_years_not_observer_strict": bool(any_temporal_white_nonwhite),
        "minimum_distinct_observers_each_morph_in_best_WC_neighborhood": int(best_wc_observers)
    }


def temporal_colour_discordance(g: pd.DataFrame) -> dict:
    """Distance-matched temporal *description*, not a causal time-selection test."""
    d = pair_distance_km(g.latitude.to_numpy(float),g.longitude.to_numpy(float))
    w = (g.morph.to_numpy() == "white")
    y = g.year.to_numpy(float)
    o = g.observer.to_numpy(str)
    n=len(g)
    ii,jj=np.triu_indices(n,k=1)
    usable=(np.isfinite(y[ii]) & np.isfinite(y[jj]) & (o[ii]!="") & (o[jj]!="")
            & (o[ii]!=o[jj]) & (d[ii,jj]<=DISTANCE_BINS[-1]))
    ii,jj=ii[usable],jj[usable]
    distance=d[ii,jj]
    diff_colour=(w[ii]!=w[jj])
    diff_year=(y[ii]!=y[jj])
    # No zeros-as-data when an opportunity cell lacks pairs.
    vals=[]
    for lo,hi in zip(DISTANCE_BINS[:-1],DISTANCE_BINS[1:]):
        binmask=(distance>=lo) & (distance < hi)
        a=binmask & diff_year
        b=binmask & ~diff_year
        if int(a.sum())==0 or int(b.sum())==0:
            continue
        vals.append((int(a.sum()),int(b.sum()),
                     float(np.mean(diff_colour[a])),
                     float(np.mean(diff_colour[b]))))
    cross_total=sum(v[0] for v in vals)
    same_total=sum(v[1] for v in vals)
    if not vals or cross_total<MIN_TEMPORAL_PAIRS_EACH or same_total<MIN_TEMPORAL_PAIRS_EACH:
        return {"temporal_pair_comparison_evaluable": False,
                "n_cross_year_pairs_matched_bins": cross_total,
                "n_same_year_pairs_matched_bins": same_total}
    # Balance same-year and different-year opportunities in each geographical bin.
    weights=np.array([min(x[0],x[1]) for x in vals],float)
    delta=np.array([x[2]-x[3] for x in vals])
    return {
        "temporal_pair_comparison_evaluable": True,
        "n_distance_bins_with_both_year_strata": len(vals),
        "n_cross_year_pairs_matched_bins": cross_total,
        "n_same_year_pairs_matched_bins": same_total,
        "cross_minus_same_year_white_chromatic_discordance": float(np.average(delta,weights=weights)),
        "distance_bin_match_weight": float(weights.sum()),
    }


def summarize_species(g: pd.DataFrame) -> dict:
    n=g.morph.value_counts()
    white=int(n.get("white",0))
    coloured=int(len(g)-white)
    white_qualified=white>=MIN_GLOBAL_MORPH and coloured>=MIN_GLOBAL_MORPH
    hue_qualified=sum(n.get(x,0)>=MIN_GLOBAL_MORPH for x in NONWHITE)>=2
    yy=g.year
    wyears=yy[g.morph=="white"].dropna().unique()
    cyears=yy[g.morph!="white"].dropna().unique()
    global_repeat=(len(wyears)>=2 and len(cyears)>=2)
    row={
        "cohort":g.cohort.iloc[0],
        "inat_taxon_id":int(g.inat_taxon_id.iloc[0]),
        "species":str(g.species.iloc[0]),
        "genus":str(g.species.iloc[0]).split()[0],
        "n_classifiable":int(len(g)),
        "n_dated":int(g.year.notna().sum()),
        "n_observer_known":int((g.observer!="").sum()),
        "n_distinct_years":int(g.year.dropna().nunique()),
        "white_plus_nonwhite_5_each":bool(white_qualified),
        "two_nonwhite_hues_5_each":bool(hue_qualified),
        "both_comparisons_available":bool(white_qualified and hue_qualified),
        "both_white_and_nonwhite_recur_across_years_specieswide":bool(white_qualified and global_repeat),
    }
    for diameter in DIAMETERS_KM:
        metrics=strict_anchor_neighborhoods(g,diameter)
        for key,value in metrics.items():
            row[f"diameter_{int(diameter)}km_{key}"]=value
    row.update(temporal_colour_discordance(g))
    return row


def bootstrap_effect(values: np.ndarray, seed: int) -> list[float] | None:
    if len(values)<5:
        return None
    rng=np.random.default_rng(seed)
    samples=values[rng.integers(0,len(values),(N_BOOT,len(values)))].mean(axis=1)
    return [float(x) for x in np.quantile(samples,(0.025,0.975))]


def describe(df: pd.DataFrame, cohort: str) -> dict:
    both=df.loc[df.both_comparisons_available]
    out={
        "n_species_classifiable":int(len(df)),
        "n_white_plus_nonwhite_5_each":int(df.white_plus_nonwhite_5_each.sum()),
        "n_two_nonwhite_hues_5_each":int(df.two_nonwhite_hues_5_each.sum()),
        "n_both_comparisons_available":int(len(both)),
        "n_white_nonwhite_global_multiyear":int(
            df.both_white_and_nonwhite_recur_across_years_specieswide.sum()),
        "n_species_any_valid_observation_year":int((df.n_dated>0).sum()),
        "n_species_at_least_two_observation_years":int((df.n_distinct_years>=2).sum()),
        "diameters_km":{},
    }
    for diameter in DIAMETERS_KM:
        pre=f"diameter_{int(diameter)}km_"
        wc=df[pre+"white_nonwhite_multi_year_two_observers"].astype(bool)
        hh=df[pre+"two_nonwhite_hues_multi_year_two_observers"].astype(bool)
        a=int((wc & hh & df.both_comparisons_available).sum())
        b=int((wc & ~hh & df.both_comparisons_available).sum())
        c=int((~wc & hh & df.both_comparisons_available).sum())
        neither=int((~wc & ~hh & df.both_comparisons_available).sum())
        p_exact=float(binomtest(min(b,c),b+c,0.5,alternative="two-sided").pvalue) if b+c>0 else None
        out["diameters_km"][str(int(diameter))]={
            "n_white_nonwhite_neighborhood_with_both_labels_at_least_3":int(
                df[pre+"white_nonwhite_same_neighborhood"].sum()),
            "n_white_nonwhite_strict_multiyear_two_observers":int(wc.sum()),
            "n_nonwhite_hue_strict_multiyear_two_observers":int(hh.sum()),
            "n_white_nonwhite_two_years_without_observer_gate":int(
                df[pre+"white_nonwhite_two_years_not_observer_strict"].sum()),
            "joint_opportunity_exact_mcnemar":{
                "n_matched_species":int(len(both)),
                "both_types_recur":a,"only_white_plus_nonwhite_recur":b,
                "only_two_nonwhite_hues_recur":c,"neither_recur":neither,
                "exact_two_sided_p_exploratory":p_exact,
                "interpretation":"paired within-species capability, not equal global colour frequencies"
            }
        }
    t=df.loc[df.temporal_pair_comparison_evaluable]
    vals=t.cross_minus_same_year_white_chromatic_discordance.to_numpy(float)
    rngseed=int.from_bytes(hashlib.sha256((str(SEED)+cohort).encode()).digest()[:8],"little")
    out["distance_matched_cross_year_vs_same_year"]={
        "n_species_evaluable":int(len(vals)),
        "mean_species_delta_discordance":float(vals.mean()) if len(vals) else None,
        "median_species_delta_discordance":float(np.median(vals)) if len(vals) else None,
        "bootstrap_95CI":bootstrap_effect(vals,rngseed),
        "fraction_positive":float((vals>0).mean()) if len(vals) else None,
        "estimand":"cross-year minus same-year white/nonwhite pair discordance, different observers; distance-binned 0-10,10-25,25-50 km",
        "limitation":"year-observer-site confounding and pair dependence remain; not selection or phenotypic stability proof",
    }
    return out


def main() -> None:
    ap=argparse.ArgumentParser()
    for cohort in COHORT_SHA256:
        ap.add_argument("--"+cohort,required=True,type=Path)
    ap.add_argument("--outdir",required=True,type=Path)
    args=ap.parse_args()
    args.outdir.mkdir(parents=True,exist_ok=True)
    output={
        "schema":"fcp_local_multiyear_morph_observation_posthoc_v1",
        "date_jst":"2026-10-08",
        "status":"complete",
        "role":"exploratory_multiobserver_multiyear_photo_evidence_not_fitness",
        "confirmatory_decisions_changed":False,
        "frozen_measured_source_sha256":COHORT_SHA256,
        "time_window_years":[1990,2026],
        "neighborhood_geometry":"anchor-centered circles of half stated diameter; photo-pair maximum <= stated diameter",
        "minimums":{
            "classifiable_per_species":MIN_CLASSIFIABLE,
            "each_morph_species_eligibility":MIN_GLOBAL_MORPH,
            "each_morph_neighborhood":MIN_ANCHOR_MORPH,
            "years_per_morph_in_neighborhood":MIN_DISTINCT_YEARS,
            "independent_observers_per_morph_in_neighborhood":MIN_DISTINCT_OBSERVERS,
            "paired_temporal_discordance_pairs_per_year_stratum":MIN_TEMPORAL_PAIRS_EACH
        },
        "hard_nonclaims":[
            "photo white is not chemically verified absence of pigment and is exposure-coupled",
            "within-diameter photos do not establish members of one naturally mating population",
            "repeated independent observers do not cure site-specific photography or garden biases",
            "photo dates do not establish age or lineage or individual independence",
            "different photographs across years are not recurrence of the same genotype",
            "year-specific photo counts are not an annual census of morph frequencies",
            "no fluctuating selection, adaptive maintenance or fitness from these data",
            "cross-cohort pattern replication remains within the same photo platform",
        ],
        "cohorts":{},
    }
    all_records=[]
    for cohort in COHORT_SHA256:
        df=load(getattr(args,cohort),cohort)
        rows=[summarize_species(group.reset_index(drop=True))
              for _,group in df.groupby("inat_taxon_id",sort=True)]
        metrics=pd.DataFrame(rows)
        output["cohorts"][cohort]=describe(metrics,cohort)
        all_records.append(metrics)
    pd.concat(all_records,ignore_index=True).to_csv(
        args.outdir/"species_multiyear_morph_evidence.csv",index=False)
    (args.outdir/"result.json").write_text(
        json.dumps(output,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(output,ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()
