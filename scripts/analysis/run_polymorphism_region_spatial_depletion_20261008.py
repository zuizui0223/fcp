#!/usr/bin/env python3
"""FCP within-region geographical sorting, across three species-disjoint cohorts.

Primary: within-species 50-km four-state local colour depletion in low/mid/high
absolute-latitude photographic regions, under each species' exact *regional*
calendar-month colour composition. Source-photo sampling and biological
populations are deliberately distinct inferential targets.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

SHA256 = {
    "discovery": "ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4",
    "validation": "0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6",
    "third": "57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186",
}
MORPHS = ("white", "yellow_orange", "red_pink", "blue_purple")
SCHEMES = {
    "primary": (
        ("low_0_30", 0.0, 30.0),
        ("middle_30_60", 30.0, 60.0),
        ("high_60_90", 60.0, 90.000001),
    ),
    "geographic_tropics_sensitivity": (
        ("tropical_0_23p5", 0.0, 23.5),
        ("extratropical_23p5_60", 23.5, 60.0),
        ("high_60_90", 60.0, 90.000001),
    ),
}
PAIR_POLICIES = ("all", "different_observer")
CONDITIONAL_TIME = ("month", "year_month")
SPATIAL_RADIUS_KM = 50.0
MIN_GLOBAL_CLASSIFIABLE = 40
MIN_REGION_DATED_PHOTOS = 20
MIN_NEAR_PAIRS = 30
MIN_FAR_PAIRS = 30
MIN_EXCHANGEABLE_ROWS = 10
MIN_REGION_IDENTIFIABLE_SPECIES = 25
MIN_REGION_GEOGRAPHIC_SPECIES = 30
PERMUTATIONS = 199
BOOTSTRAPS = 999
EARTH_RADIUS = 6371.0088
MASTER_SEED = 2026100817


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def get_seed(*values: object) -> int:
    encoded = "|".join(map(str, (MASTER_SEED,*values))).encode()
    return int.from_bytes(hashlib.sha256(encoded).digest()[:8], "little")


def boolean(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.casefold().str.strip().isin(
        {"true", "1", "yes", "y"}
    )


def load(path: Path, cohort: str) -> tuple[pd.DataFrame, dict]:
    if digest(path) != SHA256[cohort]:
        raise RuntimeError(f"{cohort}: immutable photo-table checksum mismatch")
    d = pd.read_csv(path,low_memory=False)
    required = {"inat_taxon_id","species","photo_id","latitude","longitude",
                "observed_on","observer_id","morph","global_classifiable"}
    missing = sorted(required-set(d.columns))
    if missing:
        raise RuntimeError(f"{cohort}: missing fields {missing}")
    d = d.loc[boolean(d.global_classifiable) &
              d.morph.astype(str).isin(MORPHS),sorted(required)].copy()
    d["latitude"] = pd.to_numeric(d.latitude,errors="coerce")
    d["longitude"] = pd.to_numeric(d.longitude,errors="coerce")
    d = d.loc[d.latitude.between(-90,90)&d.longitude.between(-180,180)].copy()
    eligible = d.groupby("inat_taxon_id").size()
    d = d.loc[d.inat_taxon_id.isin(
        eligible[eligible>=MIN_GLOBAL_CLASSIFIABLE].index)].copy()
    if d.photo_id.duplicated().any():
        raise RuntimeError(f"{cohort}: photo IDs repeated after filtering")
    date = pd.to_datetime(d.observed_on,errors="coerce")
    valid = date.dt.year.between(1990,2026).fillna(False)
    d["year"] = date.dt.year.where(valid).astype("Int64")
    d["month"] = date.dt.month.where(valid).astype("Int64")
    d["observer"] = d.observer_id.fillna("").astype(str).str.strip()
    d.loc[d.observer.isin(("None","nan","<NA>")), "observer"]=""
    d["abs_latitude"] = d.latitude.abs()
    d["cohort"] = cohort
    d = d.sort_values(["inat_taxon_id","photo_id"],kind="stable").reset_index(drop=True)
    return d,{
        "n_species_classifiable_40":int(d.inat_taxon_id.nunique()),
        "n_classifiable_photo_rows":int(len(d)),
        "n_dated_photo_rows":int(d.month.notna().sum()),
        "n_species_at_least_40_dated":int(
            (d.groupby("inat_taxon_id").month.count()>=40).sum()),
        "n_northern_hemisphere_photos":int((d.latitude>=0).sum()),
        "n_southern_hemisphere_photos":int((d.latitude<0).sum()),
    }


def earth_distance_matrix(lat: np.ndarray,lon: np.ndarray) -> np.ndarray:
    la = np.deg2rad(np.asarray(lat,float))
    lo = np.deg2rad(np.asarray(lon,float))
    c = np.cos(la)
    x = np.stack((c*np.cos(lo),c*np.sin(lo),np.sin(la)),axis=1)
    return np.arccos(np.clip(x@x.T,-1,1))*EARTH_RADIUS


def bands_for(lat: np.ndarray,scheme: str) -> np.ndarray:
    x = np.abs(np.asarray(lat,float))
    if np.any(~np.isfinite(x)) or np.any(x>90):
        raise ValueError("Latitudes must be valid coordinates")
    bands = np.full(len(x), -1,dtype=np.int8)
    for i,(name,lo,hi) in enumerate(SCHEMES[scheme]):
        mask=(x>=lo)&(x<hi)
        if np.any(mask&(bands>=0)):
            raise ValueError(f"{scheme}: overlapping region definitions")
        bands[mask] = i
    if np.any(bands<0):
        raise ValueError("some source latitudes are outside the defined bins")
    return bands


def time_strata(g: pd.DataFrame,mode: str) -> list[np.ndarray]:
    if mode=="month":
        s=g.month.astype(int)
    elif mode=="year_month":
        s=g.year.astype(int)*100+g.month.astype(int)
    else:
        raise ValueError("unknown calendar mode")
    return [np.asarray(v,dtype=int) for v in s.groupby(s,sort=True).indices.values()]


def count_exchangeable(labels: np.ndarray,groups: list[np.ndarray]) -> int:
    return int(sum(len(indices) for indices in groups
                   if len(indices)>=2 and len(np.unique(labels[indices]))>=2))


def shuffle_in_groups(labels: np.ndarray,groups: list[np.ndarray],
                      rng: np.random.Generator) -> np.ndarray:
    p=labels.copy()
    for idx in groups:
        if len(idx)>=2:
            p[idx]=rng.permutation(labels[idx])
    return p


def region_species_test(g: pd.DataFrame,cohort: str,scheme: str,region: str,
                        mode: str,policy: str) -> tuple[dict,np.ndarray] | None:
    g=g.loc[g.month.notna()].copy().reset_index(drop=True)
    if len(g)<MIN_REGION_DATED_PHOTOS:
        return None
    n=len(g)
    labels=pd.Categorical(g.morph,categories=MORPHS).codes.astype(np.int8)
    if np.any(labels<0):
        raise ValueError("Unmapped flower colour")
    u,v=np.triu_indices(n,k=1)
    geo=earth_distance_matrix(g.latitude.to_numpy(float),g.longitude.to_numpy(float))
    within=geo[u,v]<=SPATIAL_RADIUS_KM
    if policy=="different_observer":
        obs=g.observer.to_numpy(str)
        # Overall species-level colour composition stays based on all
        # region photos, as in existing whole-species FCP fixed-count nulls.
        within &= (obs[u]!="")&(obs[v]!="")&(obs[u]!=obs[v])
        far=geo[u,v]>SPATIAL_RADIUS_KM
        far &= (obs[u]!="")&(obs[v]!="")&(obs[u]!=obs[v])
    elif policy=="all":
        far=geo[u,v]>SPATIAL_RADIUS_KM
    else:
        raise ValueError("unrecognised pair policy")
    if int(within.sum())<MIN_NEAR_PAIRS or int(far.sum())<MIN_FAR_PAIRS:
        return None
    iu,iv=u[within],v[within]
    D_global=float(np.mean(labels[u]!=labels[v]))
    D_local=float(np.mean(labels[iu]!=labels[iv]))
    obs_delta=D_global-D_local
    groups=time_strata(g,mode)
    n_swappable=count_exchangeable(labels,groups)
    rng=np.random.default_rng(get_seed(cohort,scheme,region,mode,policy,
                                       int(g.inat_taxon_id.iloc[0])))
    null=np.empty(PERMUTATIONS,float)
    if n_swappable<MIN_EXCHANGEABLE_ROWS:
        null[:]=obs_delta
        identifiable=False
    else:
        for i in range(PERMUTATIONS):
            p=shuffle_in_groups(labels,groups,rng)
            null[i]=D_global-float(np.mean(p[iu]!=p[iv]))
        identifiable=bool(np.ptp(null)>=1e-12)
        if not identifiable:
            null[:]=obs_delta
    residual=0.0 if not identifiable else float(obs_delta-null.mean())
    return {
        "cohort":cohort,"inat_taxon_id":int(g.inat_taxon_id.iloc[0]),
        "species":str(g.species.iloc[0]),
        "genus":str(g.species.iloc[0]).split()[0],
        "scheme":scheme,"region":region,"calendar_mode":mode,"pair_policy":policy,
        "n_region_photos":n,"n_local_pairs":int(within.sum()),
        "n_distant_pairs":int(far.sum()),
        "n_years":int(g.year.nunique()),
        "n_distinct_observers":int(g.loc[g.observer!="","observer"].nunique()),
        "n_exchangeable_photos":n_swappable,
        "identifiable_after_calendar":bool(identifiable),
        "n_photo_white":int((labels==0).sum()),
        "n_photo_nonwhite":int((labels!=0).sum()),
        "regional_D_all":D_global,"regional_D_local":D_local,
        "regional_observed_local_depletion":obs_delta,
        "regional_calendar_null_mean_depletion":float(null.mean()),
        "regional_residual_depletion":residual,
    },null


def summarize(rows: list[dict],nulls: list[np.ndarray],
              cohort: str,scheme: str,region: str,mode: str,policy: str)->dict:
    if not rows:
        return {
            "n_geographically_evaluable_species":0,
            "n_calendar_identifiable_species":0,
            "n_nonidentifiable_kept":0,
            "status":"HOLD_NO_REGIONAL_GEOGRAPHIC_COVERAGE",
            "coverage_pass":False,
            "estimable":False,
        }
    d=pd.DataFrame(rows)
    nullmat=np.vstack(nulls)
    residual=d.regional_residual_depletion.to_numpy(float)
    obs=float(d.regional_observed_local_depletion.mean())
    means=nullmat.mean(axis=0)
    test_null=float(means.mean())
    n=int(len(d))
    nident=int(d.identifiable_after_calendar.sum())
    rng=np.random.default_rng(get_seed("bootstrap",cohort,scheme,region,mode,policy))
    boot=residual[rng.integers(0,n,size=(BOOTSTRAPS,n))].mean(axis=1)
    grouped=d.groupby("genus").regional_residual_depletion.mean()
    ci=[float(v) for v in np.quantile(boot,[.025,.975])]
    p=(1+int(np.count_nonzero(means>=obs)))/(1+PERMUTATIONS)
    pass_gate=(n>=MIN_REGION_GEOGRAPHIC_SPECIES and
               nident>=MIN_REGION_IDENTIFIABLE_SPECIES)
    return {
        "estimable":True,"n_geographically_evaluable_species":n,
        "n_calendar_identifiable_species":nident,
        "n_nonidentifiable_kept":int(n-nident),
        "n_genera":int(len(grouped)),
        "n_local_pairs_across_species":int(d.n_local_pairs.sum()),
        "n_distant_pairs_across_species":int(d.n_distant_pairs.sum()),
        "mean_observed_local_depletion":obs,
        "mean_calendar_null_depletion":test_null,
        "mean_residual_depletion":float(residual.mean()),
        "mean_identifiable_subset_residual":(float(d.loc[
            d.identifiable_after_calendar,"regional_residual_depletion"].mean())
             if nident else None),
        "fraction_species_residual_positive":float((residual>0).mean()),
        "species_bootstrap_95CI":ci,
        "genus_balanced_mean":float(grouped.mean()),
        "permutation_upper_p_unadjusted":float(p),
        "coverage_pass":bool(pass_gate),
        "status":"REGIONAL_ESTIMABLE_EXPLORATORY" if pass_gate
                 else "HOLD_LIMITED_REGION_OR_CONDITIONAL_SUPPORT",
        "nonidentifiable_convention":"no within-stratum reassignment => exact zero identified additional geographic effect, not biological zero",
    }


def holm(values: list[float]) -> list[float]:
    n=len(values)
    order=np.argsort(values)
    result=np.empty(n,float)
    current=0.0
    for rank,index in enumerate(order):
        current=max(current,min(1.0,values[index]*(n-rank)))
        result[index]=current
    return [float(x) for x in result]


def analyze_cohort(d: pd.DataFrame,cohort: str) -> tuple[dict,pd.DataFrame]:
    out={"schemes":{}}
    saved=[]
    for scheme,parts in SCHEMES.items():
        # Add coverage counts independent of whether a local colour effect
        # happens to be estimable, avoiding false biological absences.
        assign=bands_for(d.latitude.to_numpy(float),scheme)
        dscheme=d.copy()
        dscheme["_band"]=assign
        groups=dscheme.groupby("inat_taxon_id",sort=True)
        out["schemes"][scheme]={"regions":{}}
        for bname,lo,hi in parts:
            by_region=dscheme.loc[dscheme.abs_latitude.between(
                lo, hi, inclusive="left")].copy()
            covered=by_region.groupby("inat_taxon_id").agg(
                n=("photo_id","size"),n_dated=("month","count"))
            nany=int(len(covered))
            n20=int((covered.n_dated>=MIN_REGION_DATED_PHOTOS).sum()) if nany else 0
            opportunity={
                "n_species_with_any_region_photos":nany,
                "n_species_with_at_least_20_dated_region_photos":n20,
                "n_region_photos":int(len(by_region)),
                "n_dated_region_photos":int(by_region.month.notna().sum()),
                "minimum_region_photos_per_species":MIN_REGION_DATED_PHOTOS,
            }
            scenarios={}
            for mode in CONDITIONAL_TIME:
                for policy in PAIR_POLICIES:
                    key=f"{mode}__{policy}"
                    metrics=[];nulls=[]
                    for taxon,g in by_region.groupby("inat_taxon_id",sort=True):
                        res=region_species_test(g,cohort,scheme,bname,mode,policy)
                        if res is None:continue
                        row,nul=res
                        metrics.append(row);nulls.append(nul)
                    scenarios[key]=summarize(
                        metrics,nulls,cohort,scheme,bname,mode,policy)
                    saved.extend(metrics)
            out["schemes"][scheme]["regions"][bname]={
                "opportunity":opportunity,"scenarios":scenarios,
            }
        # Within scheme, Holm-correct all 3 regions x 2 observer
        # policies for primary calendar MONTH tests. No favourable
        # region-specific effect is selected after opening outcomes.
        ids=[];raw=[]
        for bname,_,_ in parts:
            for policy in PAIR_POLICIES:
                v=out["schemes"][scheme]["regions"][bname]["scenarios"][f"month__{policy}"]
                ids.append((bname,policy))
                raw.append(v.get("permutation_upper_p_unadjusted",1.0)
                           if v["coverage_pass"] else 1.0)
        adjusted=holm(raw)
        for (bname,policy),pv in zip(ids,adjusted):
            out["schemes"][scheme]["regions"][bname]["scenarios"][
                f"month__{policy}"]["holm_p_across_3regions_x_2observers"]=pv
    return out,pd.DataFrame(saved)


def cross_cohort_decision(output: dict) -> dict:
    decision={}
    for scheme,regions in SCHEMES.items():
        decision[scheme]={}
        for bname,lo,hi in regions:
            per=[]
            for cohort in SHA256:
                elem=output["cohorts"][cohort]["schemes"][scheme]["regions"][bname]
                per.append(elem)
            enough=all(
                all(v["scenarios"][f"month__{p}"]["coverage_pass"]
                    for p in PAIR_POLICIES)
                for v in per
            )
            supported=bool(enough and all(
                all(v["scenarios"][f"month__{p}"].get("holm_p_across_3regions_x_2observers",1)<.05 and
                    v["scenarios"][f"month__{p}"].get("mean_residual_depletion",0)>0 and
                    v["scenarios"][f"month__{p}"].get("species_bootstrap_95CI",[-1,1])[0]>0
                    for p in PAIR_POLICIES)
                for v in per
            ))
            decision[scheme][bname]={
                "cohort_geographic_photo_species":[v["scenarios"]["month__all"][
                    "n_geographically_evaluable_species"] for v in per],
                "cohort_calendar_identifiable_species":[v["scenarios"]["month__all"][
                    "n_calendar_identifiable_species"] for v in per],
                "cohort_month_conditioned_effect":[v["scenarios"]["month__all"].get(
                    "mean_residual_depletion") for v in per],
                "all_cohorts_minimum_30geo_25identified":bool(enough),
                "month_conditioned_positive_replication_supported":supported,
                "status":"SUPPORT_SAME_SOURCE_REGIONAL_GEOGRAPHIC_STRUCTURE" if supported
                         else "HOLD_REGION_COVERAGE" if not enough
                         else "NOT_SUPPORTED_AS_REPEATED_REGIONAL_STRUCTURE",
            }
    return decision


def main() -> None:
    ap=argparse.ArgumentParser()
    for cohort in SHA256:
        ap.add_argument("--"+cohort,type=Path,required=True)
    ap.add_argument("--outdir",required=True,type=Path)
    args=ap.parse_args()
    args.outdir.mkdir(parents=True,exist_ok=True)
    results={
        "schema":"fcp_latitudinal_region_spatial_depletion_posthoc_v1",
        "date_jst":"2026-10-08",
        "status":"complete_source_verified",
        "role":"retrospective_cross_species_region_conditional_geography",
        "confirmatory_decisions_changed":False,
        "source_sha256":SHA256,
        "same_provider_same_photo_colour_measurement":True,
        "species_disjoint_across_three_cohorts":True,
        "radius_km":SPATIAL_RADIUS_KM,
        "primary_bands":"low 0-30 degrees; middle 30-60; high 60-90 absolute latitude",
        "sensitivity_bands":"tropical 0-23.5; extratropical 23.5-60; high 60-90",
        "regional_species_classifiable_min_global":MIN_GLOBAL_CLASSIFIABLE,
        "minimum_dated_photos_within_region":MIN_REGION_DATED_PHOTOS,
        "minimum_pair_opportunity_each":"at least 30 <=50km and at least 30 >50km conspecific photographed pairs WITHIN the SAME region",
        "minimum_identifiable_species_per_band":MIN_REGION_IDENTIFIABLE_SPECIES,
        "minimum_geo_species_per_band":MIN_REGION_GEOGRAPHIC_SPECIES,
        "permutations":PERMUTATIONS,"species_bootstraps":BOOTSTRAPS,
        "hard_nonclaims":[
            "sampled photos and their sites are not an unbiased worldwide biological population census",
            "low-latitude 0-30 degrees includes subtropical photo sites and is not purely tropical",
            "photo-space <=50km regions are not genetically demonstrated mating populations",
            "regional comparison is within exactly the same species' regional photos, not between species in different climates",
            "photo dates and month-controlled labels do not identify flowering individual age or season in opposite hemisphere",
            "different species may occur in more than one region; do not treat between-region comparisons as independent tests",
            "a month-nonexchangeable regional null contributes zero IDENTIFIABLE extra structure, not biologically absent structure",
            "northern and southern hemisphere sampling, citizen-science effort and temperate species representation can differ",
            "white photo classification has known exposure association, with no direct tissue pigment or reflectance control",
            "neither selection, genetics, UV protection, pollinator role, nor lineage-specific evolutionary mechanisms are estimated",
            "only same-image-provider three-cohort result; does not establish all-angiosperm world generality",
        ],
        "source_coverage":{},
        "cohorts":{},
        "cross_cohort":{},
    }
    details=[]
    for cohort in SHA256:
        d,coverage=load(getattr(args,cohort),cohort)
        results["source_coverage"][cohort]=coverage
        payload,rows=analyze_cohort(d,cohort)
        results["cohorts"][cohort]=payload
        details.append(rows)
    results["cross_cohort"]=cross_cohort_decision(results)
    if any(len(t) for t in details):
        pd.concat(details,ignore_index=True).to_csv(
            args.outdir/"regional_species_conditional_effects.csv",index=False)
    else:
        pd.DataFrame().to_csv(args.outdir/"regional_species_conditional_effects.csv",index=False)
    (args.outdir/"result.json").write_text(
        json.dumps(results,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(results,ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()
