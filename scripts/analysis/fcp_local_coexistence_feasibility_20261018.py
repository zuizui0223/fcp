#!/usr/bin/env python3
"""FCP: distinguish range-wide colour spatial sorting from local morph coexistence.

This is a NEW exploratory negative/feasibility check, NOT the original frozen
50-km vertex null. Preserve each species' observed year-month colour-state
counts and all original photo coordinates; shuffle colour labels only WITHIN
the same species and year-month. Require different observers and different
original observation IDs for a spatially near photo pair.

A different photographic hue near another observation is only a *candidate*
for local polymorphism, not evidence of genetically distinct morphs, different
plants, a delimited population, balancing selection or speciation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

PHOTO_SHA256 = "8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6"
COLOUR_STATES = ("white", "yellow_orange", "red_pink", "blue_purple")
MONTH_RADIUS_KM = 5.
MONTH_MAX_POSITION_ACCURACY_M = 1000.
PERMUTATIONS = 199
SOURCE_COLOUR_COUNT = 227107
SOURCE_ROWS = 447500
SOURCE_TAXA = 4475
USECOLS = (
    "inat_taxon_id", "species", "photo_id", "observation_id",
    "observer_id", "latitude", "longitude", "positional_accuracy_m",
    "observed_on", "morph", "global_classifiable",
)


def sha(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as fp:
        for chunk in iter(lambda:fp.read(1<<20), b""):
            h.update(chunk)
    return h.hexdigest()


def require(ok: bool, msg: str) -> None:
    if not ok:
        raise ValueError("FCP_LOCAL_FALSIFICATION_INVALID: " + msg)


def classify(photo: pd.DataFrame) -> pd.DataFrame:
    keep = photo["global_classifiable"].fillna(False).astype(str).str.strip().str.lower().isin(
        ("true", "yes", "1")
    ) & photo["morph"].isin(COLOUR_STATES)
    d=photo.loc[keep].copy()
    dates=pd.to_datetime(d["observed_on"],errors="coerce")
    d["year"]=dates.dt.year
    d["year_month"]=dates.dt.year*12+dates.dt.month
    d["day"]=dates.dt.date
    return d


def great_circle_pairs(g: pd.DataFrame) -> tuple[np.ndarray,np.ndarray,np.ndarray]:
    lat=np.deg2rad(g.latitude.to_numpy(float))
    lon=np.deg2rad(g.longitude.to_numpy(float))
    c=np.cos(lat)
    xyz=np.column_stack((c*np.cos(lon),c*np.sin(lon),np.sin(lat)))
    u,v=np.triu_indices(len(g),k=1)
    distance=np.arccos(np.clip((xyz@xyz.T)[u,v],-1.,1.))*6371.0088
    return u,v,distance


def strict_pairs(g: pd.DataFrame, u: np.ndarray, v: np.ndarray,
                 km: np.ndarray, *, spatial_km: float, date_unit: str,
                 max_accuracy_m: float) -> np.ndarray:
    month=g[date_unit].to_numpy()
    different_observer=g.observer_id.to_numpy()
    original_id=g.observation_id.to_numpy()
    accuracy=g.positional_accuracy_m.to_numpy(float)
    return (
        (km<=spatial_km) &
        (month[u]==month[v]) & pd.notna(month[u]) & pd.notna(month[v]) &
        (different_observer[u]!=different_observer[v]) &
        pd.notna(different_observer[u]) & pd.notna(different_observer[v]) &
        (original_id[u]!=original_id[v]) &
        pd.notna(original_id[u]) & pd.notna(original_id[v]) &
        np.isfinite(accuracy[u]) & np.isfinite(accuracy[v]) &
        (accuracy[u]<=max_accuracy_m) & (accuracy[v]<=max_accuracy_m)
    )


def conditional_month_null(g: pd.DataFrame, u: np.ndarray,
                           v: np.ndarray, *, reps: int=PERMUTATIONS):
    labels=g.morph.to_numpy()
    observed=float(np.mean(labels[u]!=labels[v]))
    month=g.year_month.to_numpy()
    finite=np.unique(month[pd.notna(month)])
    strata=[np.flatnonzero(month==t) for t in finite]
    seed=int.from_bytes(
        hashlib.sha256(f"fcp-local-month-20261018|{int(g.inat_taxon_id.iloc[0])}".encode()).digest()[:8],
        "little"
    )
    rng=np.random.default_rng(seed)
    null=np.empty(reps,float)
    for k in range(reps):
        p=labels.copy()
        for indices in strata:
            if len(indices)>1:
                p[indices]=rng.permutation(labels[indices])
        null[k]=np.mean(p[u]!=p[v])
    return observed,null


def summarise_test(rows: list[dict], nulls: dict[int,np.ndarray], min_pairs: int) -> dict:
    chosen=[r for r in rows if r["strict_pairs_5km_same_month"]>=min_pairs]
    require(len(chosen)>0, "no adequately sampled independent locality-month pairs")
    observed=np.array([r["local_discordance"] for r in chosen],float)
    mat=np.stack([nulls[r["taxon_id"]] for r in chosen])
    require(mat.shape==(len(chosen),PERMUTATIONS), "missing monthly matched null")
    expected=mat.mean(axis=1)
    # Difference is null local pair dissimilarity minus observed local
    # dissimilarity. Selection uses photo graph only, not colour outcome.
    excess=float((expected-observed).mean())
    simulated=(expected[:,None]-mat).mean(axis=0)
    p=float((1+np.count_nonzero(simulated>=excess))/(PERMUTATIONS+1))
    return {
        "min_strict_local_pairs":min_pairs,
        "n_evaluable_species":len(chosen),
        "n_species_with_different_colour_nearby_photographs":int(np.count_nonzero(observed>0)),
        "n_species_with_nonconstant_conditioned_null":int(np.count_nonzero(mat.std(axis=1)>1e-15)),
        "mean_observed_local_dissimilarity":float(observed.mean()),
        "mean_month_conditioned_null_local_dissimilarity":float(expected.mean()),
        "mean_excess_local_similarity_after_month_conditioning":excess,
        "upper_p_199_month_conditioned_nulls":p,
        "null_q025":float(np.quantile(simulated,.025)),
        "null_q975":float(np.quantile(simulated,.975)),
        "is_confirmatory":False,
    }


def audit(photo: pd.DataFrame) -> tuple[dict,pd.DataFrame]:
    require(len(photo)==SOURCE_ROWS and photo.inat_taxon_id.nunique()==SOURCE_TAXA,
            "original photo/source species denominator differs from complete data")
    require(photo.photo_id.nunique()==SOURCE_ROWS, "non-independent original photo ID reused")
    d=classify(photo)
    require(len(d)==SOURCE_COLOUR_COUNT and d.inat_taxon_id.nunique()==SOURCE_TAXA,
            "source real-image phenotype terminal or classifiability changed")
    counts={"5km_same_month":{}, "5km_same_day":{}}
    opp={key:{"pairs":0,"different_colour_pairs":0,
              "n_species_with_local_pairs":0,
              "n_species_with_any_different_colour_pair":0,
              "n_species_with_ge5_local_pairs":0,
              "n_species_with_ge30_local_pairs":0} for key in counts}
    rows=[];nulls={}
    for taxon,g in d.groupby("inat_taxon_id",sort=False):
        g=g.reset_index(drop=True)
        if len(g)<2: continue
        u,v,km=great_circle_pairs(g)
        morph=g.morph.to_numpy()
        for name,unit in (("5km_same_month","year_month"),("5km_same_day","day")):
            yes=strict_pairs(g,u,v,km,spatial_km=MONTH_RADIUS_KM,
                             date_unit=unit,max_accuracy_m=MONTH_MAX_POSITION_ACCURACY_M)
            n=int(yes.sum())
            diff=int(np.count_nonzero(yes & (morph[u]!=morph[v])))
            c=opp[name]
            c["pairs"]+=n
            c["different_colour_pairs"]+=diff
            c["n_species_with_local_pairs"]+=int(n>0)
            c["n_species_with_any_different_colour_pair"]+=int(diff>0)
            c["n_species_with_ge5_local_pairs"]+=int(n>=5)
            c["n_species_with_ge30_local_pairs"]+=int(n>=30)
            if name=="5km_same_month" and n>=5:
                ob,world=conditional_month_null(g,u[yes],v[yes])
                rows.append({
                    "taxon_id":int(taxon), "species":str(g.species.iloc[0]),
                    "strict_pairs_5km_same_month":n,
                    "different_colour_local_pairs":diff,
                    "local_discordance":ob,
                    "month_conditioned_null_mean_dissimilarity":float(world.mean()),
                    "local_excess_similarity":float(world.mean()-ob),
                })
                nulls[int(taxon)]=world
    return {
        "schema":"fcp_real_fivecohort_local_coexistence_feasibility_month_matched_null_v1",
        "source_original_photo_rows":SOURCE_ROWS,
        "source_four_state_classifiable_rows":SOURCE_COLOUR_COUNT,
        "source_original_photo_sha256":PHOTO_SHA256,
        "year_month_stratification":"preserve exact within-species monthly colour-state multiset",
        "same_iNaturalist_provider":True,
        "all_tests_post_hoc":True,
        "null_permutations":PERMUTATIONS,
        "strict_observational_locality_rule":{
            "distance_km":MONTH_RADIUS_KM,
            "max_each_reported_positional_accuracy_m":MONTH_MAX_POSITION_ACCURACY_M,
            "independent_photo_observation_ids":True,
            "different_identified_observers":True,
            "month_or_day_matched_as_named":True,
        },
        "opportunity":opp,
        "month_stratified_5km_local_difference":{
            "at_least_five_pairs":summarise_test(rows,nulls,5),
            "at_least_thirty_pairs":summarise_test(rows,nulls,30),
        },
        "biological_limits":[
            "the photographic coexistence opportunities do not establish genetically distinct morphs",
            "photos may be different observers photographing the same plant or cultivated variants",
            "5km is not verified plant-population identity and may span distinct microhabitats",
            "absence of short-distance null significance does not prove morphs are absent",
            "month-preserving null does not control independent genetic structure or full image exposure",
            "no field frequency, reproductive fitness, neutral genomic gene flow, or reproductive isolation evidence",
            "the original 50km mean result is not changed or replaced by this local post hoc diagnostic",
        ],
    },pd.DataFrame(rows)


def main() -> None:
    ap=argparse.ArgumentParser()
    ap.add_argument("--photos",required=True,type=Path)
    ap.add_argument("--out",required=True,type=Path)
    args=ap.parse_args()
    require(args.photos.is_file() and sha(args.photos)==PHOTO_SHA256,
            "only the authenticated immutable original 447500-photo measured table can be used")
    photo=pd.read_csv(args.photos,usecols=list(USECOLS),low_memory=False)
    result,species=audit(photo)
    args.out.mkdir(parents=True,exist_ok=True)
    out=args.out/"fcp_local_coexistence_month_matched_null_result.json"
    out.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")
    species.to_csv(args.out/"fcp_local_coexistence_5km_month_eligible_taxa.csv",index=False)
    print("FCP_LOCAL_COEXISTENCE_MONTH_CONDITIONED_NULL_COMPLETED",
          json.dumps(result["month_stratified_5km_local_difference"],sort_keys=True))


if __name__=="__main__":
    main()
