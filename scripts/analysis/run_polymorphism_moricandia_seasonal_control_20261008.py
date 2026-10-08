#!/usr/bin/env python3
"""Literature-defined positive control: Moricandia arvensis seasonal floral plasticity.

2020/2024 literature: lilac spring flowers, white summer flowers on the same
genotype; floral seasonal polyphenism can be nonadaptive in trait space despite
selection favouring extended summer flowering. This is a photographic diagnosis,
not a phenotype/fitness demonstration or an independent genetic test.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

TARGET_TAXON_ID = 165527
TARGET_SPECIES = "Moricandia arvensis"
VALIDATION_SHA256 = "0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6"
MIN_SEASON_PHOTOS = 5
MIN_BOTH_MORPH_PHOTOS = 3
MIN_SITE_SEASON_PHOTOS = 2
SITE_MAX_DIAMETER_KM = 25.0
EARTH_RADIUS = 6371.0088


def digest(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for part in iter(lambda:f.read(1<<20),b""):
            h.update(part)
    return h.hexdigest()


def true_values(v: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(v):
        return v.fillna(False).astype(bool)
    return v.fillna("").astype(str).str.strip().str.casefold().isin(("true","1","yes","y"))


def classify_season(month: int) -> str:
    if month in (3,4,5):
        return "spring"
    if month in (6,7,8):
        return "summer"
    return "other"


def point_distances(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    lat=np.deg2rad(lat.astype(float))
    lon=np.deg2rad(lon.astype(float))
    c=np.cos(lat)
    vec=np.column_stack([c*np.cos(lon),c*np.sin(lon),np.sin(lat)])
    return np.arccos(np.clip(vec@vec.T,-1,1))*EARTH_RADIUS


def neighborhood_opportunity(d: pd.DataFrame) -> dict:
    """Only describes whether a near-site contains observations in both seasons."""
    if len(d) < 4:
        return {"n_anchors_with_both_seasons":0,"n_anchors_with_expected_colour_season_pair":0}
    geo=point_distances(d.latitude.to_numpy(float),d.longitude.to_numpy(float))
    season=d.season.astype(str).to_numpy()
    white=d.white.to_numpy(bool)
    years=d.year.to_numpy(float)
    obs=d.observer.to_numpy(str)
    n_both=0
    n_pair=0
    n_pair_both_observers=0
    n_pair_across_years=0
    for i in range(len(d)):
        near=geo[i]<=SITE_MAX_DIAMETER_KM/2
        a=near & (season=="spring")
        b=near & (season=="summer")
        if int(a.sum())<MIN_SITE_SEASON_PHOTOS or int(b.sum())<MIN_SITE_SEASON_PHOTOS:
            continue
        n_both+=1
        if not (np.any(a & ~white) and np.any(b & white)):
            continue
        n_pair+=1
        x=obs[a & ~white];y=obs[b & white]
        x=x[x!=""];y=y[y!=""]
        if len(x)>0 and len(y)>0 and len(set(x)|set(y))>=2:
            n_pair_both_observers+=1
        sy=years[a & ~white];su=years[b & white]
        if np.any(np.isfinite(sy)) and np.any(np.isfinite(su)) and (
            len(set(sy[np.isfinite(sy)])|set(su[np.isfinite(su)]))>=2):
            n_pair_across_years+=1
    return {
        "n_anchors_with_both_seasons":n_both,
        "n_anchors_with_expected_colour_season_pair":n_pair,
        "n_expected_pair_anchors_with_at_least_two_observers_total":n_pair_both_observers,
        "n_expected_pair_anchors_spanning_two_observation_years":n_pair_across_years,
        "note":"Overlapping anchors not independent biological populations or hypothesis tests.",
    }


def analyze(data: pd.DataFrame) -> dict:
    required={"inat_taxon_id","species","photo_id","morph","global_classifiable",
              "latitude","longitude","observer_id","observed_on"}
    missing=sorted(required-set(data.columns))
    if missing:
        raise ValueError(f"Missing source fields {missing}")
    sel=true_values(data.global_classifiable) & (
        pd.to_numeric(data.inat_taxon_id,errors="coerce")==TARGET_TAXON_ID)
    d=data.loc[sel,sorted(required)].copy()
    if not d.empty:
        d=d.loc[d.morph.astype(str).isin(("white","yellow_orange","red_pink","blue_purple"))]
    d["latitude"]=pd.to_numeric(d.latitude,errors="coerce")
    d["longitude"]=pd.to_numeric(d.longitude,errors="coerce")
    d["year"]=pd.to_datetime(d.observed_on,errors="coerce").dt.year
    d["month"]=pd.to_datetime(d.observed_on,errors="coerce").dt.month
    d["season"]=d.month.map(lambda x:classify_season(int(x)) if pd.notna(x) else "unknown")
    d["white"]=d.morph.astype(str)=="white"
    d["observer"]=d.observer_id.fillna("").astype(str).str.strip()
    if not d.empty and any(str(s)!=TARGET_SPECIES for s in d.species.unique()):
        raise ValueError("iNaturalist target taxon identity differs from fixed species")
    if d.photo_id.duplicated().any():
        raise ValueError("Duplicate photo IDs in frozen Moricandia subset")
    north=d.loc[d.latitude.between(0,90) & d.longitude.between(-180,180)].copy()
    season_north=north.loc[north.season.isin(("spring","summer"))].copy()
    counts={}
    for season in ("spring","summer","other","unknown"):
        rows=north.loc[north.season==season]
        counts[season]={
            "n_photos":int(len(rows)),
            "n_visible_white":int(rows.white.sum()),
            "n_visible_nonwhite":int((~rows.white).sum()),
            "n_distinct_years":int(rows.year.dropna().nunique()),
            "n_known_observers":int(rows.loc[rows.observer!="","observer"].nunique()),
        }
    spring=counts["spring"]; summer=counts["summer"]
    nwhite=spring["n_visible_white"]+summer["n_visible_white"]
    ncol=spring["n_visible_nonwhite"]+summer["n_visible_nonwhite"]
    support=(spring["n_photos"]>=MIN_SEASON_PHOTOS and
             summer["n_photos"]>=MIN_SEASON_PHOTOS and
             nwhite>=MIN_BOTH_MORPH_PHOTOS and
             ncol>=MIN_BOTH_MORPH_PHOTOS)
    test={"estimable":bool(support)}
    if support:
        table=[[summer["n_visible_white"],summer["n_visible_nonwhite"]],
               [spring["n_visible_white"],spring["n_visible_nonwhite"]]]
        odds,p=fisher_exact(table,alternative="greater")
        odds_corrected=((table[0][0]+0.5)*(table[1][1]+0.5) /
                        ((table[0][1]+0.5)*(table[1][0]+0.5)))
        test.update({
            "summer_white_fraction":float(table[0][0]/sum(table[0])),
            "spring_white_fraction":float(table[1][0]/sum(table[1])),
            "summer_minus_spring_white_fraction":float(
                table[0][0]/sum(table[0])-table[1][0]/sum(table[1])),
            "odds_ratio_raw":float(odds),
            "odds_ratio_haldane_0p5":float(odds_corrected),
            "fisher_one_sided_white_enriched_in_summer_p_exploratory":float(p),
            "counts":table,
        })
    return {
        "schema":"fcp_moricandia_seasonal_positive_control_posthoc_v1",
        "target_taxon_id":TARGET_TAXON_ID,
        "taxon":TARGET_SPECIES,
        "result_status":"positive_control_evaluable" if support else "not_estimable_positive_control",
        "role":"literature_defined_seasonal_plasticity_analogue_not_new_evolutionary_mechanism",
        "source_validation_sha256":VALIDATION_SHA256,
        "n_classifiable_target_photos":int(len(d)),
        "n_northern_hemisphere_georeferenced":int(len(north)),
        "n_southern_or_invalid_coordinate":int(len(d)-len(north)),
        "season_definitions":{
            "spring_months":[3,4,5],"summer_months":[6,7,8],
            "minimum_photos_in_both_seasons":MIN_SEASON_PHOTOS,
            "minimum_visible_white_and_nonwhite_in_seasonal_set":MIN_BOTH_MORPH_PHOTOS
        },
        "seasonal_counts_northern_hemisphere":counts,
        "colour_season_test":test,
        "site_opportunity_northern_hemisphere":neighborhood_opportunity(season_north),
        "literature_control":[
            {"doi":"10.1038/s41467-020-17875-1","claim":"individual spring-lilac to summer-white floral plasticity"},
            {"doi":"10.1093/evlett/qrae017","claim":"floral plasticity maladaptive; indirectly maintained by selection for summer flowering"},
        ],
        "hard_nonclaims":[
            "target is post hoc selected after inspecting FCP local-photo candidate names",
            "Moricandia floral plasticity mechanism is established prior art, not discovered by this test",
            "coarse photographic white may be exposure dependent, including seasonally",
            "photographed plants may not all be wild or from native-range populations",
            "photo-level Fisher independence is limited; p is exploratory not fitness selection inference",
            "a null may arise from missing seasonal photos and does not overturn field experiments",
            "no direct genotypes, same-plant changes, pigment chemistry or pollinator selection measured",
        ],
        "confirmatory_decisions_changed":False,
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--validation",type=Path,required=True)
    ap.add_argument("--outdir",type=Path,required=True)
    args=ap.parse_args()
    if digest(args.validation)!=VALIDATION_SHA256:
        raise RuntimeError("Validation measured-table source hash mismatch")
    r=analyze(pd.read_csv(args.validation,low_memory=False))
    args.outdir.mkdir(parents=True,exist_ok=True)
    (args.outdir/"mor_icandia_seasonal_control.json").write_text(
        json.dumps(r,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps(r,indent=2,ensure_ascii=False))


if __name__=="__main__":
    main()
