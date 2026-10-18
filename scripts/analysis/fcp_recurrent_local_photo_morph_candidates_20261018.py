#!/usr/bin/env python3
"""Find recurring SAME-SITE, SAME-CALENDAR-MONTH photographic colour contrasts.

No image pixels, per-observation identifiers or exact locality coordinates are
published. A colour contrast is an audit opportunity, NOT evidence of different
individual genotypes, stable breeding-population morphs, balancing selection, or
incipient speciation. This is a descriptive, post hoc screen; no p-values.

For every taxon and each distance ceiling, use all original *classifiable*
genuine photo terminals. A contrast pair requires distinct observation IDs,
distinct nonmissing observers, exactly matched YEAR-MONTH, within radius and
positional accuracy <=20% of radius for EACH photo. A recurring contrast uses
the SAME unordered colour pair in the SAME calendar MONTH in two different
years with FOUR different observers and observation IDs; the complete-link
maximum distance between all FOUR photos is <= the spatial radius (no DBSCAN
chain or coarse-bin artefact). Internal ROI colour-fraction thresholds are
SENSITIVITY screens of photo-class assignment, not validation of flower morphs.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

SOURCE_CSV_SHA256 = "8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6"
SOURCE_PHOTOS = 447500
SOURCE_CLASSIFIED_PHOTOS = 227107
SOURCE_TAXA = 4475
MORPHS = ("white", "yellow_orange", "red_pink", "blue_purple")
RADIUS_KM = (0.25, 0.5, 1.0, 2.0, 5.0)
PURITY_THRESHOLDS = (0.50, 0.55, 0.60, 0.65, 0.70, 0.75)
PALETTE = {
    "white": ("flower_fraction_white",),
    "yellow_orange": ("flower_fraction_yellow", "flower_fraction_orange"),
    "red_pink": ("flower_fraction_red", "flower_fraction_pink", "flower_fraction_magenta"),
    "blue_purple": ("flower_fraction_blue", "flower_fraction_purple"),
}
OTHER = (
    "inat_taxon_id", "species", "photo_id", "observation_id", "observer_id",
    "latitude", "longitude", "positional_accuracy_m", "observed_on",
    "morph", "global_classifiable",
)
COLUMNS = tuple(dict.fromkeys((*OTHER, *(col for cols in PALETTE.values() for col in cols))))


def check(condition, reason):
    if not condition:
        raise ValueError("FCP_RECURRENT_MORPH_OPPORTUNITY_UNVERIFIED: " + reason)


def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as src:
        for block in iter(lambda:src.read(1<<20), b""):
            h.update(block)
    return h.hexdigest()


def read_and_validate(path, *, require_original=True):
    if require_original:
        check(digest(path) == SOURCE_CSV_SHA256,
              "original source-photo SHA256 mismatch")
    d = pd.read_csv(path, usecols=list(COLUMNS), low_memory=False)
    if require_original:
        check(len(d) == SOURCE_PHOTOS and d.inat_taxon_id.nunique() == SOURCE_TAXA
              and d.photo_id.nunique() == SOURCE_PHOTOS,
              "complete photo, taxon or unique photo-ID denominator mismatched")
    d = d.loc[d.global_classifiable.fillna(False).astype(str).str.lower()
              .isin(("true", "1", "yes")) & d.morph.isin(MORPHS)].copy()
    if require_original:
        check(len(d) == SOURCE_CLASSIFIED_PHOTOS,
              "true four-state measured photo support mismatch")
    date = pd.to_datetime(d.observed_on, errors="coerce")
    d = d.loc[date.notna()].copy()
    d["year"] = date[date.notna()].dt.year.astype("int16")
    d["month"] = date[date.notna()].dt.month.astype("int8")
    d["year_month"] = d.year.astype("int32") * 12 + d.month
    for morph, fields in PALETTE.items():
        d["fraction_"+morph] = d[list(fields)].sum(axis=1)
    arr = d[["fraction_"+m for m in MORPHS]].to_numpy(float)
    mask = np.column_stack([d.morph.to_numpy()!=m for m in MORPHS])
    main = np.column_stack([d.morph.to_numpy()==m for m in MORPHS])
    d["primary_fraction"] = np.sum(np.where(main, arr, 0.), axis=1)
    d["max_other_fraction"] = np.max(np.where(mask, arr, -1.), axis=1)
    check(np.isfinite(d.primary_fraction).all() and
          np.isfinite(d.max_other_fraction).all(), "nonfinite ROI colour components")
    return d


def xyz_from_lat_lon(g):
    lat=np.deg2rad(g.latitude.to_numpy(float))
    lon=np.deg2rad(g.longitude.to_numpy(float))
    co=np.cos(lat)
    return np.column_stack((co*np.cos(lon), co*np.sin(lon), np.sin(lat)))


def max_diameter_km(four_xyz):
    u,v=np.triu_indices(4,1)
    return float(np.max(np.arccos(np.clip((four_xyz@four_xyz.T)[u,v],-1.,1.))
                        * 6371.0088))


def strict_colour_pairs(g, xyz, radius):
    u,v=np.triu_indices(len(g),1)
    distance=np.arccos(np.clip((xyz@xyz.T)[u,v],-1.,1.))*6371.0088
    ym=g.year_month.to_numpy()
    observer=g.observer_id.to_numpy()
    obsid=g.observation_id.to_numpy()
    accuracy=g.positional_accuracy_m.to_numpy(float)
    morph=g.morph.to_numpy()
    accmax=min(1000., 200.*radius)
    valid = (
        (distance <= radius) &
        (ym[u] == ym[v]) &
        (observer[u] != observer[v]) & pd.notna(observer[u]) & pd.notna(observer[v]) &
        (obsid[u] != obsid[v]) & pd.notna(obsid[u]) & pd.notna(obsid[v]) &
        np.isfinite(accuracy[u]) & np.isfinite(accuracy[v]) &
        (accuracy[u] <= accmax) & (accuracy[v] <= accmax)
    )
    contrasts=valid & (morph[u]!=morph[v])
    return u[valid],v[valid],u[contrasts],v[contrasts]


def best_recurrence(g, xyz, u, v, radius):
    year=g.year.to_numpy()
    month=g.month.to_numpy()
    observer=g.observer_id.to_numpy()
    obsid=g.observation_id.to_numpy()
    morph=g.morph.to_numpy()
    purity=g.primary_fraction.to_numpy(float)
    secondary=g.max_other_fraction.to_numpy(float)
    by_month_and_pair=defaultdict(list)
    for a,b in zip(u,v):
        key=(int(month[a]),tuple(sorted((str(morph[a]),str(morph[b])))))
        by_month_and_pair[key].append((int(a),int(b)))
    records=[]
    for (calmonth, hues), pairs in by_month_and_pair.items():
        eligible=[]
        for index,(a,b) in enumerate(pairs):
            for c,d in pairs[index+1:]:
                if int(year[a]) == int(year[c]):
                    continue
                four=[a,b,c,d]
                if len(set(obsid[four]))<4:
                    continue
                # Four different observers across the two years, not just
                # different observers within each year.
                if len(set(observer[four])) != 4:
                    continue
                diameter=max_diameter_km(xyz[four])
                if diameter>radius+1e-8:
                    continue
                eligible.append({
                    "taxon_id":int(g.inat_taxon_id.iloc[0]),
                    "species":str(g.species.iloc[0]),
                    "radius_km":float(radius),
                    "calendar_month":int(calmonth),
                    "coarse_colour_pair":" / ".join(hues),
                    "example_year_1":min(int(year[a]),int(year[c])),
                    "example_year_2":max(int(year[a]),int(year[c])),
                    "distinct_observers_in_example":4,
                    "complete_link_diameter_km":diameter,
                    "min_photo_primary_fraction":float(np.min(purity[four])),
                    "max_photo_other_colour_fraction":float(np.max(secondary[four])),
                })
        if not eligible:
            continue
        # Different photo pairs from the same vicinity are NOT independent
        # replicates. Choose one best evidence set per taxon/radius/month/hue.
        eligible.sort(key=lambda x:(
            -x["min_photo_primary_fraction"], x["complete_link_diameter_km"],
            x["example_year_1"], x["example_year_2"]))
        row=eligible[0]
        row["candidate_evidence_set_count_not_independent"]=len(eligible)
        row["possible_photo_hue_mixing"]=True
        for threshold in PURITY_THRESHOLDS:
            row[f"some_four_photo_event_purity_ge_{int(100*threshold)}"] = any(
                e["min_photo_primary_fraction"] >= threshold for e in eligible
            )
        row["some_four_photo_event_alt_le20"] = any(
            e["max_photo_other_colour_fraction"] <= .20 for e in eligible
        )
        records.append(row)
    return records


def analyse(photo):
    tally={float(r):{
        "taxa_with_month_matched_independent_observer_pair":set(),
        "taxa_with_different_colours_nearby":set(),
        "taxa_with_ge5_local_pairs":set(),
        "taxa_with_ge30_local_pairs":set(),
        "independent_observer_pair_instances_not_independent":0,
        "different_colour_pair_instances_not_independent":0
    } for r in RADIUS_KM}
    sites=[]
    for taxon,g in photo.groupby("inat_taxon_id",sort=False):
        g=g.reset_index(drop=True)
        if len(g)<2:
            continue
        xyz=xyz_from_lat_lon(g)
        for radius in RADIUS_KM:
            u,v,a,b = strict_colour_pairs(g,xyz,radius)
            n=len(u)
            if not n:
                continue
            t=tally[radius]
            t["taxa_with_month_matched_independent_observer_pair"].add(int(taxon))
            t["independent_observer_pair_instances_not_independent"]+=n
            t["different_colour_pair_instances_not_independent"]+=len(a)
            if len(a):
                t["taxa_with_different_colours_nearby"].add(int(taxon))
                sites.extend(best_recurrence(g,xyz,a,b,radius))
            if n>=5:
                t["taxa_with_ge5_local_pairs"].add(int(taxon))
            if n>=30:
                t["taxa_with_ge30_local_pairs"].add(int(taxon))
    table=pd.DataFrame(sites)
    check(len(table)>0,"no independent multi-year recurring colour contrasts")
    forbidden={"latitude","longitude","photo_id","observation_id","observer_id","observer",
               "locality","geohash","site_coordinate","site_lat","site_lon"}
    check(not forbidden.intersection(set(table.columns)),
          "sensitive precise locality or source observation ID leaked")
    output={}
    for radius,t in tally.items():
        s=table.loc[table.radius_km==radius]
        output[str(radius)] = {
            "photo_opportunity_species":len(t["taxa_with_month_matched_independent_observer_pair"]),
            "photographic_colour_contrast_species":len(t["taxa_with_different_colours_nearby"]),
            "species_with_ge5_independent_observer_pairs":len(t["taxa_with_ge5_local_pairs"]),
            "species_with_ge30_independent_observer_pairs":len(t["taxa_with_ge30_local_pairs"]),
            "photo_pair_instances_not_independent":t["independent_observer_pair_instances_not_independent"],
            "contrasting_photo_pair_instances_not_independent":t["different_colour_pair_instances_not_independent"],
            "interannual_recurrent_four_observer_species":int(s.taxon_id.nunique()),
            "interannual_species_calendar_month_hue_groups":len(s),
            "interannual_species_purity_threshold_sensitivity":{
                str(int(100*q)):int(
                    s.loc[s[f"some_four_photo_event_purity_ge_{int(100*q)}"],"taxon_id"].nunique()
                ) for q in PURITY_THRESHOLDS
            },
            "interannual_species_with_some_alt_fraction_le20":int(
                s.loc[s.some_four_photo_event_alt_le20,"taxon_id"].nunique()
            ),
        }
    return {
        "schema":"fcp_fivecohort_interannual_local_photo_colour_contrast_candidate_v1",
        "actual_original_measured_photos":SOURCE_PHOTOS,
        "original_four_state_classified_photos":SOURCE_CLASSIFIED_PHOTOS,
        "original_source_taxa":SOURCE_TAXA,
        "genuine_source_sha256":SOURCE_CSV_SHA256,
        "disclaimer":"exploratory candidate opportunities only; not established inherited morphs",
        "independent_observers_across_four_images_required":4,
        "independent_observation_ids_across_four_images_required":4,
        "same_named_hue_pair_and_calendar_month_repeated_in_different_years":True,
        "site_aggregation":"all four original photo locations mutually <= radius; no chaining",
        "position_accuracy_max_fraction_of_radius":0.2,
        "photo_purity":"fraction of segmented flower ROI in the algorithmically assigned coarse colour group, not a morph-certainty score",
        "photo_purity_thresholds_post_hoc":list(PURITY_THRESHOLDS),
        "precise_coordinates_and_source_observation_identifiers_publicly_withheld":True,
        "radius_sensitivity":output,
        "not_evidence_of":{
            "independent_genetic_morphs":True, "same_breeding_population":True,
            "stable_fitness_maintained_FCP":True, "speciation":True,
            "absence_of_bicolour_or_ontogenetic_colour_change":True,
        },
    },table


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--photos",type=Path,required=True)
    parser.add_argument("--output_dir",type=Path,required=True)
    args=parser.parse_args()
    photo=read_and_validate(args.photos)
    result,candidates=analyse(photo)
    args.output_dir.mkdir(parents=True,exist_ok=True)
    (args.output_dir/"fcp_local_recurrence_quality_summary.json").write_text(
        json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    candidates.to_csv(
        args.output_dir/"fcp_repeated_local_photo_contrast_candidates_NO_EXACT_COORDINATES.csv",
        index=False)
    print("FCP_REPEATED_LOCAL_PHOTO_CONTRAST_QA_COMPLETE",
          json.dumps(result["radius_sensitivity"],sort_keys=True))


if __name__=="__main__":
    main()
