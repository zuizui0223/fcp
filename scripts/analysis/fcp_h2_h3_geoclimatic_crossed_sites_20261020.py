#!/usr/bin/env python3
"""Outcome-blind photo-SITE opportunity screen for nested H2/H3 field selection.

Identify geographic/photo-supported localities where the same photographed
species offers both (a) far but two-variable-climate-similar locality pairs and
(b) nearby but two-variable-climate-different locality pairs. No genetic
ancestry, history, true breeding populations, local adaptation or fitness is
inferred. Original coordinates and observation IDs are NEVER exported.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import zipfile

import numpy as np
import pandas as pd
from PIL import Image

PHOTO_SHA = "8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6"
BIO_SHA = "00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf"
EXPECTED_ROWS, EXPECTED_TAXA, EXPECTED_CLASSIFIED = 447500, 4475, 227107
GEOCELL_DEGREES = (1, 2)
SITE_MIN_PHOTOS = (2, 3, 5)
SITE_MIN_OBSERVERS = {2: 2, 3: 2, 5: 3}
BIO_INDICES = (5, 12)
COHORTS = ("discovery", "validation", "prospective", "wave4", "new_U100")
ANCHOR_NAMES = (
    "Campanula punctata", "Hesperis matronalis", "Anemone coronaria",
    "Iris lutescens", "Leptosiphon parviflorus", "Ipomoea purpurea",
    "Cirsium palustre", "Campanula latifolia",
)
USECOLS = (
    "inat_taxon_id", "species", "fcp_source_cohort",
    "photo_id", "observation_id", "observer_id", "latitude", "longitude",
    "positional_accuracy_m", "morph", "global_classifiable",
)
ENV_SCALES = (3.0, 400.0)  # degrees Celsius and mm/year; for site matching only
EARTH_KM = 6371.0088


def check(value, reason):
    if not value:
        raise RuntimeError("FCP_HISTORY_ENV_MATCHING_SOURCE_INVALID: " + reason)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for buf in iter(lambda: f.read(1 << 20), b""):
            h.update(buf)
    return h.hexdigest()


def source_photos(path):
    check(digest(path) == PHOTO_SHA, "unmodified photo measurement SHA mismatch")
    d = pd.read_csv(path, usecols=list(USECOLS), low_memory=False)
    check(len(d) == EXPECTED_ROWS and d.inat_taxon_id.nunique() == EXPECTED_TAXA
          and d.photo_id.nunique() == EXPECTED_ROWS, "full authentic photo denominator mismatch")
    check(d.groupby("inat_taxon_id").fcp_source_cohort.nunique().eq(1).all(),
          "species are not source-cohort disjoint")
    m = d.global_classifiable.fillna(False).astype(str).str.strip().str.lower().isin(("true", "yes", "1"))
    check(int(m.sum()) == EXPECTED_CLASSIFIED, "source classified denominator changed")
    classes = d.loc[m & d.morph.isin(("white", "yellow_orange", "red_pink", "blue_purple"))]
    classified = classes.groupby("inat_taxon_id").morph.value_counts().unstack(fill_value=0)
    classified = classified.reindex(columns=["white", "yellow_orange", "red_pink", "blue_purple"], fill_value=0)
    classified["n_classifiable"] = classified.sum(axis=1)
    classified["n_white"] = classified.white
    classified["n_nonwhite"] = classified.n_classifiable - classified.n_white
    # Photo-label suitability is assessed for site-planning summaries AFTER
    # sites/pairs are selected without using colour labels.
    d["latitude"] = pd.to_numeric(d.latitude, errors="coerce")
    d["longitude"] = pd.to_numeric(d.longitude, errors="coerce")
    d["positional_accuracy_m"] = pd.to_numeric(d.positional_accuracy_m, errors="coerce")
    d = d.loc[d.latitude.between(-90, 90)
              & d.longitude.between(-180, 180)
              & d.positional_accuracy_m.between(0, 5000)
              & d.observation_id.notna() & d.observer_id.notna()].copy()
    d = d.sort_values(["inat_taxon_id", "observation_id", "photo_id"], kind="mergesort")
    d = d.drop_duplicates(["inat_taxon_id", "observation_id"]).reset_index(drop=True)
    return d, classified


def climate_sample(photo, bio_zip):
    check(digest(bio_zip) == BIO_SHA, "frozen WorldClim 2.1 archive SHA mismatch")
    lat = photo.latitude.to_numpy(float)
    lon = photo.longitude.to_numpy(float)
    rr = np.floor((90-lat)*6).astype(np.int64)
    cc = np.floor((lon+180)*6).astype(np.int64)
    valid = (rr>=0)&(rr<1080)&(cc>=0)&(cc<2160)
    with zipfile.ZipFile(bio_zip) as z:
        for bi in BIO_INDICES:
            member = f"wc2.1_10m_bio_{bi}.tif"
            check(member in z.namelist(), f"{member} missing from WorldClim archive")
            with Image.open(io.BytesIO(z.read(member))) as im:
                check(im.size == (2160,1080), "WorldClim raster grid drift")
                meta = im.tag_v2
                check(np.allclose(meta.get(33550,()),(1/6,1/6,0), atol=1e-14,rtol=0),
                      "WorldClim cell scale changed")
                check(np.allclose(meta.get(33922,()),(0,0,0,-180,90,0),atol=1e-14,rtol=0),
                      "WorldClim tiepoint changed")
                check(42113 in meta, "no WorldClim NODATA")
                val = np.array(im)
                no = np.asarray(float(str(meta[42113]).strip("\x00")), dtype=val.dtype).item()
            sample = np.full(len(photo), np.nan)
            raw = val[rr[valid], cc[valid]]
            good = (raw != no)&np.isfinite(raw)
            sample[valid] = np.where(good,raw.astype(float),np.nan)
            if bi == 5:
                check(not np.any(np.isfinite(sample)&((sample < -100)|(sample > 100))),
                      "unmasked impossible BIO5 temperatures")
            if bi == 12:
                check(not np.any(np.isfinite(sample)&((sample < 0)|(sample > 20000))),
                      "unmasked impossible BIO12 annual precipitation")
            photo[f"BIO{bi}"] = sample
    photo = photo.dropna(subset=["BIO5","BIO12"]).copy()
    return photo


def make_sites(photo, deg, minimum):
    p = photo.copy()
    ix = np.clip(np.floor((p.longitude.to_numpy(float)+180)/deg).astype(int),0,360//deg-1)
    iy = np.clip(np.floor((p.latitude.to_numpy(float)+90)/deg).astype(int),0,180//deg-1)
    p["sitecode"] = iy*(360//deg)+ix
    q = (p.groupby(["inat_taxon_id", "sitecode"], sort=True)
           .agg(species=("species","first"),
                cohort=("fcp_source_cohort","first"),
                n_obs=("observation_id","size"),
                n_observer=("observer_id","nunique"),
                lat=("latitude","mean"),
                lon=("longitude","mean"),
                BIO5=("BIO5","mean"),
                BIO12=("BIO12","mean")).reset_index())
    eligible = q[(q.n_obs>=minimum)&(q.n_observer>=SITE_MIN_OBSERVERS[minimum])].copy()
    return eligible


def great_circle_km(lat,lon):
    la = np.deg2rad(np.asarray(lat,float))
    lo = np.deg2rad(np.asarray(lon,float))
    x = np.stack([np.cos(la)*np.cos(lo),np.cos(la)*np.sin(lo),np.sin(la)],axis=1)
    dist = np.clip(x@x.T,-1,1)
    return np.arccos(dist)*EARTH_KM


def crossed_pair_opportunity(g):
    n = len(g)
    if n < 4:
        return {"n_potential_sites":n,"n_pairs":n*(n-1)//2,
                "n_far_similar":0,"n_near_different":0,
                "n_disjoint_geoclimate_crossings":0,"status":"<4_photo_supported_sites"}
    u,v = np.triu_indices(n,k=1)
    dist = great_circle_km(g.lat.to_numpy(float),g.lon.to_numpy(float))[u,v]
    env = g[["BIO5","BIO12"]].to_numpy(float)
    dd = env[u]-env[v]
    ed = np.sqrt(np.sum((dd/np.asarray(ENV_SCALES)[None,:])**2, axis=1))
    check(np.isfinite(dist).all() and np.isfinite(ed).all(),"nonfinite site distances")
    if np.ptp(dist)<1e-9 or np.ptp(ed)<1e-9:
        return {"n_potential_sites":n,"n_pairs":len(u),
                "n_far_similar":0,"n_near_different":0,
                "n_disjoint_geoclimate_crossings":0,
                "status":"distance_or_climate_rank_degenerate"}
    # Site comparisons are defined entirely in geo-climate space.
    # No colour labels are used to define opportunity or pair ranks.
    near_q,far_q=np.quantile(dist,[1/3,2/3])
    sim_q,diff_q=np.quantile(ed,[1/3,2/3])
    far=np.where((dist>=far_q)&(ed<=sim_q))[0]
    near=np.where((dist<=near_q)&(ed>=diff_q))[0]
    n_cross=0
    example=None
    # Only opportunity existence, not all correlated pairs, is the estimand.
    for a in far:
        for b in near:
            if len({int(u[a]),int(v[a]),int(u[b]),int(v[b])})==4:
                n_cross=1
                example=(int(a),int(b))
                break
        if example is not None:break
    strict_far=np.where((dist>=400)&(ed<=1))[0]
    strict_near=np.where((dist<=200)&(ed>=1.5))[0]
    strict_cross=any(len({int(u[a]),int(v[a]),int(u[b]),int(v[b])})==4
                     for a in strict_far for b in strict_near)
    return {
       "n_potential_sites":int(n), "n_pairs":int(len(u)),
       "n_far_similar":int(len(far)),"n_near_different":int(len(near)),
       "n_disjoint_geoclimate_crossings":int(n_cross),
       "strict_n_geographically_far_climate_similar":int(len(strict_far)),
       "strict_n_geographically_near_climate_different":int(len(strict_near)),
       "strict_disjoint_crossing":bool(strict_cross),
       "geo_min_km":float(np.min(dist)),"geo_max_km":float(np.max(dist)),
       "bio5_range_C":float(np.ptp(env[:,0])),
       "bio12_range_mm":float(np.ptp(env[:,1])),
       "status":"crossed_photo_locality_opportunity" if n_cross else "no_crossed_photo_locality_opportunity",
    }


def synthetic_tests():
    a = great_circle_km([0,0],[0,1])
    check(abs(a[0,1]-111.2)<1,"great-circle km arithmetic")
    sites = pd.DataFrame({
      "lat":[0,0,0,0], "lon":[0,1,8,9],
      "BIO5":[20,25,20,25],"BIO12":[1000,1700,1000,1700]
    })
    # Far pairs at similar climate and near pairs at divergent climate.
    o = crossed_pair_opportunity(sites)
    check(o["n_disjoint_geoclimate_crossings"]==1,"crossed-sites synthetic feasibility")
    sites["BIO5"]=20
    sites["BIO12"]=1000
    empty = crossed_pair_opportunity(sites)
    check(empty["n_disjoint_geoclimate_crossings"]==0,"zero climate contrast isn't a crossed design")
    print("FCP_H2_H3_GEOGRAPHY_CLIMATE_CROSSED_SITE_FIXTURES_PASS")


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--photos",type=Path)
    p.add_argument("--worldclim-zip",type=Path)
    p.add_argument("--out",type=Path)
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    synthetic_tests()
    if a.self_test:return
    check(a.photos and a.worldclim_zip and a.out,"input and output paths required")
    photo,cl = source_photos(a.photos)
    before=len(photo)
    photo=climate_sample(photo,a.worldclim_zip)
    original_taxon_names = cl.index
    output={}
    panel=[]
    for deg in GEOCELL_DEGREES:
        for minimum in SITE_MIN_PHOTOS:
            sites=make_sites(photo,deg,minimum)
            byspecies=[]
            for tid,g in sites.groupby("inat_taxon_id",sort=True):
                d=crossed_pair_opportunity(g)
                if len(g)<4:continue
                count=cl.loc[tid] if tid in cl.index else None
                row={
                    "inat_taxon_id":int(tid),"species":str(g.species.iloc[0]),
                    "source_cohort":str(g.cohort.iloc[0]),
                    "classifiable_photo_count":int(count.n_classifiable) if count is not None else 0,
                    "white_photo_count":int(count.n_white) if count is not None else 0,
                    "nonwhite_photo_count":int(count.n_nonwhite) if count is not None else 0,
                    "white_nonwhite_photo_opportunity_ge5":bool(count is not None and count.n_white>=5 and count.n_nonwhite>=5),
                    **d,
                    "species_and_geography_not_native_or_genetically_verified":True,
                }
                byspecies.append(row)
            df=pd.DataFrame(byspecies)
            key=f"{deg}degree_min{minimum}independentobservations"
            if df.empty:
                output[key]={"n_species_ge4_sites":0,"n_crossed_species":0,"n_strict_crossed_species":0}
                continue
            output[key]={
              "n_original_unique_observation_rows_with_bio5_and_bio12":int(len(photo)),
              "n_taxa_with_at_least_one_site":int(sites.inat_taxon_id.nunique()),
              "n_taxa_with_at_least_4_sites":int(len(df)),
              "n_taxa_with_both_crossed_pair_types_disjoint":int(df.n_disjoint_geoclimate_crossings.sum()),
              "n_taxa_with_absolute_threshold_crossed_pairs":int(df.strict_disjoint_crossing.sum()),
              "n_crossed_taxa_with_white_and_nonwhite_photo_labels_ge5":int((df.n_disjoint_geoclimate_crossings.eq(1)&df.white_nonwhite_photo_opportunity_ge5).sum()),
              "n_true_genetic_population_or_historical_connectivity_tests":0,
              "site_density":"at least "+str(minimum)+" independent source observation IDs and at least "+str(SITE_MIN_OBSERVERS[minimum])+" observers in a "+str(deg)+"-degree photo locality",
              "selected_using_colour_outcome":False,
              "primary_comparison_definition":"far distance upper-third and BIO5+BIO12 climate distance lower-third vs near distance lower-third and climate upper-third within species, 4 distinct photo cells",
              "absolute_sensitivity":"far≥400km and climate distance≤1; near≤200km and climate distance≥1.5; climate distance sqrt((dBIO5/3C)^2+(dBIO12/400mm)^2)",
            }
            if deg==1 and minimum==3:
                panel=sorted(byspecies,key=lambda x:(
                    -int(x["n_disjoint_geoclimate_crossings"]),
                    -int(x["strict_disjoint_crossing"]),
                    -int(x["white_nonwhite_photo_opportunity_ge5"]),
                    -int(x["n_potential_sites"]),
                    x["species"],
                ))
    anchor=[]
    selected = {x["species"]:x for x in panel}
    for name in ANCHOR_NAMES:
        x=selected.get(name)
        anchor.append({"species":name,
           "in_1degree_3obs_crossed_site_panel":x is not None,
           "n_potential_sites":int(x["n_potential_sites"]) if x is not None else None,
           "crossed_geo_climate_pairs_4_distinct_sites":int(x["n_disjoint_geoclimate_crossings"]) if x is not None else None,
           "absolute_threshold_crossed_sites":bool(x["strict_disjoint_crossing"]) if x is not None else None,
           "genotype_or_demographic_history_available_in_this_comparison":False,
           "native_wild_status_unverified":True})
    a.out.mkdir(parents=True,exist_ok=True)
    frame=pd.DataFrame(panel)
    if len(frame):
        frame.to_csv(a.out/"global_taxa_crossed_geo_climate_opportunity_NO_COORDINATES.csv",index=False)
    result={
     "schema":"fcp_actual_447500_photo_sites_H2_H3_geoclimatic_crossing_feasibility_v1",
     "original_photo_measurement_sha256":PHOTO_SHA,
     "worldclim_bioclim_archive_sha256":BIO_SHA,
     "original_photo_rows":EXPECTED_ROWS,
     "original_species":EXPECTED_TAXA,
     "original_classifiable_photo_rows":EXPECTED_CLASSIFIED,
     "independent_source_observation_rows_before_bioclim":int(before),
     "independent_source_observation_rows_with_bioclim":int(len(photo)),
     "climate_variables":"WorldClim2.1 long-term BIO5 max warmest-month temperature and BIO12 annual precipitation, 1970-2000",
     "label_outcome_not_used_to_choose_cells_or_pairings":True,
     "spatial_cell_is_a_source_photo_locality_not_a_breeding_population":True,
     "source_era_observer_and_photo_biased":True,
     "no_native_filter_or_historical_genomic_or_adaptation_measurement":True,
     "wild_or_cultivated_and_native_or_introduced_status_unverified":True,
     "source_frame_diagnostics":output,
     "named_anchor_diagnostics":anchor,
     "sensitive_coordinates_and_photo_identifiers_exported":False,
     "claim_ceiling":"sampling/geoclimatic identifiability of field-pair contrasts only; no relative historical vs environmental causal contribution",
    }
    (a.out/"fcp_global_H2_H3_climate_crossing_source_result.json").write_text(
        json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("FCP_REAL_H2_H3_OPPORTUNITY_RESULT",json.dumps(result,sort_keys=True,allow_nan=False))
    print("FCP_H2_H3_FIELD_CONTRAST_CANDIDATE_TOP",json.dumps([
       {"species":x["species"],"n_potential_sites":x["n_potential_sites"],
        "n_disjoint_geoclimate_crossings":x["n_disjoint_geoclimate_crossings"],
        "strict_disjoint_crossing":x["strict_disjoint_crossing"],
        "white_nonwhite_photo_opportunity_ge5":x["white_nonwhite_photo_opportunity_ge5"]}
       for x in panel[:20]],sort_keys=True,allow_nan=False))


if __name__=="__main__":
    main()
