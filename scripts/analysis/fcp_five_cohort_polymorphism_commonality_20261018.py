#!/usr/bin/env python3
"""Within-species floral-colour polymorphism commonality across five REAL FCP cohorts.

Interpretation is photograph-conditional, not population morph frequencies.
All geographic nulls use the original per-species 199 label permutations,
historical fixed seed, classifiable/50-km local-pair gates, and observed site graph.
Post hoc morphology-blind majority/other, D-strata, genus, observer/year and
white-exclusion comparisons cannot become prospective mechanistic claims.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.run_polymorphism_distributed_polymorphism_posthoc_20261007 import (
    MORPHS, MORPH_CODE, PERMUTATIONS, PRIMARY_RADIUS, MIN_LOCAL_PAIRS,
    as_bool, analyse_species, stable_seed, pairwise_geo_km,
)

PHOTOS_SHA = "8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6"
SPATIAL_SHA = "960e8c4bae8a9ecb01fefbef06b6b88cd8fdc63d86af3e4c354cab4eb1aedd35"
SPATIAL_RESULT_SHA = "79aeab0555ec2efd865606a9bdc3e06574a5168c38029ee16e22b945e1a9eacf"
COHORTS = ("discovery", "validation", "prospective", "wave4", "new_U100")
SEED_COHORT = {"discovery": "discovery", "validation": "validation",
               "prospective": "third", "wave4": "third", "new_U100": "third"}
PHOTOS_PER_SPECIES = 100


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as fd:
        for chunk in iter(lambda: fd.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def demand(flag, message):
    if not flag:
        raise ValueError("FCP_COMMONALITY_SOURCE_OR_NULL_INVALID: " + message)


def mean_null(observed, null):
    obs = np.asarray(observed, float)
    z = np.asarray(null, float)
    demand(obs.ndim == 1 and len(obs) > 1 and z.shape == (len(obs), PERMUTATIONS),
           "invalid independent species-level null matrix")
    demand(np.isfinite(obs).all() and np.isfinite(z).all(), "nonfinite observed or null values")
    world = z.mean(axis=0)
    value = float(obs.mean())
    return {
        "n_species": len(obs),
        "mean_depletion": value,
        "median_depletion": float(np.median(obs)),
        "positive_species_fraction": float(np.mean(obs > 0)),
        "null_mean": float(world.mean()),
        "null_q025": float(np.quantile(world, 0.025)),
        "null_q975": float(np.quantile(world, 0.975)),
        "p_upper": float((1 + np.count_nonzero(world >= value)) / (PERMUTATIONS + 1)),
    }


def all_species_counts(photo):
    valid = photo.loc[as_bool(photo.global_classifiable) & photo.morph.isin(MORPHS)]
    freq = valid.groupby(["inat_taxon_id", "species", "fcp_source_cohort", "morph"]).size().unstack(fill_value=0)
    for morph in MORPHS:
        if morph not in freq:
            freq[morph] = 0
    counts = freq[MORPHS]
    n = counts.sum(axis=1)
    shares = counts.div(n, axis=0)
    data = pd.DataFrame({
        "n_classifiable": n, "n_observed_states": (counts > 0).sum(axis=1),
        "D": 1 - (shares ** 2).sum(axis=1),
        "majority_fraction": shares.max(axis=1),
        "majority_morph": shares.idxmax(axis=1),
        "genus": [index[1].split()[0] for index in counts.index],
    }, index=counts.index).reset_index()
    return data


def species_profile(data):
    d = data.D.to_numpy(float)
    m = data.majority_fraction.to_numpy(float)
    return {
        "n_species": len(data), "n_genera": data.genus.nunique(),
        "mean_D": float(d.mean()), "median_D": float(np.median(d)),
        "D_q25": float(np.quantile(d, .25)), "D_q75": float(np.quantile(d, .75)),
        "median_photo_majority_fraction": float(np.median(m)),
        "photo_majority_ge75pct": float(np.mean(m >= .75)),
        "photo_majority_ge90pct": float(np.mean(m >= .9)),
        "photo_majority_ge95pct": float(np.mean(m >= .95)),
        "n_observed_states": {str(k): int((data.n_observed_states == k).sum()) for k in (1,2,3,4)},
        "majority_state_counts": {c: int((data.majority_morph == c).sum()) for c in MORPHS},
    }


def test_group(rows, values, nulls):
    indices = rows.index.to_numpy(dtype=int)
    return mean_null(values[indices], nulls[indices])


def fixed_label_perm(g, cohort):
    # The old code sorts within species on photo ID before RNG permutation.
    g = g.sort_values("photo_id", kind="stable").reset_index(drop=True)
    labels = np.array([MORPH_CODE[str(x)] for x in g.morph], dtype=np.int8)
    # Only advance ONE frozen RNG stream, identical to analyse_species().
    rng = np.random.default_rng(stable_seed(cohort, int(g.inat_taxon_id.iloc[0])))
    null = np.stack([rng.permutation(labels) for _ in range(PERMUTATIONS)])
    return g, labels, null


def obs_and_null_over_edges(labels, permutations, u, v, pair_diversity, edges):
    if int(edges.sum()) < MIN_LOCAL_PAIRS:
        return None
    uu, vv = u[edges], v[edges]
    observed = pair_diversity - float(np.mean(labels[uu] != labels[vv]))
    null = pair_diversity - (permutations[:, uu] != permutations[:, vv]).mean(axis=1)
    return observed, null


def calculate(photo, spatial, original):
    eligible_40 = all_species_counts(photo).loc[lambda d: d.n_classifiable >= 40].copy()
    demand(len(eligible_40) == 3297, "original four-state species eligibility drift")
    demand(int(as_bool(photo.global_classifiable).sum()) == 227107, "real classifiable photo total drift")
    demand(photo.inat_taxon_id.nunique() == 4475 and len(photo) == 447500, "full original pixels missing")
    combined = eligible_40.merge(
        spatial[["inat_taxon_id", "D_pair", "depletion_50km", "local_pairs_50km", "n_classifiable"]],
        on="inat_taxon_id", validate="one_to_one", suffixes=("", "_old"))
    demand(np.all(combined.n_classifiable == combined.n_classifiable_old), "species classified rows changed")
    require_d = combined.D_pair * (combined.n_classifiable - 1) / combined.n_classifiable
    demand(np.allclose(combined.D, require_d, atol=1e-11), "original finite-sample D not reproduced")
    valid_spatial = spatial.loc[spatial.depletion_50km.notna()]
    demand(len(valid_spatial) == 1761, "true 50km eligible species changed")
    ids = set(valid_spatial.inat_taxon_id)
    by_taxon = valid_spatial.set_index("inat_taxon_id")
    photo_sub = photo.loc[photo.inat_taxon_id.isin(ids) &
                          as_bool(photo.global_classifiable) & photo.morph.isin(MORPHS)]
    observed, binary, null_original, null_binary = [], [], [], []
    records = []
    filters = {"different_observer": [], "different_year": []}

    for taxon, raw in photo_sub.groupby("inat_taxon_id", sort=True):
        source = str(raw.fcp_source_cohort.iloc[0])
        g, labels, permutations = fixed_label_perm(raw, SEED_COHORT[source])
        old = by_taxon.loc[taxon]
        row, old_worlds = analyse_species(g, SEED_COHORT[source])
        demand(int(row["local_pairs_50km"]) == int(old.local_pairs_50km) and
               abs(row["depletion_50km"] - old.depletion_50km) < 1e-11 and
               abs(row["D_pair"] - old.D_pair) < 1e-11 and
               abs(old_worlds[50.].mean() - old.null_mean_depletion_50km) < 1e-11,
               "original per-species colour/geometry/permutation result diverged")
        u, v = np.triu_indices(len(labels), k=1)
        distance = pairwise_geo_km(g.latitude.to_numpy(float), g.longitude.to_numpy(float))[u,v]
        local = distance <= 50
        counted = np.bincount(labels, minlength=4)
        dominant = int(np.argmax(counted))   # tie break uses fixed original morph order
        majority = (labels == dominant)
        pair_binary = float(np.mean(majority[u] != majority[v]))
        binary_result = obs_and_null_over_edges(
            majority, permutations == dominant, u, v, pair_binary, local)
        demand(binary_result is not None, "binary support cannot be missing")
        observed.append(float(row["depletion_50km"]))
        null_original.append(old_worlds[50.])
        binary.append(binary_result[0])
        null_binary.append(binary_result[1])
        records.append({"taxon_id":int(taxon), "cohort":source, "genus":str(g.species.iloc[0]).split()[0],
                        "n_classifiable":len(g), "four_state_depletion":float(row["depletion_50km"]),
                        "majority_vs_other_depletion":float(binary_result[0]), "n_local_pairs":int(local.sum())})
        obsid = g.observer_id.to_numpy()
        ovalid = ~pd.isna(obsid)
        years = pd.to_datetime(g.observed_on, errors="coerce").dt.year.to_numpy(float)
        yvalid = np.isfinite(years)
        for key, edges in (
            ("different_observer", local & (obsid[u] != obsid[v]) & ovalid[u] & ovalid[v]),
            ("different_year", local & (years[u] != years[v]) & yvalid[u] & yvalid[v]),
        ):
            result = obs_and_null_over_edges(labels, permutations, u, v, float(row["D_pair"]), edges)
            if result is not None:
                filters[key].append((source, result[0], result[1]))

    enriched = pd.DataFrame(records).merge(
        combined[["inat_taxon_id","D","majority_fraction","majority_morph","n_observed_states"]],
        left_on="taxon_id", right_on="inat_taxon_id", validate="one_to_one").drop(columns="inat_taxon_id")
    four = np.asarray(observed,float)
    bb = np.asarray(binary,float)
    z4 = np.stack(null_original)
    zb = np.stack(null_binary)
    pooled = mean_null(four,z4)
    baseline = original["pooled_original_50km"]
    demand(abs(pooled["mean_depletion"]-baseline["mean_depletion"]) < 1e-11
           and pooled["p_upper"]==baseline["mean_depletion_p_upper"], "old full mean permutation not reproducible")
    by_cohort = {}
    for c in COHORTS:
        sub = enriched.loc[enriched.cohort == c]
        by_cohort[c] = {
            "four_state": test_group(sub,four,z4),
            "majority_vs_other": test_group(sub,bb,zb),
        }
        frozen = original["all_sources"][c]["50"]
        demand(abs(by_cohort[c]["four_state"]["mean_depletion"] - frozen["mean_depletion"])<1e-11
               and by_cohort[c]["four_state"]["p_upper"]==frozen["mean_depletion_p_upper"],
               "old cohort 50km statistics changed")
    by_colour = {}
    for morph in MORPHS:
        sub = enriched.loc[enriched.majority_morph == morph]
        by_colour[morph] = {"four_state":test_group(sub,four,z4),
                            "majority_vs_other":test_group(sub,bb,zb)}
    enriched["D_quartile"] = pd.qcut(enriched.D.rank(method="first"),4,labels=["Q1","Q2","Q3","Q4"]).astype(str)
    quartiles = {}
    for group, sub in enriched.groupby("D_quartile"):
        quartiles[group] = {"D_range":[float(sub.D.min()),float(sub.D.max())],
                            "four_state":test_group(sub,four,z4),
                            "majority_vs_other":test_group(sub,bb,zb)}
    by_genus = list(enriched.groupby("genus").indices.values())
    genus = {
        "n_genera":len(by_genus),
        "four_state":mean_null(np.array([four[ix].mean() for ix in by_genus]),
                               np.vstack([z4[ix].mean(axis=0) for ix in by_genus])),
        "majority_vs_other":mean_null(np.array([bb[ix].mean() for ix in by_genus]),
                                     np.vstack([zb[ix].mean(axis=0) for ix in by_genus])),
    }
    restricted = enriched.loc[enriched.n_classifiable >= 60]
    high_support = {"four_state":test_group(restricted,four,z4),
                    "majority_vs_other":test_group(restricted,bb,zb)}
    observed_controls = {}
    for key, rows in filters.items():
        gobs = np.array([r[1] for r in rows],float)
        gnull = np.stack([r[2] for r in rows])
        observed_controls[key] = {"pooled":mean_null(gobs,gnull), "by_cohort":{}}
        for c in COHORTS:
            ii = [i for i,r in enumerate(rows) if r[0]==c]
            observed_controls[key]["by_cohort"][c] = mean_null(gobs[ii],gnull[ii])

    # The nonwhite-only analysis has a *different* eligible population and
    # three coarse states; no claim that it reproduces each old nonwhite null.
    nonwhite_records = []
    nonwhite_photo = photo.loc[as_bool(photo.global_classifiable) &
                               photo.morph.isin(MORPHS) & photo.morph.ne("white")]
    for tid, raw in nonwhite_photo.groupby("inat_taxon_id",sort=True):
        if len(raw) < 40:
            continue
        raw = raw.dropna(subset=["latitude","longitude"])
        if len(raw) < 40:
            continue
        c=str(raw.fcp_source_cohort.iloc[0])
        g = raw.sort_values("photo_id",kind="stable").reset_index(drop=True)
        row, null = analyse_species(g,SEED_COHORT[c])
        if not np.isfinite(row["depletion_50km"]):
            continue
        nonwhite_records.append((c,float(row["depletion_50km"]),null[50.]))
    require = len(nonwhite_records) > 3
    demand(require,"no nonwhite species")
    nobs = np.array([z[1] for z in nonwhite_records])
    nnul = np.stack([z[2] for z in nonwhite_records])
    observed_controls["nonwhite_only"] = {
        "pooled":mean_null(nobs,nnul),
        "by_cohort":{c:mean_null(nobs[[i for i,row in enumerate(nonwhite_records) if row[0]==c]],
                                 nnul[[i for i,row in enumerate(nonwhite_records) if row[0]==c]])
                     for c in COHORTS},
    }
    return {
        "schema":"fcp_real_five_cohort_intraspecific_polymorphism_commonality_v1",
        "primary_inferential_status":"post_hoc_same_iNaturalist_source_and_unchanged_original_null",
        "frozen_original_four_state_50km_reproduced":True,
        "source_photo_sha256":PHOTOS_SHA, "source_spatial_sha256":SPATIAL_SHA,
        "photo_100_rows_species":4475, "original_photo_rows":447500,
        "classifiable_photo_rows":227107, "null_worlds":199,
        "minimum_classifiable_photos":40, "min_50km_local_pairs":30,
        "radius_50km_post_hoc":True,
        "photo_state_profile":{"all_ge40":species_profile(eligible_40),
                               "ge60":species_profile(eligible_40.loc[eligible_40.n_classifiable >= 60]),
                               "by_cohort":{c:species_profile(eligible_40.loc[eligible_40.fcp_source_cohort==c])
                                            for c in COHORTS}},
        "spatial_frozen_four_state":pooled,
        "spatial_majority_vs_other":mean_null(bb,zb),
        "by_cohort":by_cohort, "by_photo_majority_colour":by_colour,
        "by_frozen_species_global_D_quartile":quartiles,
        "genus_equal_weight":genus, "ge60_classifiable":high_support,
        "technical_falsifications":observed_controls,
        "hard_nonclaims":[
            "photographed state frequencies are not wild-population morph frequencies",
            "same photographic source is not an independent data-provider replication",
            "geographical photo graph is not a delimited genetic population census",
            "no experimental measure of drift, fitness, local adaptation or pollinator selection",
            "dominant morph identity, quartile thresholds and this commonality analysis are post hoc",
            "different-observer, different-year and white exclusion do not eliminate all photographic biases",
            "nonwhite-only results need not be supported in each cohort",
            "multiple exploratory subgroup p-values are not universally multiplicity-corrected",
        ],
    }, enriched, eligible_40


def main():
    p=argparse.ArgumentParser()
    for arg in ("photos","spatial","spatial_result","receipt","out"):
        p.add_argument("--"+arg,type=Path,required=True)
    args=p.parse_args()
    demand(sha(args.photos)==PHOTOS_SHA, "actual merged photo bytes SHA mismatch")
    demand(sha(args.spatial)==SPATIAL_SHA and sha(args.spatial_result)==SPATIAL_RESULT_SHA,
           "original spatial table/null JSON SHA mismatch")
    receipt=json.loads(args.receipt.read_text())
    demand(receipt.get("merged_photo_table_sha256")==PHOTOS_SHA and
           receipt.get("total_original_pixel_terminal_photo_rows")==447500, "source receipt invalid")
    use=["inat_taxon_id","species","fcp_source_cohort","photo_id","latitude",
         "longitude","morph","global_classifiable","observer_id","observed_on"]
    photo=pd.read_csv(args.photos,usecols=use,low_memory=False)
    species=pd.read_csv(args.spatial,low_memory=False)
    prior=json.loads(args.spatial_result.read_text())
    result,spatial,profile=calculate(photo,species,prior)
    args.out.mkdir(parents=True,exist_ok=True)
    (args.out/"fcp_five_cohort_polymorphism_commonality.json").write_text(
        json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")
    spatial.to_csv(args.out/"fcp_five_cohort_commonality_spatial_species.csv",index=False)
    profile.to_csv(args.out/"fcp_five_cohort_commonality_photo_state_species.csv",index=False)
    print("FCP_REAL_FIVE_COHORT_COMMONALITY_COMPLETED",
          json.dumps({k:result[k] for k in ("spatial_frozen_four_state",
                   "spatial_majority_vs_other","genus_equal_weight","technical_falsifications")},sort_keys=True))


if __name__=="__main__":
    main()
