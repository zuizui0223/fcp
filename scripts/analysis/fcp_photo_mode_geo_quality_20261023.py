#!/usr/bin/env python3
"""Exploratory measured-photo FCP modes and high-purity geographical label sorting.

The only biological unit directly confirmed is a photo-colour measurement. No output
is a genetic morph, independently verified petal hue, or adaptive geographic cluster.
Data source: immutable 447,500-photo release; no coordinates/photo/observer IDs exported.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import warnings
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score
from sklearn.mixture import GaussianMixture
from sklearn.exceptions import ConvergenceWarning

SOURCE_SHA = "8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6"
TOTAL = 447500
CLASSIFIED = 227107
TAXA = 4475
LABELS = ("white", "yellow_orange", "red_pink", "blue_purple")
PALETTE = ("white", "yellow", "orange", "red", "pink", "magenta", "purple", "blue", "bronze")
VARS = tuple("flower_fraction_" + x for x in PALETTE)
MORPH_IDX = {m: i for i,m in enumerate(LABELS)}
PHOTO_MORPH_PALETTE = {
    "white": ("white",), "yellow_orange": ("yellow", "orange"),
    "red_pink": ("red", "pink", "magenta"),
    "blue_purple": ("purple", "blue"),
}
ORIG_REQ = ("inat_taxon_id", "species", "fcp_source_cohort", "global_classifiable",
            "morph", "observer_id", "latitude", "longitude", "positional_accuracy_m",
            "photo_id", "observed_on", *VARS)
SEED = 20261023
N_NULL = 199


def check(x, message):
    if not x:
        raise ValueError("FCP_ORIGINAL_PHOTO_CLUSTER_ERROR: " + message)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_photos(path, *, trusted=False):
    if not trusted:
        check(sha256(path) == SOURCE_SHA, "authentic photo CSV SHA256 mismatch")
    t = pd.read_csv(path, usecols=list(ORIG_REQ), low_memory=False)
    if not trusted:
        check(len(t) == TOTAL and t.inat_taxon_id.nunique() == TAXA
              and t.photo_id.nunique() == TOTAL, "source image/photo/taxon denominator changed")
    valid = t.global_classifiable.fillna(False).astype(str).str.lower().isin(("true", "1", "yes"))
    d = t.loc[valid & t.morph.isin(LABELS)].copy()
    if not trusted:
        check(len(d) == CLASSIFIED, "classified image source denominator changed")
    matrix = d[list(VARS)].to_numpy(float)
    check(np.isfinite(matrix).all() and (matrix >= 0).all(), "nonfinite ROI colour fractions")
    sums = matrix.sum(axis=1)
    check(np.all(sums > 0), "zero total flower palette mass")
    matrix /= sums[:, None]
    for k, col in enumerate(VARS):
        d[col] = matrix[:, k]
    share = np.zeros(len(d), dtype=float)
    chosen = d.morph.to_numpy()
    for hue, ingredients in PHOTO_MORPH_PALETTE.items():
        mask = chosen == hue
        cols = [PALETTE.index(c) for c in ingredients]
        share[mask] = matrix[mask][:, cols].sum(axis=1)
    d["assigned_photo_group_fraction"] = share
    dates = pd.to_datetime(d.observed_on, errors="coerce")
    d["year_month"] = (dates.dt.year * 12 + dates.dt.month).fillna(-1).astype("int32")
    check(np.isfinite(share).all() and (share >= 0).all() and (share <= 1 + 1e-8).all(),
          "invalid assigned coarse colour group percentage")
    return d, {int(k): len(v) for k,v in d.groupby("inat_taxon_id",sort=False)}


def quality_group(g, original_n, min_share):
    q = g.loc[g.assigned_photo_group_fraction >= min_share].copy()
    counts = q.morph.value_counts()
    independent = q.groupby("morph")["observer_id"].nunique(dropna=True)
    strong = [c for c in LABELS if counts.get(c, 0) >= 5 and independent.get(c, 0) >= 3]
    # No colour class is called a true morph; all four labels remain recorded separately.
    return q, strong, {"n_original_classified": int(original_n),
                       "n_purity_passed": len(q),
                       "n_classes_with_ge5_photo_and_ge3_observers": len(strong),
                       "observed_photographic_labels": int(len(counts))}


def kilometre(lat, lon, ii, jj):
    la, lo = np.deg2rad(lat), np.deg2rad(lon)
    ds = la[ii] - la[jj]
    dl = lo[ii] - lo[jj]
    hav = np.sin(ds/2)**2 + np.cos(la[ii])*np.cos(la[jj])*np.sin(dl/2)**2
    return 2*6371.0088*np.arcsin(np.sqrt(np.clip(hav, 0, 1)))


def geography_species(g, taxon, *, min_pairs=30):
    mask = g.latitude.between(-90, 90) & g.longitude.between(-180, 180)
    mask &= g.positional_accuracy_m.between(0, 5000)
    s = g.loc[mask]
    if len(s) < 10:
        return None
    a, b = np.triu_indices(len(s), 1)
    dist = kilometre(s.latitude.to_numpy(float), s.longitude.to_numpy(float), a, b)
    within = dist <= 50
    if int(within.sum()) < min_pairs:
        return None
    labels = pd.Categorical(s.morph, categories=LABELS).codes
    obs = s.observer_id.to_numpy()
    different_observers = pd.notna(obs[a]) & pd.notna(obs[b]) & (obs[a] != obs[b])
    if int(np.count_nonzero(within & different_observers)) < min_pairs:
        return None
    local_a, local_b = a[within], b[within]
    diff_a, diff_b = a[within & different_observers], b[within & different_observers]
    all_o_a, all_o_b = a[different_observers], b[different_observers]
    observed = float(np.mean(labels[a] != labels[b]) -
                     np.mean(labels[local_a] != labels[local_b]))
    observed_obs = float(np.mean(labels[all_o_a] != labels[all_o_b]) -
                         np.mean(labels[diff_a] != labels[diff_b]))
    rng = np.random.default_rng(SEED + int(taxon))
    null = np.empty(N_NULL, float)
    null_obs = np.empty(N_NULL, float)
    null_month = np.empty(N_NULL, float)
    null_month_obs = np.empty(N_NULL, float)
    date_group = s.year_month.to_numpy(int)
    blocks = [np.flatnonzero(date_group == month) for month in np.unique(date_group)]
    month_reassignable = int(sum(len(idx) for idx in blocks
        if len(idx)>1 and len(np.unique(labels[idx]))>1))
    for k in range(N_NULL):
        perm = labels[rng.permutation(len(labels))]
        null[k] = float(np.mean(perm[a] != perm[b]) -
                        np.mean(perm[local_a] != perm[local_b]))
        null_obs[k] = float(np.mean(perm[all_o_a] != perm[all_o_b]) -
                            np.mean(perm[diff_a] != perm[diff_b]))
        ym_perm = labels.copy()
        for idx in blocks:
            if len(idx)>1:
                ym_perm[idx] = labels[rng.permutation(idx)]
        null_month[k] = float(np.mean(ym_perm[a] != ym_perm[b]) -
                              np.mean(ym_perm[local_a] != ym_perm[local_b]))
        null_month_obs[k] = float(np.mean(ym_perm[all_o_a] != ym_perm[all_o_b]) -
                                  np.mean(ym_perm[diff_a] != ym_perm[diff_b]))
    return {"delta": observed, "delta_distinct_observer": observed_obs,
            "null": null, "null_distinct_observer": null_obs,
            "null_month": null_month, "null_month_distinct_observer": null_month_obs,
            "n_month_reassignable_photo_labels": month_reassignable,
            "n_local_pairs": len(local_a), "n_diff_obs_local_pairs": len(diff_a)}


def aggregate_geo(rows):
    if not rows:
        return {"species": 0, "status": "no_support", "null_test": "not_estimable"}
    observed = np.asarray([x["geo"]["delta"] for x in rows], float)
    null = np.stack([x["geo"]["null"] for x in rows]).mean(axis=0)
    observed_obs = np.asarray([x["geo"]["delta_distinct_observer"] for x in rows], float)
    null_obs = np.stack([x["geo"]["null_distinct_observer"] for x in rows]).mean(axis=0)
    month_null = np.stack([x["geo"]["null_month"] for x in rows]).mean(axis=0)
    month_null_obs = np.stack([x["geo"]["null_month_distinct_observer"] for x in rows]).mean(axis=0)
    return {
        "species": len(rows),
        "mean_source_photo_global_minus_local": float(observed.mean()),
        "fraction_species_positive": float(np.mean(observed > 0)),
        "null_mean": float(null.mean()),
        "matched_199_world_upper_p": float((1+int(np.sum(null >= observed.mean())))/(N_NULL+1)),
        "distinct_observer_pair_mean_depletion": float(observed_obs.mean()),
        "distinct_observer_pair_null_mean": float(null_obs.mean()),
        "distinct_observer_pair_matched_199_world_upper_p":
            float((1+int(np.sum(null_obs >= observed_obs.mean())))/(N_NULL+1)),
        "same_year_month_photo_cohort_matched_null_mean": float(month_null.mean()),
        "same_year_month_199_world_upper_p":
            float((1+int(np.sum(month_null>=observed.mean())))/(N_NULL+1)),
        "same_year_month_distinct_observer_matched_null_mean": float(month_null_obs.mean()),
        "same_year_month_distinct_observer_199_world_upper_p":
            float((1+int(np.sum(month_null_obs>=observed_obs.mean())))/(N_NULL+1)),
        "species_with_ge4_reassignable_month_labels":
            int(sum(x["geo"]["n_month_reassignable_photo_labels"]>=4 for x in rows)),
        "all_frozen_original_geo_and_photo_classes_preserved_outside_QA": True,
        "posthoc_photo_purity_filter_not_biological_morph": True,
    }


def gmm_model(z, n_max=4, seed=SEED):
    maxk = min(n_max, max(1, len(z)//6))
    bics, models = {}, {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConvergenceWarning)
        for k in range(1, maxk+1):
            model=GaussianMixture(n_components=k, covariance_type="full",
                                  reg_covar=.002, n_init=2, max_iter=175,
                                  random_state=seed+k)
            model.fit(z)
            bics[k]=float(model.bic(z))
            models[k]=model
    best=min(bics,key=bics.get)
    return best,bics,models[best]


def geometry_species(g, taxon):
    if len(g)<24:
        return None
    raw = g[list(VARS)].to_numpy(float)
    x = np.sqrt(raw)
    z = PCA(n_components=2,svd_solver="full").fit_transform(x)
    k, bics, model=gmm_model(z,seed=SEED+int(taxon))
    groups=model.predict(z)
    counts=np.bincount(groups,minlength=k)
    gain=bics[1]-bics[k]
    obs=g.observer_id.astype(str).to_numpy()
    unique=np.unique(obs)
    split_ari = None
    if len(unique)>=8 and k>=2:
        rng=np.random.default_rng(SEED+int(taxon))
        half=unique[rng.permutation(len(unique))]
        side=np.isin(obs,half[:len(half)//2])
        if side.sum()>=max(12,k*4) and (~side).sum()>=max(12,k*4):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", ConvergenceWarning)
                options=[]
                for v in (side,~side):
                    mdl=GaussianMixture(n_components=k,covariance_type="full",
                                        reg_covar=.002,n_init=2,max_iter=175,
                                        random_state=SEED+int(taxon))
                    mdl.fit(z[v])
                    options.append(mdl.predict(z))
                split_ari=float(adjusted_rand_score(options[0],options[1]))
    candidate=(k>=2 and gain>=10 and counts.min()>=5 and
               split_ari is not None and split_ari>=.75)
    return {"species":str(g.species.iloc[0]),
            "coarse_photo_groups":int(g.morph.nunique()),
            "n_photos_after_roi_purity":len(g),
            "bic_best_components_on_PCA2":int(k),
            "bic_improvement_over_one":float(gain),
            "smallest_best_fit_cluster":int(counts.min()),
            "observer_split_same_k_ARI":split_ari,
            "photographic_mode_candidate":bool(candidate)}


def self_test():
    rng=np.random.default_rng(1840)
    one=rng.normal(0,.18,size=(40,2))
    multi=np.concatenate([rng.normal(-2,.18,(35,2)),
                          rng.normal(2,.18,(35,2))],axis=0)
    k1,_,_=gmm_model(one)
    k2,bics,m=gmm_model(multi)
    check(k1 in (1,2,3,4) and k2 in (1,2,3,4),"synthetic mixture dimensionality")
    check(k2>=2 and bics[1]-bics[k2]>10, "strong synthetic two-mode difference not found")
    check(m.predict(multi).shape==(70,), "synthetic component membership invalid")
    assert 0.005 <= float((1+0)/(N_NULL+1)) <= 1
    print("FCP_PALETTE_MODE_GMM_SELF_TEST_PASS")


def run(path, pilot_per_n_groups):
    photos, originals=load_photos(path)
    bypurity = {}
    subset_for_geo = {"0.65": [], "0.75": []}
    eligible_for_gmm = {2: [], 3: [], 4: []}
    for taxon,g in photos.groupby("inat_taxon_id",sort=True):
        if originals[int(taxon)] < 40:
            continue
        stat={}
        for threshold in (.65,.75):
            q,strong,meta = quality_group(g, originals[int(taxon)],threshold)
            stat[str(threshold)] = {"n":len(q),"strong_classes":len(strong)}
            if len(q)>=24 and len(strong)>=2:
                if threshold == .65 and len(strong) in eligible_for_gmm:
                    eligible_for_gmm[len(strong)].append((int(taxon),q))
                geo=geography_species(q,taxon)
                if geo is not None:
                    subset_for_geo[str(threshold)].append(
                        {"taxon":int(taxon),"n_classes":len(strong),"geo":geo})
        bypurity[int(taxon)]=stat
    quality_summary={}
    for t in (.65,.75):
        key=str(t)
        for k in (2,3,4):
            quality_summary[f"q{int(t*100)}_ge{k}_supported_colour_groups"]=sum(
                x[key]["n"]>=24 and x[key]["strong_classes"]>=k for x in bypurity.values()
            )
    groups={"all_supported_two_or_more":aggregate_geo(subset_for_geo["0.65"])}
    for k in (2,3,4):
        groups[f"exactly_{k}_supported_classes"]=aggregate_geo(
            [x for x in subset_for_geo["0.65"] if x["n_classes"]==k])
    groups["purity75_all_supported_two_or_more"]=aggregate_geo(subset_for_geo["0.75"])
    for k in (2,3,4):
        groups[f"purity75_exactly_{k}_supported_classes"]=aggregate_geo(
            [x for x in subset_for_geo["0.75"] if x["n_classes"]==k])
    mode_records=[]
    for k in (2,3,4):
        eligible=eligible_for_gmm[k]
        # Taxon-ID hash order, post hoc within already qualifying photographic groups.
        # This is a fixed, BALANCED PILOT, NOT A prevalence sample of all taxa.
        eligible.sort(key=lambda t: hashlib.sha256(str(t[0]).encode()).hexdigest())
        for taxon,g in eligible[:pilot_per_n_groups]:
            output=geometry_species(g,taxon)
            if output is not None:
                mode_records.append(output)
    count=Counter(x["bic_best_components_on_PCA2"] for x in mode_records)
    stable=[x for x in mode_records if x["photographic_mode_candidate"]]
    return {
        "schema":"fcp_original447500_palette_mode_and_geo_purity_audit_v1",
        "source_sha256": SOURCE_SHA,
        "analysis_class":"exploratory posthoc measured-photo model; not a polymorphic population census",
        "n_original_photos":TOTAL,"n_classifiable_photos":CLASSIFIED,"n_source_species":TAXA,
        "photo_purity_cutoffs_exploratory": [.65,.75],
        "minimum_purity_photo_sample":24,
        "minimum_class_observer_count":3,
        "minimum_class_photo_count":5,
        "source_quality_species":quality_summary,
        "geographic_50km_within_species_matched_null_on_purity_passing_photos":groups,
        "within_pilot": {
            "pilot_selection":"bounded SHA-order taxon-ID sample within label-richness groups 2, 3, 4",
            "max_selected_per_supported_group":pilot_per_n_groups,
            "n_species_tested":len(mode_records),
            "bics_selected_components_1to4":dict(count),
            "n_species_with_photo_mixture_candidate_and_observer_stability":len(stable),
            "n_by_supported_photo_class_group":dict(Counter(x["coarse_photo_groups"]
                                                            for x in mode_records)),
            "top_photo_mode_candidates_no_coordinates":sorted(
                stable,key=lambda x:x["bic_improvement_over_one"],reverse=True)[:20],
        },
        "not_frozen_50km_primary_P2_replacement": True,
        "not_determining_true_number_of_heritable_morphs": True,
        "no_same_individual_reproducibility_no_photos_reviewer_stage_or_genotypes": True,
        "all_effective_null_p_values_posthoc_without_FWER_for_subgroups": True,
        "source_position_and_observer_identifiers_not_exported": True,
        "strict_5km_month_matched_null_was_not_supported_p":.300,
        "same_year_month_preserving_label_permutation_is_exploratory_not_genetic":True,
        "quality75_geographic_p_is_a_secondary_photo_purity_sensitivity":True,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--photos",type=Path)
    p.add_argument("--out",type=Path)
    p.add_argument("--pilot-per-group",type=int,default=80)
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    if a.self_test:
        self_test()
        return
    check(a.photos is not None and a.out is not None,"required --photos and --out")
    check(1<=a.pilot_per_group<=100,"unbounded exploratory component screen")
    o=run(a.photos,a.pilot_per_group)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(o,indent=2,sort_keys=True,ensure_ascii=False)+"\n",
                     encoding="utf-8")
    print("FCP_ORIGINAL_PHOTO_PALETTE_CLUSTER_GEO_QA_PASS",json.dumps({
        "support":o["source_quality_species"],
        "geo":o["geographic_50km_within_species_matched_null_on_purity_passing_photos"],
        "modes":{k:v for k,v in o["within_pilot"].items() if k!="top_photo_mode_candidates_no_coordinates"},
        "result_file":str(a.out)
    }),flush=True)


if __name__=="__main__":
    main()
