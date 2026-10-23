#!/usr/bin/env python3
"""Colour-specific 50-km photographic spatial autocorrelation with matched nulls.

Ecological question: for each coarse colour label separately, are two nearby
conspecific photos more often BOTH that colour than expected after fixing the
species' four-state totals, geographic vertices and (for the seasonal null)
its recorded calendar-month colour composition?

This is a post hoc photo-label analysis, NOT DBSCAN, genetic morph evidence,
ecological selection, native-range validation or a true flower-colour confusion
matrix. Never infer exact locality, observer identity or private photo IDs in
published output.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

PHOTO_SHA = "8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6"
COLOURS = ("white", "yellow_orange", "red_pink", "blue_purple")
REPS = 199
KM = 50.
R_EARTH = 6371.0088
SOURCE_ROWS = 447500
SOURCE_CLASSIFIED = 227107
SOURCE_TAXA = 4475
PRIMARY_STRICT_TAXA = 1761
PRIMARY_STRICT_DELTA = 0.0176149021685411
USECOLS = (
    "inat_taxon_id", "photo_id", "observer_id", "observed_on",
    "morph", "global_classifiable", "latitude", "longitude",
    "positional_accuracy_m",
)


def check(x, msg):
    if not x:
        raise ValueError("FCP_COLOUR_CLUSTER_AUDIT: " + msg)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as fp:
        for chunk in iter(lambda: fp.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_source(path, authentic=True):
    if authentic:
        check(sha(path) == PHOTO_SHA, "source 447500-photo checksum changed")
    t = pd.read_csv(path, usecols=USECOLS, low_memory=False)
    if authentic:
        check(len(t) == SOURCE_ROWS, "source original-photo row count drift")
        check(t.photo_id.nunique() == SOURCE_ROWS, "duplicate original photo IDs")
        check(t.inat_taxon_id.nunique() == SOURCE_TAXA, "source taxon frame drift")
    flag = t.global_classifiable.fillna(False).astype(str).str.lower().isin(("1", "true", "yes"))
    t = t.loc[flag].copy()
    check(t.morph.isin(COLOURS).all(), "unexpected coarse colour assignment")
    if authentic:
        check(len(t) == SOURCE_CLASSIFIED, "source four-state count changed")
        check(t.inat_taxon_id.nunique() == SOURCE_TAXA, "taxa dropped when classified")
    for k in ("latitude", "longitude", "positional_accuracy_m"):
        t[k] = pd.to_numeric(t[k], errors="coerce")
    t["month"] = pd.to_datetime(t.observed_on, errors="coerce").dt.month.fillna(0).astype("int8")
    return t


def location_graph(group):
    good = (group.latitude.between(-90, 90) & group.longitude.between(-180, 180)
            & group.positional_accuracy_m.between(0, 5000))
    g = group.loc[good]
    if len(g) < 2:
        return None
    la = np.deg2rad(g.latitude.to_numpy(float))
    lo = np.deg2rad(g.longitude.to_numpy(float))
    a, b = np.triu_indices(len(g), 1)
    z = np.sin((la[a] - la[b]) / 2) ** 2
    z += np.cos(la[a]) * np.cos(la[b]) * np.sin((lo[a] - lo[b]) / 2) ** 2
    d = 2 * R_EARTH * np.arcsin(np.sqrt(np.clip(z, 0, 1)))
    near = d <= KM
    if not near.any():
        return None
    y = pd.Categorical(g.morph, categories=COLOURS).codes.astype(np.int8)
    check(np.all(y >= 0), "invalid categorical code")
    months = g.month.to_numpy(np.int8)
    return y, a[near], b[near], int(near.sum()), months


def same_colour_pair_fraction(y, i, j):
    equal = y[i] == y[j]
    return np.bincount(y[i][equal], minlength=4).astype(float) / len(i)


def month_groups(months):
    return [np.flatnonzero(months == m) for m in np.unique(months) if np.count_nonzero(months == m) > 1]


def one_species(group, taxon_id, *, one_observer=False):
    original_n = len(group)
    if one_observer:
        group = (group.dropna(subset=["observer_id"])
                 .sort_values(["observer_id", "photo_id"], kind="mergesort")
                 .drop_duplicates("observer_id").sort_index(kind="mergesort"))
    if original_n < 40 or len(group) < 40:
        return None
    graph = location_graph(group)
    if graph is None:
        return None
    y, i, j, n_pairs, months = graph
    if n_pairs < 30:
        return None
    counts = np.bincount(y, minlength=4)
    admitted = (counts >= 5) & ((len(y) - counts) >= 5)
    observed = same_colour_pair_fraction(y, i, j)
    prob_global_same = counts * (counts - 1) / (len(y) * (len(y) - 1))
    all_mismatch = float(1 - np.sum(prob_global_same))
    observed_mismatch = float(1 - np.sum(observed))
    groups = month_groups(months)
    month_swaps_possible = any(len(np.unique(y[z])) > 1 for z in groups)
    random = np.random.default_rng(20261023 + int(taxon_id) * 17 + int(one_observer))
    standard = np.empty((REPS, 4), float)
    seasonal = np.empty((REPS, 4), float)
    for k in range(REPS):
        p = y[random.permutation(len(y))]
        standard[k] = same_colour_pair_fraction(p, i, j)
        q = y.copy()
        for z in groups:
            q[z] = random.permutation(y[z])
        seasonal[k] = same_colour_pair_fraction(q, i, j)
    return {
        "admitted": admitted, "observed": observed, "standard": standard,
        "seasonal": seasonal, "month_swaps_possible": month_swaps_possible,
        "mismatch_delta": all_mismatch - observed_mismatch,
        "n_local_pairs": n_pairs, "n_classifiable": original_n,
        "n_valid": len(y),
    }


def holm_by_family(values):
    order = np.argsort(values, kind="stable")
    out = np.zeros(len(values), float)
    current = 0.
    for rank, k in enumerate(order):
        current = max(current, min(1., (len(values) - rank) * float(values[k])))
        out[k] = current
    return out.tolist()


def evaluate_taxa(t, *, one_observer=False):
    n_species_strict = 0
    original_strict_effects = []
    by_col = [list() for _ in COLOURS]
    count_month_swappable = np.zeros(4, int)
    for ix, (tid, group) in enumerate(t.groupby("inat_taxon_id", sort=True)):
        # For source-reproduction, the strict denominator uses pre-thinning photos.
        if not one_observer:
            original_graph = location_graph(group)
            if len(group) >= 40 and original_graph is not None and original_graph[3] >= 30:
                n_species_strict += 1
                y, i, j, *_ = original_graph
                n=len(y)
                c=np.bincount(y, minlength=4)
                global_diff=float(1-np.sum(c*(c-1))/(n*(n-1)))
                local_diff=float(np.mean(y[i]!=y[j]))
                original_strict_effects.append(global_diff-local_diff)
        computed = one_species(group, tid, one_observer=one_observer)
        if computed is None:
            continue
        for col in range(4):
            if computed["admitted"][col]:
                by_col[col].append(computed)
                count_month_swappable[col] += int(computed["month_swaps_possible"])
        if (ix + 1) % 1000 == 0:
            print("FCP_COLOUR_CLUSTER_TAXA", "observer_thinned" if one_observer else "original",
                  ix+1, flush=True)
    if not one_observer:
        check(n_species_strict == PRIMARY_STRICT_TAXA,
              f"original strict 50km n changed {n_species_strict}")
        check(abs(np.mean(original_strict_effects) - PRIMARY_STRICT_DELTA) < 1e-10,
              "historical 50km result differs from source-locked value")

    family = {m: [] for m in ("global_composition", "calendar_month_composition")}
    raw_p = {m: [] for m in family}
    for colour, observed_list in zip(COLOURS, by_col):
        n=len(observed_list)
        if n == 0:
            for label in family:
                family[label].append({
                    "colour":colour, "n_species":0, "p_upper_199":None,
                    "mean_excess":None, "mean_observed_same_pair_fraction":None,
                })
                raw_p[label].append(1.)
            continue
        obs=np.array([x["observed"][COLOURS.index(colour)] for x in observed_list])
        for label, source in (("global_composition", "standard"),
                              ("calendar_month_composition", "seasonal")):
            worlds=np.mean(np.stack([x[source][:,COLOURS.index(colour)] for x in observed_list]), axis=0)
            delta=float(obs.mean()-worlds.mean())
            p=float((1+np.count_nonzero(worlds>=obs.mean()))/(REPS+1))
            raw_p[label].append(p)
            family[label].append({
                "colour":colour, "n_species":n,
                "n_species_with_nontrivial_month_colour_swaps":int(count_month_swappable[COLOURS.index(colour)]),
                "mean_observed_same_pair_fraction":float(obs.mean()),
                "mean_matched_null_fraction":float(worlds.mean()),
                "mean_excess":delta,
                "p_upper_199":p,
                "null_world_q025":float(np.quantile(worlds,.025)),
                "null_world_q975":float(np.quantile(worlds,.975)),
                "species_equal":True,
                "support_gate":"n_classifiable>=40; n_geo_valid with >=30 near pairs; colour>=5, other_colour>=5",
            })
    for label, rows in family.items():
        adjusted=holm_by_family(raw_p[label])
        for k, row in enumerate(rows):
            row["holm_p_four_colours"] = adjusted[k] if row["n_species"]>0 else None
    return {
        "n_original_strict_species_reproduced":n_species_strict if not one_observer else None,
        "original_strict_mean_reproduced":float(np.mean(original_strict_effects)) if not one_observer else None,
        "colour_specific_null_families":family,
        "observer_thinned_colour_blind_by_photo_id":one_observer,
    }


def self_test():
    # Geographic clustering of both colours; no year-month morph validation inferred.
    rows=[]
    for tid in (31, 32):
        for label in range(2):
            for i in range(24):
                rows.append({
                    "inat_taxon_id":tid, "photo_id":tid*1000+label*30+i,
                    "observer_id":tid*1000+label*30+i, "observed_on":"2024-06-01",
                    "morph":COLOURS[label], "global_classifiable":True,
                    "latitude":35. if not label else 38.,
                    "longitude":135. + .005*(i%5),
                    "positional_accuracy_m":10., "month":6,
                })
    d=pd.DataFrame(rows)
    r=one_species(d[d.inat_taxon_id==31],31)
    check(r is not None and r["n_local_pairs"]>=30,"synthetic opportunity lost")
    check(all(r["admitted"][:2]),"synthetic photo-colour supports lost")
    # Calendar-month null is permissive here because both morphs share June.
    check(r["observed"][0]>r["standard"][:,0].mean(),"synthetic white clustering not recovered")
    check(r["observed"][1]>r["seasonal"][:,1].mean(),"synthetic month-null clustering lost")
    d.loc[d.morph=="blue_purple","morph"]="yellow_orange" if False else "blue_purple"
    # A colour label segregated *only* by month is fixed under calendar-month
    # conditioning and cannot be called a month-adjusted spatial cluster.
    d.loc[d.morph==COLOURS[0], "month"]=5
    d.loc[d.morph==COLOURS[1], "month"]=6
    rr=one_species(d[d.inat_taxon_id==31],31)
    check(np.allclose(rr["seasonal"][:,0],rr["observed"][0]),
          "month-conditioned null failed to block phenology-confounded cluster")
    return True


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--photos",type=Path)
    parser.add_argument("--out",type=Path)
    parser.add_argument("--self-test",action="store_true")
    z=parser.parse_args()
    if z.self_test:
        print("FCP_COLOUR_CLUSTER_SELF_TEST_PASS" if self_test() else "FCP_COLOUR_CLUSTER_FAIL")
        return
    check(z.photos and z.out,"photos and output location required")
    photos=read_source(z.photos,authentic=True)
    out={
        "schema":"fcp_photo_colour_specific_50km_matched_permutation_and_month_v1",
        "date_jst":"2026-10-23",
        "original_source_photo_sha256":PHOTO_SHA,
        "n_source_taxa":SOURCE_TAXA,
        "n_source_original_photos":SOURCE_ROWS,
        "n_classified_source_photos":SOURCE_CLASSIFIED,
        "same_photo_provider":True,
        "posthoc":True,
        "permutations":REPS,
        "radius_km":KM,
        "colour_order":list(COLOURS),
        "description":"Both-local-photos-of-one-colour excess by taxon, not estimated discrete clusters",
        "source_null":"shuffle original labelled photo vertices within species holding all four counts fixed",
        "seasonal_null":"shuffle original labelled vertices within each species and calendar month",
        "important_limit":"Not an individual-petal validation, true genotypic morph cluster, selection or native-range audit",
        "original":evaluate_taxa(photos,one_observer=False),
        "one_per_observer":evaluate_taxa(photos,one_observer=True),
    }
    z.out.mkdir(parents=True,exist_ok=True)
    (z.out/"fcp_colour_specific_50km_month_control_results_20261023.json").write_text(
        json.dumps(out,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("FCP_COLOUR_CLUSTER_RESULTS",json.dumps(
        {mode:{
            t:[(x["colour"],x["n_species"],x["mean_excess"],x["holm_p_four_colours"])
               for x in data["colour_specific_null_families"][t]]
            for t in ("global_composition","calendar_month_composition")
        } for mode,data in (("original",out["original"]),("one_per_observer",out["one_per_observer"]))}
    ),flush=True)


if __name__=="__main__":
    main()
