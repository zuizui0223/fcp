#!/usr/bin/env python3
"""Source-locked robustness audit of photographed flower-colour labels.

Two stages: permanent 4,475-taxon four-colour counts, then authenticated
447,500-photo observer/observation evidence and matched spatial permutations.
Neither stage identifies inherited morphs, image errors or selection.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

LEDGER_SHA = "23d3daadc9799c2e7b5867d54ccc1ed98369ba79a81b72dc7f831e5197e63c49"
PHOTO_SHA = "8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6"
COLOURS = ("white", "yellow_orange", "red_pink", "blue_purple")
COLUMNS = tuple("miss_coarse_photo_" + c for c in COLOURS)
PERMS = 199
SEED = 20261021


def require(ok, msg):
    if not ok:
        raise RuntimeError("FCP_COLOUR_LABEL_AUDIT: " + msg)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def audit_ledger(d):
    counts = d[list(COLUMNS)].to_numpy(dtype=int)
    n = d.n_classified_photo.to_numpy(dtype=int)
    require(len(d) == 4475 and d.inat_taxon_id.nunique() == 4475, "wrong taxon frame")
    require((counts.sum(axis=1) == n).all(), "colour counts inconsistent")
    require(d.n_source_photo.eq(100).all() and n.sum() == 227107, "photo count drift")
    summary = {
        "n_species": int(len(d)),
        "n_original_photo_rows": 447500,
        "n_four_colour_classifiable_photo_rows": int(n.sum()),
        "n_unclassified_photo_rows": int(447500 - n.sum()),
        "observed_n_colour_classes": {
            str(k): int(((counts > 0).sum(axis=1) == k).sum())
            for k in (1, 2, 3, 4)
        },
        "minimum_misclassified_photos_required_to_be_monochrome": {
            str(k): int((n - counts.max(axis=1) > k).sum())
            for k in (0, 1, 2, 3, 5, 10, 15, 20)
        },
        "colour_photo_count_tiers": {},
        "nclass_bins": {},
    }
    for min_n in (8, 20, 40, 60, 80):
        t = {"n_species": int((n >= min_n).sum())}
        for photos_per_colour in (1, 2, 3, 5, 10):
            eligible = (counts >= photos_per_colour).sum(axis=1)
            for k in (2, 3, 4):
                t[f"at_least_{k}_colours_each_ge_{photos_per_colour}"] = int(
                    ((n >= min_n) & (eligible >= k)).sum()
                )
        summary["colour_photo_count_tiers"][str(min_n)] = t
    for lo, hi in ((8, 19), (20, 39), (40, 59), (60, 79), (80, 100)):
        idx = (n >= lo) & (n <= hi)
        summary["nclass_bins"][f"{lo}-{hi}"] = {
            "n_species": int(idx.sum()),
            "at_least_2_colours_1_photo": int((idx & ((counts > 0).sum(axis=1) >= 2)).sum()),
            "at_least_2_colours_5_photos": int((idx & ((counts >= 5).sum(axis=1) >= 2)).sum()),
            "at_least_3_colours_1_photo": int((idx & ((counts > 0).sum(axis=1) >= 3)).sum()),
            "at_least_4_colours_1_photo": int((idx & ((counts > 0).sum(axis=1) >= 4)).sum()),
        }
    require(summary["observed_n_colour_classes"] ==
            {"1": 286, "2": 1743, "3": 1632, "4": 814},
            "permanent source photograph-colour classification drift")
    require(summary["colour_photo_count_tiers"]["8"]["at_least_2_colours_each_ge_5"] == 2400,
            "two-photo-colour support threshold drift")
    require(summary["colour_photo_count_tiers"]["8"]["at_least_3_colours_each_ge_5"] == 435,
            "three-photo-colour support threshold drift")
    require(summary["colour_photo_count_tiers"]["8"]["at_least_4_colours_each_ge_5"] == 19,
            "four-photo-colour support threshold drift")
    return summary


def haversine(lat, lon, a, b):
    la = np.deg2rad(lat.astype(float))
    lo = np.deg2rad(lon.astype(float))
    q = np.sin((la[a] - la[b]) / 2)**2
    q += np.cos(la[a]) * np.cos(la[b]) * np.sin((lo[a] - lo[b]) / 2)**2
    return 2 * 6371.0088 * np.arcsin(np.sqrt(np.clip(q, 0, 1)))


def distance_null(g, taxon, seed):
    v = (g.latitude.between(-90, 90) & g.longitude.between(-180, 180) &
         g.positional_accuracy_m.between(0, 5000))
    g = g.loc[v]
    n = len(g)
    if n < 2:
        return None
    a, b = np.triu_indices(n, 1)
    dist = haversine(g.latitude.to_numpy(float), g.longitude.to_numpy(float), a, b)
    within = dist <= 50
    if not within.any():
        return None
    y = pd.Categorical(g.morph, categories=COLOURS).codes.astype(np.int8)
    i, j = a[within], b[within]
    all_mismatch = float((y[a] != y[b]).mean())
    local_mismatch = float((y[i] != y[j]).mean())
    null = np.empty(PERMS, dtype=float)
    rng = np.random.default_rng(seed + int(taxon))
    for k in range(PERMS):
        p = y[rng.permutation(n)]
        null[k] = all_mismatch - float((p[i] != p[j]).mean())
    return all_mismatch - local_mismatch, null, int(within.sum())


def matched_null(items):
    if not items:
        return {"n_species": 0, "mean": None, "upper_p_199": None}
    observed = np.asarray([x[0] for x in items], dtype=float)
    null = np.stack([x[1] for x in items]).mean(axis=0)
    return {
        "n_species": len(items),
        "mean_delta": float(observed.mean()),
        "median_delta": float(np.median(observed)),
        "fraction_positive": float((observed > 0).mean()),
        "mean_null": float(null.mean()),
        "upper_p_199": float((1 + np.count_nonzero(null >= observed.mean())) / (PERMS + 1)),
        "permutations": PERMS,
    }


def audit_original_photos(path, ledger):
    require(sha256(path) == PHOTO_SHA, "authenticated 447500-photo SHA mismatch")
    cols = ["inat_taxon_id", "photo_id", "observation_id", "observer_id",
            "latitude", "longitude", "positional_accuracy_m", "global_classifiable", "morph"]
    t = pd.read_csv(path, usecols=cols, low_memory=False)
    require(len(t) == 447500 and t.photo_id.nunique() == 447500, "photo IDs or count drift")
    require(t.inat_taxon_id.nunique() == 4475, "taxon count drift")
    flag = t.global_classifiable.fillna(False).astype(str).str.lower().isin(("true", "1", "yes"))
    require(flag.sum() == 227107, "four-colour photo denominator drift")
    d = t.loc[flag].copy()
    require(d.morph.isin(COLOURS).all(), "unrecognized photograph colour state")
    for name in ("latitude", "longitude", "positional_accuracy_m"):
        d[name] = pd.to_numeric(d[name], errors="coerce")
    idx = ledger.set_index("inat_taxon_id")
    c = (d.groupby(["inat_taxon_id", "morph"], observed=True).size()
           .unstack(fill_value=0).reindex(index=idx.index, columns=COLOURS, fill_value=0))
    require(np.array_equal(c.to_numpy(int), idx[list(COLUMNS)].to_numpy(int)),
            "per-photograph colour counts disagree with permanent species ledger")
    observers = (d.groupby(["inat_taxon_id", "morph"], observed=True).observer_id.nunique()
                 .unstack(fill_value=0).reindex(index=idx.index, columns=COLOURS, fill_value=0))
    observations = (d.groupby(["inat_taxon_id", "morph"], observed=True).observation_id.nunique()
                    .unstack(fill_value=0).reindex(index=idx.index, columns=COLOURS, fill_value=0))
    max_per_observer = int(t.dropna(subset=["observer_id"])
                           .groupby(["inat_taxon_id", "observer_id"]).size().max())
    require(max_per_observer <= 2, "source sampling observer cap violated")
    count5 = c.to_numpy() >= 5
    supported = count5 & (observers.to_numpy() >= 3) & (observations.to_numpy() >= 3)
    g = d.dropna(subset=["observation_id"]).groupby(["inat_taxon_id", "observation_id"])
    same_obs = g.agg(n=("morph", "size"), k=("morph", "nunique"))
    duplicate = same_obs.n >= 2
    conflict = duplicate & (same_obs.k >= 2)
    out = {
        "n_source_classifiable_photos": len(d),
        "max_original_photos_per_observer_per_taxon": max_per_observer,
        "n_species_two_colours_each_ge5_photos": int((count5.sum(axis=1) >= 2).sum()),
        "n_species_two_colours_each_ge5_photos_and_3_observers_and_3_observations":
            int((supported.sum(axis=1) >= 2).sum()),
        "n_species_three_colours_each_ge5_photos_and_3_observers_and_3_observations":
            int((supported.sum(axis=1) >= 3).sum()),
        "n_species_four_colours_each_ge5_photos_and_3_observers_and_3_observations":
            int((supported.sum(axis=1) >= 4).sum()),
        "same_observation_multi_photo_comparisons": int(duplicate.sum()),
        "same_observation_conflicting_colour_label_groups": int(conflict.sum()),
        "taxa_with_same_observation_conflicting_colour": int(
            same_obs.loc[conflict].index.get_level_values(0).nunique()),
    }
    tiers = {name: [] for name in (
        "full_original", "photo_5_each_n40_pairs30", "photo_10_each_n40_pairs30",
        "photo_5_each_3types_n40_pairs30", "one_per_observer",
        "one_per_observer_photo_5_each_n40_pairs30",
        "one_per_observer_photo_10_each_n40_pairs30")}
    per_species = []
    for k, (tid, group) in enumerate(d.groupby("inat_taxon_id", sort=True)):
        tid = int(tid)
        original = distance_null(group, tid, SEED)
        expected = idx.loc[tid, "global_vs_local_50km_colour_mismatch_delta"]
        if original is None:
            require(pd.isna(expected), "source ledger local-pair NA inconsistency")
        else:
            require(np.isfinite(expected) and abs(expected - original[0]) < 1e-9,
                    "original spatial delta mismatch for taxon " + str(tid))
            tiers["full_original"].append(original)
            source_counts = c.loc[tid].to_numpy(int)
            if len(group) >= 40 and original[2] >= 30:
                if (source_counts >= 5).sum() >= 2:
                    tiers["photo_5_each_n40_pairs30"].append(original)
                if (source_counts >= 10).sum() >= 2:
                    tiers["photo_10_each_n40_pairs30"].append(original)
                if (source_counts >= 5).sum() >= 3:
                    tiers["photo_5_each_3types_n40_pairs30"].append(original)
        # Photo selection depends on observer and photo ID, NOT on measured colour.
        independent = (group.dropna(subset=["observer_id"])
                       .sort_values(["observer_id", "photo_id"], kind="mergesort")
                       .drop_duplicates("observer_id").sort_index(kind="mergesort"))
        indep_effect = distance_null(independent, tid, 20261023)
        if indep_effect is not None:
            tiers["one_per_observer"].append(indep_effect)
            independent_counts = independent.morph.value_counts().reindex(COLOURS, fill_value=0)
            if len(independent) >= 40 and indep_effect[2] >= 30:
                if (independent_counts >= 5).sum() >= 2:
                    tiers["one_per_observer_photo_5_each_n40_pairs30"].append(indep_effect)
                if (independent_counts >= 10).sum() >= 2:
                    tiers["one_per_observer_photo_10_each_n40_pairs30"].append(indep_effect)
        per_species.append({
            "inat_taxon_id": tid,
            "n_classified_photos": int(len(group)),
            "n_photo_colours": int((c.loc[tid] > 0).sum()),
            "n_photo_colours_supported_5photos_3observers_3observations":
                int(supported[idx.index.get_loc(tid)].sum()),
            "n_one_photo_per_observer": int(len(independent)),
            "n_source_50km_pairs": original[2] if original else 0,
            "n_one_per_observer_50km_pairs": indep_effect[2] if indep_effect else 0,
        })
        if (k + 1) % 1000 == 0:
            print("FCP_COLOUR_LABEL_AUDIT_TAXA_COMPLETE", k + 1, flush=True)
    out["species_conditioned_spatial_null_sensitivity"] = {
        k: matched_null(v) for k, v in tiers.items()
    }
    primary = out["species_conditioned_spatial_null_sensitivity"]["full_original"]
    require(primary["n_species"] == 3600, "original 50-km analysis denominator drift")
    require(abs(primary["mean_delta"] - 0.016115883798023662) < 1e-10,
            "original 50-km effect drift")
    require(primary["upper_p_199"] == 0.005, "original permutation p drift")
    return out, pd.DataFrame(per_species)


def self_test():
    x = pd.DataFrame({
        "n_source_photo": [100, 100], "n_classified_photo": [4, 6],
        "miss_coarse_photo_white": [2, 5],
        "miss_coarse_photo_yellow_orange": [2, 1],
        "miss_coarse_photo_red_pink": [0, 0],
        "miss_coarse_photo_blue_purple": [0, 0],
        "inat_taxon_id": [1, 2],
    })
    # Audit ledger requires exact 4,475 source taxa; test core arithmetic independently.
    counts = x[list(COLUMNS)].to_numpy(int)
    require(np.array_equal(counts.sum(axis=1), x.n_classified_photo), "synthetic counts")
    q = pd.DataFrame({
        "morph": ["white", "blue_purple", "white"],
        "latitude": [0., 0., 25.], "longitude": [0., .02, 0.],
        "positional_accuracy_m": [100., 100., 100.],
    })
    z = distance_null(q, 1, SEED)
    require(z is not None and z[2] == 1 and abs(z[0] - (2 / 3 - 1)) < 1e-12,
            "synthetic local spatial mismatch wrong")
    print("FCP_COLOUR_LABEL_AUDIT_SYNTHETIC_PASS", flush=True)


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--species-ledger", type=Path,
        default=Path("results/fcp_np_all4475_unified_species_ledger_20261021.csv"))
    a.add_argument("--photos", type=Path)
    a.add_argument("--out", type=Path)
    a.add_argument("--self-test", action="store_true")
    z = a.parse_args()
    self_test()
    if z.self_test:
        return
    require(sha256(z.species_ledger) == LEDGER_SHA, "source species ledger SHA changed")
    ledger = pd.read_csv(z.species_ledger, low_memory=False)
    summary = {
        "schema": "fcp_photographic_colour_label_robustness_v1",
        "classification_and_spatial_inference_only_not_genetic_morphs": True,
        "exploratory_after_primary_results": True,
        "source_photo_sha256": PHOTO_SHA,
        "source_ledger_sha256": LEDGER_SHA,
        "species_counts": audit_ledger(ledger),
    }
    panel = None
    if z.photos:
        summary["original_photo_and_observer_check"], panel = audit_original_photos(z.photos, ledger)
    result = json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if z.out:
        z.out.mkdir(parents=True, exist_ok=True)
        (z.out / "fcp_photo_label_robustness_summary_20261023.json").write_text(result)
        if panel is not None:
            panel.to_csv(z.out / "fcp_observer_colour_support_taxa_NO_OBSERVER_IDS_20261023.csv",
                         index=False)
    print("FCP_COLOUR_LABEL_AUDIT_RESULT", json.dumps(summary, sort_keys=True,
          separators=(",", ":"), allow_nan=False), flush=True)


if __name__ == "__main__":
    main()
