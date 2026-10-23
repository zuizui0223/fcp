#!/usr/bin/env python3
"""Test colour-specific photographic geographic associations, not genetic morphs.

Works on the original SHA-verified 447,500-photo FCP source. Test units are
species, and observed photographs never leave the runner. Matched label shuffles
preserve within-species coarse-colour counts and geographic vertices. Four
colour-specific tests receive a within-family maxT adjustment. The comparison
against calendar-month-preserving shuffles is a *separate* temporal sensitivity,
not a pre-registered or source-independent confirmation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

PHOTO_SHA = "8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6"
LEDGER_SHA = "23d3daadc9799c2e7b5867d54ccc1ed98369ba79a81b72dc7f831e5197e63c49"
COLORS = ("white", "yellow_orange", "red_pink", "blue_purple")
FRACTIONS = {
    "white": ("flower_fraction_white",),
    "yellow_orange": ("flower_fraction_yellow", "flower_fraction_orange"),
    "red_pink": ("flower_fraction_red", "flower_fraction_pink", "flower_fraction_magenta"),
    "blue_purple": ("flower_fraction_blue", "flower_fraction_purple"),
}
BASE_COLUMNS = (
    "inat_taxon_id", "photo_id", "observer_id", "latitude", "longitude",
    "positional_accuracy_m", "observed_on", "global_classifiable", "morph",
)
COLUMNS = tuple(dict.fromkeys(
    (*BASE_COLUMNS, *(v for fields in FRACTIONS.values() for v in fields))
))
REPS = 199
SEED = 20261023
RADIUS_KM = 50.0
CASES = (
    ("two_colors_5photos", 5, False, False, False),
    ("one_photo_per_observer_5photos", 5, True, False, False),
    ("two_colors_10photos", 10, False, False, False),
    ("one_photo_per_observer_10photos", 10, True, False, False),
    ("two_colors_5photos_month_preserving_null", 5, False, False, True),
    ("dominant_pixel_fraction65_5photos", 5, False, True, False),
    ("one_observer_dominant_fraction65_5photos", 5, True, True, False),
)


def assert_truth(value, message):
    if not value:
        raise AssertionError("FCP_COLOUR_CLUSTER_AUDIT: " + message)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_source(path, ledger_path):
    assert_truth(sha256(path) == PHOTO_SHA, "original photographic source SHA changed")
    assert_truth(sha256(ledger_path) == LEDGER_SHA, "species summary SHA changed")
    d = pd.read_csv(path, usecols=list(COLUMNS), low_memory=False)
    ledger = pd.read_csv(ledger_path, usecols=[
        "inat_taxon_id", "n_classified_photo",
        *["miss_coarse_photo_" + c for c in COLORS]
    ])
    assert_truth(len(d) == 447500 and d.photo_id.nunique() == 447500,
                 "photo-row or photo-ID denominator")
    assert_truth(d.inat_taxon_id.nunique() == 4475, "source taxon population")
    classified = d.global_classifiable.fillna(False).astype(str).str.lower().isin(
        ("true", "yes", "1")
    )
    assert_truth(int(classified.sum()) == 227107, "four-colour support")
    d = d.loc[classified].copy()
    assert_truth(d.morph.isin(COLORS).all(), "original category mismatch")
    for x in ("latitude", "longitude", "positional_accuracy_m"):
        d[x] = pd.to_numeric(d[x], errors="coerce")
    d = d.loc[
        d.latitude.between(-90, 90) & d.longitude.between(-180, 180) &
        d.positional_accuracy_m.between(0, 5000)
    ].copy()
    # The original sampling filter already demanded coordinate accuracy <=5km.
    assert_truth(len(d) == 227107, "unexpected removal of original classified geo points")
    pd_counts = (d.groupby(["inat_taxon_id", "morph"], observed=True).size()
                   .unstack(fill_value=0)
                   .reindex(index=ledger.inat_taxon_id, columns=COLORS, fill_value=0))
    source_counts = ledger.set_index("inat_taxon_id")[
        ["miss_coarse_photo_" + c for c in COLORS]
    ].to_numpy(dtype=int)
    assert_truth(np.array_equal(pd_counts.to_numpy(dtype=int), source_counts),
                 "original photo labels disagree with validated species ledger")
    date = pd.to_datetime(d.observed_on, errors="coerce")
    d["month_of_year"] = date.dt.month.fillna(-1).astype("int8")
    fraction_names = []
    for c in COLORS:
        name = "colour_fraction_" + c
        d[name] = d[list(FRACTIONS[c])].sum(axis=1)
        fraction_names.append(name)
    f = d[fraction_names].to_numpy(dtype=float)
    y = pd.Categorical(d.morph, categories=COLORS).codes
    assert_truth((y >= 0).all(), "missing category integer")
    d["primary_colour_pixel_fraction"] = f[np.arange(len(d)), y]
    f[np.arange(len(d)), y] = -np.inf
    d["other_max_colour_pixel_fraction"] = f.max(axis=1)
    assert_truth(d.primary_colour_pixel_fraction.notna().all(),
                 "ROI colour fraction missing")
    return d


def haversine(lat, lon, a, b):
    lat = np.deg2rad(np.asarray(lat, dtype=float))
    lon = np.deg2rad(np.asarray(lon, dtype=float))
    q = np.sin((lat[a] - lat[b]) / 2) ** 2
    q += np.cos(lat[a]) * np.cos(lat[b]) * np.sin((lon[a] - lon[b]) / 2) ** 2
    return 2 * 6371.0088 * np.arcsin(np.sqrt(np.clip(q, 0, 1)))


def choose_observers(g):
    # No flower-hue labels enter the photo identity sorting step.
    return (g.dropna(subset=["observer_id"])
            .sort_values(["observer_id", "photo_id"], kind="mergesort")
            .drop_duplicates("observer_id")
            .sort_index(kind="mergesort"))


def prepare(g, min_each=5, one_per_observer=False, dominant=False):
    if dominant:
        g = g.loc[
            (g.primary_colour_pixel_fraction >= .65) &
            ((g.primary_colour_pixel_fraction -
              g.other_max_colour_pixel_fraction) >= .20)
        ]
    if one_per_observer:
        g = choose_observers(g)
    n = len(g)
    if n < 40:
        return None
    y = pd.Categorical(g.morph, categories=COLORS).codes.astype(np.int8)
    counts = np.bincount(y, minlength=4)
    if int((counts >= min_each).sum()) < 2:
        return None
    a, b = np.triu_indices(n, 1)
    within = haversine(g.latitude.to_numpy(), g.longitude.to_numpy(), a, b) <= RADIUS_KM
    if int(within.sum()) < 30:
        return None
    i, j = a[within], b[within]
    expected = counts * (counts - 1) / (n * (n - 1))
    observed = same_colour_frequency(y, i, j) - expected
    return {
        "labels": y, "counts": counts, "expected": expected,
        "i": i, "j": j, "months": g.month_of_year.to_numpy(),
        "observed": observed, "n_photos": n, "n_pairs": len(i),
        "richness": int((counts >= min_each).sum()),
    }


def same_colour_frequency(y, i, j):
    matches = y[i] == y[j]
    return np.bincount(y[i[matches]], minlength=4).astype(float) / len(i)


def permute_labels(y, months, rng, month_preserving):
    if not month_preserving:
        return y[rng.permutation(len(y))]
    out = y.copy()
    for m in np.unique(months):
        idx = np.flatnonzero(months == m)
        if len(idx) > 1:
            out[idx] = y[rng.permutation(idx)]
    return out


class Accumulator:
    def __init__(self):
        self.n_species = 0
        self.richness = {str(k): 0 for k in (2, 3, 4)}
        self.obs_all = np.zeros(4)
        self.null_all = np.zeros((REPS, 4))
        self.n_colours = np.zeros(4, dtype=int)
        self.obs_supported = np.zeros(4)
        self.null_supported = np.zeros((REPS, 4))
        self.positive_supported = np.zeros(4, dtype=int)
        self.month_movable_species = 0

    def add(self, panel, seed, min_each, month_preserving):
        self.n_species += 1
        self.richness[str(panel["richness"])] += 1
        y, i, j = panel["labels"], panel["i"], panel["j"]
        counts, expected = panel["counts"], panel["expected"]
        obs = panel["observed"]
        selected = counts >= min_each
        self.obs_all += obs
        self.n_colours += selected
        self.obs_supported += np.where(selected, obs, 0.)
        self.positive_supported += selected & (obs > 0.)
        months = panel["months"]
        if month_preserving:
            self.month_movable_species += int(any(
                len(np.unique(y[months == m])) > 1
                for m in np.unique(months)
            ))
        rng = np.random.default_rng(seed)
        for rep in range(REPS):
            yp = permute_labels(y, months, rng, month_preserving)
            delta = same_colour_frequency(yp, i, j) - expected
            self.null_all[rep] += delta
            self.null_supported[rep] += np.where(selected, delta, 0.)

    def finish(self, name, month_preserving):
        if not self.n_species:
            return {"n_species": 0, "case": name,
                    "warning": "insufficient source opportunity; no p-values"}
        n = self.n_species
        obs_all = self.obs_all / n
        null_all = self.null_all / n
        mean_delta = float(obs_all.sum())
        mean_null = float(null_all.sum(axis=1).mean())
        p_total = float((1 + (null_all.sum(axis=1) >= mean_delta).sum()) / (REPS + 1))
        colour_obs = np.divide(self.obs_supported, self.n_colours,
                               out=np.full(4, np.nan), where=self.n_colours > 0)
        colour_null = np.divide(self.null_supported, self.n_colours,
                                out=np.full_like(self.null_supported, np.nan),
                                where=self.n_colours[None, :] > 0)
        n_colour = self.n_colours
        ok = n_colour >= 20
        means = np.nanmean(np.where(ok, colour_null, np.nan), axis=0)
        sd = np.nanstd(colour_null, axis=0, ddof=1)
        sd = np.where(np.isfinite(sd) & (sd > 0), sd, np.nan)
        maxima = np.nanmax(np.where(ok, (colour_null - means) / sd, np.nan), axis=1) \
            if np.any(ok) else np.full(REPS, np.nan)
        per_colour = {}
        for k, c in enumerate(COLORS):
            if not ok[k]:
                per_colour[c] = {
                    "supported_species": int(n_colour[k]),
                    "status": "not_estimable_less_than_20_species"
                }
                continue
            z = float((colour_obs[k] - means[k]) / sd[k])
            per_colour[c] = {
                "supported_species": int(n_colour[k]),
                "observed_local_same_colour_excess": float(colour_obs[k]),
                "matched_null_mean": float(means[k]),
                "matched_null_q025_q975": list(map(float,
                    np.quantile(colour_null[:, k], [.025, .975]))),
                "observed_minus_matched_null": float(colour_obs[k] - means[k]),
                "matched_null_positive_tail_p_199": float(
                    (1 + np.count_nonzero(colour_null[:, k] >= colour_obs[k])) / (REPS + 1)
                ),
                "family4_maxT_upper_p_199": float(
                    (1 + np.count_nonzero(maxima >= z)) / (REPS + 1)
                ),
                "fraction_species_local_same_colour_exceeds_global": float(
                    self.positive_supported[k] / n_colour[k]
                ),
            }
        return {
            "case": name,
            "n_species": n,
            "n_original_minority_supported_richness": self.richness,
            "n_month_permutable_species": self.month_movable_species
               if month_preserving else None,
            "overall_observed_mean_global_minus_local_mismatch": mean_delta,
            "overall_matched_null_mean": mean_null,
            "overall_upper_p_199": p_total,
            "overall_adds_exactly_from_four_colour_components": bool(
                np.isclose(float(obs_all.sum()), mean_delta)
            ),
            "colour_components_all_included_species_means": {
                c: float(obs_all[k]) for k, c in enumerate(COLORS)
            },
            "colour_specific_supported_purity_not_genotype": per_colour,
        }


def analyze(photo):
    accum = {name: Accumulator() for name, *_ in CASES}
    counter = 0
    for taxon, g in photo.groupby("inat_taxon_id", sort=True):
        tid = int(taxon)
        for case_id, (name, support, observer, dominant, temporal) in enumerate(CASES):
            panel = prepare(g, support, observer, dominant)
            if panel is None:
                continue
            accum[name].add(panel, SEED + 101 * tid + 104729 * case_id,
                            support, temporal)
        counter += 1
        if counter % 1000 == 0:
            print("FCP_COLOUR_COMPONENT_CLUSTER_SPECIES_PROCESSED", counter, flush=True)
    results = {
        name: accum[name].finish(name, temporal)
        for name, _, _, _, temporal in CASES
    }
    for case, expected in (
        ("two_colors_5photos", 813),
        ("one_photo_per_observer_5photos", 727),
        ("two_colors_10photos", 462),
        ("one_photo_per_observer_10photos", 406),
        ("two_colors_5photos_month_preserving_null", 813),
    ):
        assert_truth(results[case]["n_species"] == expected,
                     "old authenticated source support tier mismatch for " + case)
    return {
        "schema": "fcp_authentic_447500_photo_colour_specific_50km_cluster_null_v1",
        "original_photo_sha256": PHOTO_SHA,
        "source_ledger_sha256": LEDGER_SHA,
        "n_original_photos": 447500,
        "n_classified_photos": 227107,
        "n_original_taxa": 4475,
        "spatial_radius_km": RADIUS_KM,
        "permutation_worlds": REPS,
        "source_photo_ids_and_exact_coordinates_exported": False,
        "familywise_hue_test": "four-colour one-sided maxT, studentised across shared permutations",
        "nulls": ["species-label/composition with fixed geo vertices",
                  "month-of-year preserving species-label/composition; exploratory"],
        "source_sampling": "every taxon had 100 photos; original max 2 photos per observer",
        "case_rules": {
           "base": ">=40 retained classifiable photo records, >=30 50-km local pairs, at least two hues with 5 or 10 photos",
           "observer": "single photo per nonmissing observer selected by observer/photo ID before colour-label interpretation",
           "dominant": ">=65% of sampled ROI palette belongs to labelled group and >=20 percentage point margin over next group; only compositional QC, not morph certainty",
           "month": "source photo labels shuffled within species AND calendar month; month-conditioned P does not establish same-year or same-breeding-population effect"
        },
        "results": results,
        "biological_limitations": [
            "Source labels are coarse photographed colour categories, not genotype-verified morphs.",
            "Colours photographed in one corolla or at different anthesis stages can make false biological clusters.",
            "Photographic geographic sorting could be caused by systematic regional camera/observer/petal-organ selection.",
            "Species equal weighting does not correct species selection, trait missingness or native/cultivated confusion.",
            "The tests are post hoc and use one source provider; min attainable P=0.005.",
            "No data identify maintenance costs, fitness, neutral genetic structure or local interbreeding populations.",
        ],
    }


def self_test():
    d = pd.DataFrame({
        "morph": ["white", "white", "blue_purple", "blue_purple"],
        "latitude": [0., 0., 20., 20.],
        "longitude": [0., .001, 0., .001],
        "positional_accuracy_m": [200.] * 4,
        "observer_id": [1, 2, 3, 4],
        "photo_id": [11, 12, 13, 14],
        "month_of_year": [1, 1, 2, 2],
        "primary_colour_pixel_fraction": [.9] * 4,
        "other_max_colour_pixel_fraction": [.05] * 4,
    })
    y = pd.Categorical(d.morph, categories=COLORS).codes.astype(np.int8)
    a, b = np.triu_indices(4, 1)
    close = haversine(d.latitude.to_numpy(), d.longitude.to_numpy(), a, b) <= 50
    i, j = a[close], b[close]
    stat = same_colour_frequency(y, i, j)
    assert_truth(len(i) == 2 and np.allclose(stat[[0, 3]], [.5, .5]),
                 "synthetic clustering wrong")
    counts = np.bincount(y, minlength=4)
    expected = counts * (counts - 1) / (4 * 3)
    assert_truth(np.isclose((stat - expected).sum(), 2 / 3),
                 "within species colour decomposition must sum to global-local mismatch")
    rng = np.random.default_rng(1)
    z = permute_labels(y, np.array([1, 1, 2, 2]), rng, True)
    assert_truth(np.array_equal(z, y), "month-specific class counts were not preserved")
    zz = permute_labels(y, np.array([1, 1, 2, 2]), rng, False)
    assert_truth(sorted(zz) == sorted(y), "global permutation alters class counts")
    assert_truth(len(choose_observers(d)) == 4, "observer thinning drops unrelated photos")
    print("FCP_COLOUR_SPECIFIC_SPATIAL_CLUSTER_SELF_TEST_PASS", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--photos", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    self_test()
    if args.self_test:
        return
    assert_truth(args.photos and args.out, "photos and out are required")
    ledger = Path("results/fcp_np_all4475_unified_species_ledger_20261021.csv")
    photo = read_source(args.photos, ledger)
    summary = analyze(photo)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "colour_specific_spatial_clusters_20261023.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True,
                   allow_nan=False) + "\n", encoding="utf-8"
    )
    print("FCP_COLOUR_SPECIFIC_SPATIAL_CLUSTER_SOURCE_VERIFIED", json.dumps({
        "source_sha256": summary["original_photo_sha256"],
        "n_source_photos": summary["n_original_photos"],
        "n_classified": summary["n_classified_photos"],
        "cases": {k: {
            "n_species": v["n_species"],
            "observed_mean": v.get("overall_observed_mean_global_minus_local_mismatch"),
            "upper_p_199": v.get("overall_upper_p_199"),
            "per_colour": {c: {
                "supported_species": w.get("supported_species"),
                "effect": w.get("observed_minus_matched_null"),
                "family4_p": w.get("family4_maxT_upper_p_199")
            } for c,w in v.get("colour_specific_supported_purity_not_genotype", {}).items()}
        } for k,v in summary["results"].items()}
    }, sort_keys=True, allow_nan=False), flush=True)


if __name__ == "__main__":
    main()
