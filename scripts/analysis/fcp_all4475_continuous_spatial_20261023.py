#!/usr/bin/env python3
"""Nine-palette continuous-colour geography in the authenticated FCP photo atlas.

Post hoc analyses of photographs, not genetic morph or population inference.
Reuses the source-locked normalization and base-2 Jensen-Shannon divergence.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "analysis"))
from fcp_full4475_hypothesis_census_20261021 import load_original, SOURCE_SHA, PHOTO_LABELS
from run_polymorphism_distributed_polymorphism_robustness_20261007 import pairwise_jsd

RADII = np.array([5., 25., 50., 100., 250.])
REPS = 199
EARTH = 6371.0088


def require(ok, msg):
    if not ok:
        raise RuntimeError("FCP_CONTINUOUS_SPATIAL: " + msg)


def one_taxon(group, tid, one_observer=False, reps=REPS):
    original_n = len(group)
    if original_n < 40:
        return None
    if one_observer:
        group = (group.dropna(subset=["observer_id"])
                 .sort_values(["observer_id", "photo_id"], kind="stable")
                 .drop_duplicates("observer_id").sort_index(kind="stable"))
        if len(group) < 40:
            return None
    valid = (group.latitude.between(-90, 90) & group.longitude.between(-180, 180)
             & group.positional_accuracy_m.between(0, 5000))
    group = group.loc[valid].reset_index(drop=True)
    n = len(group)
    if n < 2:
        return None
    a, b = np.triu_indices(n, k=1)
    lat = np.deg2rad(group.latitude.to_numpy(float))
    lon = np.deg2rad(group.longitude.to_numpy(float))
    z = (np.sin((lat[a] - lat[b]) / 2)**2
         + np.cos(lat[a])*np.cos(lat[b])*np.sin((lon[a] - lon[b]) / 2)**2)
    geo = 2*EARTH*np.arcsin(np.sqrt(np.clip(z, 0, 1)))
    near = np.flatnonzero(geo <= RADII[-1])
    order = near[np.argsort(geo[near], kind="stable")]
    counts = np.searchsorted(geo[order], RADII, side="right")
    if counts[2] < 30:
        return None

    palette = np.stack(group.palette_norm9.to_numpy())
    dis = pairwise_jsd(palette)
    full = float(dis[a, b].mean())
    u, v = a[order], b[order]
    obs_cum = np.cumsum(dis[u, v])
    observed = np.array([
        full - obs_cum[c-1]/c if c >= 30 else np.nan for c in counts
    ])
    months = group.date.dt.month.fillna(0).astype(int).to_numpy()
    strata = [idx for m in np.unique(months)
              if len(idx := np.flatnonzero(months == m)) > 1]
    rng = np.random.default_rng(20261023 + 31*int(tid) + int(one_observer))
    null = np.full((2, reps, len(RADII)), np.nan)
    for iteration in range(reps):
        shuffled = rng.permutation(n)
        month_shuffled = np.arange(n)
        for idx in strata:
            month_shuffled[idx] = rng.permutation(idx)
        for family, perm in enumerate((shuffled, month_shuffled)):
            cum = np.cumsum(dis[perm[u], perm[v]])
            for k, count in enumerate(counts):
                if count >= 30:
                    null[family, iteration, k] = full - cum[count-1]/count

    labels = pd.Categorical(group.morph, categories=PHOTO_LABELS).codes
    require(np.all(labels >= 0), "unexpected photographic label")
    coarse_50 = float(
        np.mean(labels[a] != labels[b])
        - np.mean(labels[u[:counts[2]]] != labels[v[:counts[2]]])
    )
    frequencies = group.morph.value_counts().reindex(PHOTO_LABELS, fill_value=0).to_numpy()
    row = {
        "taxon_id": int(tid), "species": str(group.species.iloc[0]),
        "n_classified_original": original_n, "n_valid_locations": n,
        "k_observed_valid": int((frequencies > 0).sum()),
        "k_five_valid": int((frequencies >= 5).sum()),
        "coarse_50km_delta": coarse_50,
        "n_months_swappable_colour": sum(
            np.unique(palette[idx], axis=0).shape[0] > 1 for idx in strata
        ),
    }
    for k, rad in enumerate(RADII):
        row[f"pairs_{int(rad)}km"] = int(counts[k])
        if np.isfinite(observed[k]):
            row[f"raw_jsd_depletion_{int(rad)}km"] = float(observed[k])
            for family, name in enumerate(("unrestricted", "month")):
                row[f"excess_jsd_{int(rad)}km_{name}"] = float(
                    observed[k] - null[family, :, k].mean()
                )
        else:
            row[f"raw_jsd_depletion_{int(rad)}km"] = None
            row[f"excess_jsd_{int(rad)}km_unrestricted"] = None
            row[f"excess_jsd_{int(rad)}km_month"] = None
    return row, null


def holm(p):
    order = np.argsort(p)
    result = np.zeros(len(p))
    max_p = 0.
    for rank, idx in enumerate(order):
        max_p = max(max_p, min(1., (len(p)-rank)*p[idx]))
        result[idx] = max_p
    return result.tolist()


def summarise(outcomes, reps):
    output = {}
    for family, name in enumerate(("unrestricted", "month")):
        radii = []
        for k, rad in enumerate(RADII):
            selected = [
                (r, arr[family, :, k]) for r, arr in outcomes
                if r[f"raw_jsd_depletion_{int(rad)}km"] is not None
            ]
            if not selected:
                radii.append({"radius_km": float(rad), "n_species": 0, "p_upper": 1.})
                continue
            obs = np.array([r[f"raw_jsd_depletion_{int(rad)}km"] for r, _ in selected])
            worlds = np.vstack([v for _, v in selected]).mean(axis=0)
            radii.append({
                "radius_km": float(rad),
                "n_species": len(selected),
                "mean_raw_depletion": float(obs.mean()),
                "mean_null_depletion": float(worlds.mean()),
                "mean_excess": float(obs.mean() - worlds.mean()),
                "p_upper": float((1 + np.count_nonzero(worlds >= obs.mean())) / (reps+1)),
                "n_species_swappable_months": sum(
                    r["n_months_swappable_colour"] > 0 for r, _ in selected
                ),
            })
        for row, adj in zip(radii, holm([r["p_upper"] for r in radii])):
            row["holm_p_five_radii"] = adj
        output[name] = radii

    descriptive = {}
    for key in ("k_observed_valid", "k_five_valid"):
        descriptive[key] = {}
        for k in (1, 2, 3, 4):
            selected = [r for r, _ in outcomes if r[key] == k]
            descriptive[key][str(k)] = {
                "n_species": len(selected),
                "mean_month_adjusted_50km_excess": (
                    float(np.mean([r["excess_jsd_50km_month"] for r in selected]))
                    if selected else None
                ),
            }
    return {
        "n_species_50km": len(outcomes),
        "tests": output,
        "richness_descriptive_not_morph_count": descriptive,
    }


def analyze(frame, one_observer=False, reps=REPS):
    outcomes = []
    for i, (tid, group) in enumerate(frame.groupby("inat_taxon_id", sort=True)):
        hit = one_taxon(group, tid, one_observer=one_observer, reps=reps)
        if hit is not None:
            outcomes.append(hit)
        if (i + 1) % 1000 == 0:
            print("FCP_CONTINUOUS_TAXA", one_observer, i+1, flush=True)
    require(bool(outcomes), "no estimable species")
    if not one_observer:
        require(len(outcomes) == 1761, "frozen strict 50-km species gate drift")
        coarse = np.mean([r["coarse_50km_delta"] for r, _ in outcomes])
        require(abs(coarse - .0176149021685411) < 1e-10,
                "frozen coarse four-state 50-km statistic not reproduced")
    return summarise(outcomes, reps), pd.DataFrame([r for r, _ in outcomes])


def self_test():
    rows = []
    for colour, latitude in ((0, 35.), (1, 39.)):
        for j in range(24):
            p = np.zeros(9)
            p[0 if colour == 0 else 6] = 1
            rows.append({
                "inat_taxon_id": 77, "species": "Example plant",
                "observer_id": 200*colour+j, "photo_id": 200*colour+j,
                "morph": PHOTO_LABELS[colour], "date": pd.Timestamp("2024-06-01"),
                "latitude": latitude, "longitude": 135.+.001*j,
                "positional_accuracy_m": 50., "palette_norm9": p,
            })
    data = pd.DataFrame(rows)
    a, _ = one_taxon(data, 77, reps=49)
    require(a["excess_jsd_50km_month"] > .2,
            "same-month geographic colour patches were missed")
    data.loc[data.morph == PHOTO_LABELS[0], "date"] = pd.Timestamp("2024-05-01")
    b, null = one_taxon(data, 77, reps=49)
    require(abs(b["excess_jsd_50km_month"]) < 1e-12
            and np.ptp(null[1, :, 2]) < 1e-12,
            "month-separated phenotype improperly declared spatially structured")
    data["date"] = pd.Timestamp("2024-06-01")
    data["morph"] = PHOTO_LABELS[0]
    c, _ = one_taxon(data, 77, reps=49)
    require(c["k_observed_valid"] == 1 and c["excess_jsd_50km_month"] > .2,
            "continuous geographic structure within a coarse state missed")
    print("FCP_ALL4475_CONTINUOUS_SPATIAL_SYNTHETIC_PASS")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--photos", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    self_test()
    if args.self_test:
        return
    require(args.photos and args.out, "photos and output directory required")
    all_source, classified = load_original(args.photos)
    results = {}
    args.out.mkdir(parents=True, exist_ok=True)
    for thinned in (False, True):
        name = "one_per_observer" if thinned else "original"
        summary, taxa = analyze(classified, one_observer=thinned)
        taxa.to_csv(args.out / f"{name}_per_taxon.csv", index=False)
        results[name] = summary
    record = {
        "schema": "fcp_all4475_multiscale_continuous_nine_palette_v1",
        "source_sha256": SOURCE_SHA,
        "source_original_photos": len(all_source),
        "source_four_state_classifiable": len(classified),
        "radii_km": RADII.tolist(),
        "permutations": REPS,
        "distance": "base-2 Jensen-Shannon divergence, not its square root",
        "null": "permute whole nine-hue vectors within species; in month null preserve species-calendar-month composition",
        "same_photo_provider": True,
        "posthoc": True,
        "not_inherited_morphs_populations_or_selection": True,
        "results": results,
    }
    (args.out / "summary.json").write_text(
        json.dumps(record, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print("FCP_ALL4475_CONTINUOUS_SPATIAL_FINISHED",
          json.dumps({k: v["tests"]["month"][2] for k, v in results.items()}))


if __name__ == "__main__":
    main()
