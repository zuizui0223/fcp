#!/usr/bin/env python3
"""Source-locked nine-palette flower-colour geography on the FCP 4,475-taxon frame.

Uses ORIGINAL photographed/ROI-assigned hue fractions, not invented CIELAB
or newly estimated genetic morphs. Preserves earlier P2 outputs untouched.
Three within-taxon complete-vector nulls assess geographic similarity at
5, 25, 50, 100 and 250 km:
  1. any photo vertex (species palette composition and positions fixed),
  2. only within calendar month (seasonal hue composition fixed), and
  3. only within calendar-month x frozen four-state photo label (a residual
     continuous hue pattern, above-and-beyond the original coarse classes).
Permutations, colors, and coordinates are never exported at the photo level.
All newly reported tests are exploratory post hoc, not fitness or genetics.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

PHOTO_SHA256 = "8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6"
PALETTE = ("white", "yellow", "orange", "red", "pink", "magenta", "purple", "blue", "bronze")
COARSE = ("white", "yellow_orange", "red_pink", "blue_purple")
FRACTIONS = tuple(f"flower_fraction_{c}" for c in PALETTE)
COLUMNS = ("inat_taxon_id", "species", "photo_id", "observer_id", "observed_on",
           "global_classifiable", "morph", "latitude", "longitude",
           "positional_accuracy_m", *FRACTIONS)
RADII_KM = (5., 25., 50., 100., 250.)
NULLS = ("full_species_shuffle", "within_month", "within_month_and_coarse_label")
PERMUTATIONS = 199
SEED = 20261023
EARTH_RADIUS_KM = 6371.0088
EXPECTED_ROWS = 447500
EXPECTED_CLASSIFIED = 227107
EXPECTED_TAXA = 4475
P2_ANY_PAIRS_N = 3600
P2_ANY_PAIRS_DELTA = 0.016115883798023662
P2_STRICT_N = 1761
P2_STRICT_DELTA = 0.0176149021685411


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(f"FCP_9PALETTE_GEOGRAPHY: {message}")


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_source(path: Path, authenticated: bool = True) -> pd.DataFrame:
    if authenticated:
        require(file_sha256(path) == PHOTO_SHA256, "immutable source photo bytes changed")
    original = pd.read_csv(path, usecols=list(COLUMNS), low_memory=False)
    if authenticated:
        require(len(original) == EXPECTED_ROWS, "original photo number changed")
        require(original.photo_id.nunique() == EXPECTED_ROWS, "source photo IDs not unique")
        require(original.inat_taxon_id.nunique() == EXPECTED_TAXA, "taxon frame changed")
        require(original.groupby("inat_taxon_id").size().eq(100).all(), "source sample not exactly 100/taxon")
    flag = original.global_classifiable.fillna(False).astype(str).str.strip().str.lower().isin(("true", "1", "yes"))
    selected = original.loc[flag].copy()
    if authenticated:
        require(len(selected) == EXPECTED_CLASSIFIED, "source classification denominator changed")
        require(selected.inat_taxon_id.nunique() == EXPECTED_TAXA, "taxon lost at classification")
    require(selected.morph.isin(COARSE).all(), "unexpected coarse source labels")
    for key in ("latitude", "longitude", "positional_accuracy_m", *FRACTIONS):
        selected[key] = pd.to_numeric(selected[key], errors="coerce")
    palette = selected[list(FRACTIONS)].to_numpy(float)
    require(np.isfinite(palette).all(), "nonfinite source palette components")
    require((palette >= 0).all() and (palette.sum(axis=1) > 0).all(), "invalid source nine-palette proportions")
    selected["calendar_month"] = pd.to_datetime(selected.observed_on, errors="coerce").dt.month.fillna(0).astype("int8")
    return selected


def jsd_base2_matrix(vectors: np.ndarray) -> np.ndarray:
    """Jensen-Shannon *divergence* base 2 in [0,1], no sqrt."""
    x = np.asarray(vectors, dtype=float)
    require(x.ndim == 2 and x.shape[1] == 9, "expected exactly nine hue channels")
    require(np.isfinite(x).all() and (x >= 0).all(), "invalid continuous hues")
    total = x.sum(axis=1)
    require((total > 0).all(), "zero-mass hue row")
    x = x / total[:, None]
    a = x[:, None, :]
    b = x[None, :, :]
    m = 0.5 * (a + b)
    with np.errstate(divide="ignore", invalid="ignore"):
        left = np.where(a > 0, a * np.log2(a / m), 0).sum(axis=2)
        right = np.where(b > 0, b * np.log2(b / m), 0).sum(axis=2)
    z = np.clip(0.5 * (left + right), 0, 1)
    return (z + z.T) * 0.5


def distances_km(lat: np.ndarray, lon: np.ndarray, a: np.ndarray, b: np.ndarray) -> np.ndarray:
    lat = np.deg2rad(np.asarray(lat, dtype=float))
    lon = np.deg2rad(np.asarray(lon, dtype=float))
    u = np.sin((lat[a] - lat[b]) / 2) ** 2
    u += np.cos(lat[a]) * np.cos(lat[b]) * np.sin((lon[a] - lon[b]) / 2) ** 2
    return 2.0 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(np.clip(u, 0, 1)))


def stratified_groups(months: np.ndarray, labels: np.ndarray) -> tuple[list[np.ndarray], list[np.ndarray]]:
    only_month = []
    month_label = []
    for month in np.unique(months):
        ix = np.flatnonzero(months == month)
        if len(ix) > 1:
            only_month.append(ix)
    tokens = months.astype("int16") * 4 + labels.astype("int16")
    for token in np.unique(tokens):
        ix = np.flatnonzero(tokens == token)
        if len(ix) > 1:
            month_label.append(ix)
    return only_month, month_label


def taxon_result(frame: pd.DataFrame, taxon_id: int, *, one_per_observer: bool = False,
                 permutations: int = PERMUTATIONS):
    """An eligible taxon row, permutation matrix; no photo IDs/coordinates leave this function."""
    source_n = len(frame)
    frame = frame.sort_values("photo_id", kind="stable")
    if one_per_observer:
        frame = (frame.dropna(subset=["observer_id"])
                 .sort_values(["observer_id", "photo_id"], kind="stable")
                 .drop_duplicates("observer_id").sort_index(kind="stable"))
    after_thinning_n = len(frame)
    geo = (frame.latitude.between(-90, 90) & frame.longitude.between(-180, 180)
           & frame.positional_accuracy_m.between(0, 5000))
    g = frame.loc[geo].reset_index(drop=True)
    n = len(g)
    if n < 2:
        return None
    a, b = np.triu_indices(n, k=1)
    km = distances_km(g.latitude.to_numpy(float), g.longitude.to_numpy(float), a, b)
    candidates = np.flatnonzero(km <= RADII_KM[-1])
    if len(candidates) == 0:
        return None
    ordered = candidates[np.argsort(km[candidates], kind="stable")]
    counts = np.searchsorted(km[ordered], RADII_KM, side="right")
    n50 = int(counts[2])
    if n50 == 0:
        return None
    ui, uj = a[ordered], b[ordered]
    whole = jsd_base2_matrix(g[list(FRACTIONS)].to_numpy(float))
    global_jsd = float(np.mean(whole[a, b]))
    observed_cs = np.cumsum(whole[ui, uj])
    raw = np.full(len(RADII_KM), np.nan)
    for j, count in enumerate(counts):
        if count:
            raw[j] = global_jsd - float(observed_cs[count - 1] / count)
    y = pd.Categorical(g.morph, categories=COARSE).codes.astype(np.int8)
    require(np.all(y >= 0), "invalid source coarse category")
    coarse_global = float(np.mean(y[a] != y[b]))
    coarse_delta = coarse_global - float(np.mean(y[ui[:n50]] != y[uj[:n50]]))
    months = g.calendar_month.to_numpy("int8")
    month_groups, month_label_groups = stratified_groups(months, y)
    rand = np.random.default_rng(SEED + int(taxon_id) * 31 + int(one_per_observer))
    null = np.full((len(NULLS), permutations, len(RADII_KM)), np.nan)
    for iteration in range(permutations):
        unrestricted = rand.permutation(n)
        monthly = np.arange(n)
        for ix in month_groups:
            monthly[ix] = rand.permutation(ix)
        within_label_month = np.arange(n)
        for ix in month_label_groups:
            within_label_month[ix] = rand.permutation(ix)
        for f, permutation in enumerate((unrestricted, monthly, within_label_month)):
            shuffled_cs = np.cumsum(whole[permutation[ui], permutation[uj]])
            for j, count in enumerate(counts):
                if count:
                    null[f, iteration, j] = global_jsd - float(shuffled_cs[count - 1] / count)
    counts_by_label = np.bincount(y, minlength=4)
    row = {
        "inat_taxon_id": int(taxon_id), "species": str(g.species.iloc[0]),
        "genus": str(g.species.iloc[0]).split(" ")[0],
        "n_source_classified": int(source_n), "n_after_thinning": int(after_thinning_n),
        "n_geo_valid": int(n),
        "coarse_label_richness_observed": int(np.count_nonzero(counts_by_label)),
        "coarse_label_richness_5photos": int(np.count_nonzero(counts_by_label >= 5)),
        "global_mean_jsd": global_jsd, "global_coarse_mismatch": coarse_global,
        "coarse_delta_50km": coarse_delta,
        "n_month_strata_with_diverse_coarse_labels": sum(len(np.unique(y[ix])) > 1 for ix in month_groups),
        "n_coarse_month_strata_with_palette_variation": sum(
            len(np.unique(g[list(FRACTIONS)].to_numpy(float)[ix], axis=0)) > 1
            for ix in month_label_groups),
    }
    for j, radius in enumerate(RADII_KM):
        key = f"{int(radius)}km"
        row[f"n_pairs_{key}"] = int(counts[j])
        row[f"raw_depletion_{key}"] = float(raw[j]) if np.isfinite(raw[j]) else None
        for f, label in enumerate(NULLS):
            row[f"excess_{key}_{label}"] = (float(raw[j] - null[f, :, j].mean())
                                                 if np.isfinite(raw[j]) else None)
    return row, null


def genus_bootstrap_ci(rows: list[dict], value_key: str, reps: int = 499, seed: int = SEED):
    groups = {}
    for row in rows:
        groups.setdefault(row["genus"], []).append(row[value_key])
    if len(groups) < 2:
        return [None, None]
    count = np.array([len(v) for v in groups.values()], dtype=int)
    total = np.array([np.sum(v) for v in groups.values()], dtype=float)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(groups), size=(reps, len(groups)))
    values = total[idx].sum(axis=1) / count[idx].sum(axis=1)
    return [float(v) for v in np.quantile(values, [.025, .975])]


def holm_adjust(pvalues: list[float]) -> list[float]:
    a = np.asarray(pvalues, float)
    order = np.argsort(a, kind="stable")
    adjusted = np.zeros(len(a), float)
    floor = 0.
    for place, index in enumerate(order):
        floor = max(floor, min(1., (len(a) - place) * a[index]))
        adjusted[index] = floor
    return adjusted.tolist()


def make_summary(results: list[tuple[dict, np.ndarray]], permutations: int, *,
                 strict: bool, min_pairs: int = 30, do_bootstrap: bool = True,
                 one_per_observer: bool = False) -> dict:
    """One test per radius/null; strict cohort must pass source 40/30 at 50km."""
    eligible = []
    for row, null in results:
        if strict and not (row["n_source_classified"] >= 40 and row["n_pairs_50km"] >= 30
                           and (not one_per_observer or row["n_after_thinning"] >= 40)):
            continue
        eligible.append((row, null))
    summary = {"n_species_with_50km_pair": len(eligible), "statistics": {},
               "thresholds": {"min_species_classified": 40 if strict else 2,
                              "min_local_pairs_per_radius": min_pairs}}
    across_family_p = []
    across_family_records = []
    for f, family in enumerate(NULLS):
        fam = []
        for j, r in enumerate(RADII_KM):
            name = f"{int(r)}km"
            supported = [(x, n[f, :, j]) for x, n in eligible if x[f"n_pairs_{name}"] >= min_pairs]
            if not supported:
                fam.append({"radius_km": float(r), "n_species": 0,
                            "mean_raw_depletion": None, "mean_excess": None,
                            "p_upper_199": None, "genus_bootstrap_excess_95ci": None})
                continue
            vals = np.array([x[f"raw_depletion_{name}"] for x, _ in supported])
            means = np.stack([n for _, n in supported]).mean(axis=0)
            observed = float(vals.mean())
            xs = [x for x, _ in supported]
            excess_key = f"excess_{name}_{family}"
            excess = float(np.mean([x[excess_key] for x in xs]))
            p_upper = float((1 + int(np.count_nonzero(means >= observed))) / (permutations + 1))
            ci = (genus_bootstrap_ci(xs, excess_key, seed=SEED+101*f+17*j)
                  if do_bootstrap else None)
            entry = {"radius_km": float(r), "n_species": len(xs),
                     "n_genera": len(set(x["genus"] for x in xs)),
                     "mean_raw_depletion": observed,
                     "mean_null_depletion": float(means.mean()), "mean_excess": excess,
                     "null_q025_q975": [float(v) for v in np.quantile(means, [.025, .975])],
                     "genus_bootstrap_excess_95ci": ci,
                     "fraction_positive_raw": float(np.mean(vals > 0)),
                     "p_upper_199": p_upper}
            fam.append(entry)
            across_family_p.append(p_upper)
            across_family_records.append(entry)
        valid_p = [r["p_upper_199"] for r in fam if r["p_upper_199"] is not None]
        adjusted = iter(holm_adjust(valid_p))
        for r in fam:
            r["holm_within_five_radii"] = (next(adjusted) if r["p_upper_199"] is not None else None)
        summary["statistics"][family] = fam
    for r, p in zip(across_family_records, holm_adjust(across_family_p)):
        r["holm_across_all_nulls_and_radii"] = p
    return summary


def describe_richness(results: list[tuple[dict, np.ndarray]]) -> dict:
    eligible = [(x, n) for x, n in results if x["n_source_classified"] >= 40 and x["n_pairs_50km"] >= 30]
    output = {}
    for field in ("coarse_label_richness_observed", "coarse_label_richness_5photos"):
        output[field] = {}
        for k in (1, 2, 3, 4):
            subset = [x for x, _ in eligible if x[field] == k]
            output[field][str(k)] = {
                "n_species": len(subset),
                "mean_month_excess_50km": (float(np.mean([r["excess_50km_within_month"] for r in subset])) if subset else None),
                "mean_coarse_month_residual_excess_50km": (
                    float(np.mean([r["excess_50km_within_month_and_coarse_label"] for r in subset])) if subset else None),
            }
    return output


def analyze(frame: pd.DataFrame, *, one_per_observer: bool, permutations: int):
    results = []
    n = frame.inat_taxon_id.nunique()
    for at, (taxon, g) in enumerate(frame.groupby("inat_taxon_id", sort=True), 1):
        z = taxon_result(g, int(taxon), one_per_observer=one_per_observer, permutations=permutations)
        if z is not None:
            results.append(z)
        if at % 1000 == 0:
            print("FCP_9PALETTE_TAXA", "thinned" if one_per_observer else "original", at, flush=True)
    require(n == EXPECTED_TAXA, "analysis taxa differ from all 4475 source roster")
    if not one_per_observer:
        require(len(results) == P2_ANY_PAIRS_N, "old 3,600-species 50km denominator differs")
        base = float(np.mean([r["coarse_delta_50km"] for r, _ in results]))
        require(abs(base-P2_ANY_PAIRS_DELTA)<1e-10, "old 3,600-species delta differs")
        strict = [r for r, _ in results if r["n_source_classified"]>=40 and r["n_pairs_50km"]>=30]
        require(len(strict)==P2_STRICT_N, "old 1,761-species strict denominator differs")
        strict_mean = float(np.mean([r["coarse_delta_50km"] for r in strict]))
        require(abs(strict_mean-P2_STRICT_DELTA)<1e-10, "old strict 50km delta differs")
    else:
        base = strict_mean = None
    return {
        "n_total_species_source": n,
        "n_with_50km_opportunity": len(results),
        "n_original_strict_gate": len([1 for r, _ in results if r["n_source_classified"]>=40 and r["n_pairs_50km"]>=30]),
        "old_original_3600_coarse_delta_reproduced": base,
        "old_original_1761_coarse_delta_reproduced": strict_mean,
        "strict_40photos_30pairs": make_summary(results, permutations, strict=True,
                                                  one_per_observer=one_per_observer),
        "broad_any_photo_pairs_50km": make_summary(results, permutations, strict=False, min_pairs=1,
                                                    do_bootstrap=False),
        "richness_descriptive_only": describe_richness(results),
    }, pd.DataFrame([r for r, _ in results])


def self_test():
    z = np.zeros((2, 9), float)
    z[0, 0] = z[1, 6] = 1
    require(abs(jsd_base2_matrix(z)[0, 1]-1) < 1e-12, "JSD unit scale mismatch")
    require(abs(jsd_base2_matrix(z)[0, 0]) < 1e-12, "JSD identity mismatch")
    rows = []
    for patch, latitude in ((0, 35.), (1, 39.)):
        for j in range(24):
            row = {"inat_taxon_id": 101, "species": "Test species", "photo_id": patch*100+j,
                   "observer_id": patch*100+j, "observed_on": "2024-06-01", "calendar_month": 6,
                   "morph": COARSE[patch], "latitude": latitude,
                   "longitude": 135.+.001*j, "positional_accuracy_m": 25.}
            row.update({f: 1.0 if i == (0 if patch == 0 else 6) else 0.0
                        for i, f in enumerate(FRACTIONS)})
            rows.append(row)
    s = pd.DataFrame(rows)
    result = taxon_result(s, 101, permutations=99)
    require(result is not None, "synthetic opportunity absent")
    info, null = result
    require(info["excess_50km_within_month"] > 0.2, "genuine month-matched spatial signal lost")
    require(abs(info["excess_50km_within_month_and_coarse_label"])<1e-12,
            "class-only geography wrongly counted as within-class residual")
    s.loc[s.morph==COARSE[0], "calendar_month"] = 5
    s.loc[s.morph==COARSE[1], "calendar_month"] = 6
    a, n = taxon_result(s, 101, permutations=49)
    require(abs(a["excess_50km_within_month"]) < 1e-12,
            "month-only pattern survived month conditional null")
    require(np.ptp(n[1, :, 2])<1e-12, "fixed-month null should be degenerate")
    s["calendar_month"] = 6
    s["morph"] = COARSE[0]
    b, _ = taxon_result(s, 101, permutations=49)
    require(b["coarse_label_richness_observed"] == 1, "synthetic source labels changed")
    require(b["excess_50km_within_month_and_coarse_label"]>0.2,
            "genuine within-label geographic colour variation missed")
    require(np.isclose(holm_adjust([0.01, 0.02, 0.03]), [0.03,0.04,0.04]).all(),
            "Holm correction broken")
    print("FCP_9PALETTE_GEOGRAPHY_SYNTHETIC_PASS", flush=True)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--photos", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--permutations", type=int, default=PERMUTATIONS)
    parser.add_argument("--self-test", action="store_true")
    z=parser.parse_args()
    self_test()
    if z.self_test:
        return
    require(z.photos is not None and z.out is not None, "photos and output required")
    require(z.permutations >= 19, "at least 19 permutation worlds required")
    source=read_source(z.photos)
    results={}
    z.out.mkdir(parents=True, exist_ok=True)
    for one_per_observer in (False, True):
        mode="one_per_observer" if one_per_observer else "original"
        summary, panel=analyze(source, one_per_observer=one_per_observer,
                               permutations=z.permutations)
        panel.to_csv(z.out / f"{mode}_taxon_metrics.csv", index=False)
        results[mode]=summary
    receipt={
        "schema": "fcp_all4475_original_nine_palette_multiscale_geography_posthoc_v1",
        "status": "measured_original_photo_palette_with_matched_nulls",
        "source_photo_sha256": PHOTO_SHA256,
        "source_classified_photos": EXPECTED_CLASSIFIED, "source_original_photos": EXPECTED_ROWS,
        "source_taxa": EXPECTED_TAXA, "radii_km": list(RADII_KM),
        "colour_channels": list(PALETTE), "original_four_states": list(COARSE),
        "distance": "base2 Jensen-Shannon divergence, not square-root JSD",
        "all_new_inference_posthoc": True,
        "source_photo_provider": "iNaturalist, one provider and one original ROI classifier",
        "null_families": list(NULLS), "n_permutations": z.permutations,
        "strict_gate": "n_source_classified>=40 and 50km pairs>=30; observer-thinned additionally n_after_thinning>=40; other radii >=30 pairs",
        "broad_gate": "n_classified_photo>=2 and n_50km_pairs>=1; other radii >=1 pair",
        "coordinate_gate": "source positional accuracy between 0 and 5000m",
        "limitations": ["Photo colour not inherited morphs", "No photo-level petal-error validation",
                        "Calendar month not within-year flowering stage control",
                        "Cultivated/non-native source not independently filtered",
                        "Source four-colour classification has 220393 unclassified photos",
                        "Original strict 5km same-year/month distinct-observer result remains unsupported",
                        "No biological spatial clusters or local genetic polymorphism inferred"],
        "results":results,
    }
    (z.out / "summary.json").write_text(json.dumps(receipt, indent=2, allow_nan=False)+"\n", encoding="utf8")
    print("FCP_9PALETTE_GEOGRAPHY_COMPLETE", json.dumps({
        mode:{"n_any50":r["n_with_50km_opportunity"],
              "strict_50month":r["strict_40photos_30pairs"]["statistics"]["within_month"][2],
              "strict_50_coarse_month":r["strict_40photos_30pairs"]["statistics"]["within_month_and_coarse_label"][2]}
        for mode,r in results.items()},allow_nan=False),flush=True)


if __name__ == "__main__":
    main()
