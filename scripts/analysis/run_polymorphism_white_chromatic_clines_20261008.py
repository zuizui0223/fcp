#!/usr/bin/env python3
"""Post-outcome screen: visible white/chromatic clines and local co-occurrence.

This analysis is intentionally separate from the frozen FCP manuscript claims.
Visible colours are photographic states, not measured anthocyanin or UV absorbance.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio

INPUT_SHA256 = {
    "discovery": "ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4",
    "validation": "0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6",
    "third": "57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186",
}
ELEVATION_ZIP_SHA256 = "6a9f4e9a37f289d594ed1567b63be344a4dacf247f0b69756e81df031087e3bd"
MORPHS = ("white", "yellow_orange", "red_pink", "blue_purple")
MIN_SPECIES_PHOTOS = 40
MIN_MORPH_PHOTOS = 5
MIN_LOCAL_PAIRS = 30
LOCAL_RADIUS_KM = 50.0
MIN_AXIS_SPAN = {"abs_lat": 1.0, "elevation_m": 100.0}
N_PERM = 499
N_BOOT = 1999
SEED = 2026100801
EARTH_RADIUS_KM = 6371.0088


def file_sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for piece in iter(lambda: f.read(1 << 20), b""):
            h.update(piece)
    return h.hexdigest()


def seed_for(*items: object) -> int:
    payload = "|".join(str(x) for x in (SEED, *items))
    return int.from_bytes(hashlib.sha256(payload.encode()).digest()[:8], "little")


def as_bool(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False)
    return values.fillna("").astype(str).str.lower().str.strip().isin({"true", "1", "yes", "y"})


def read_cohort(path: Path, cohort: str, *, check_hash: bool = True) -> pd.DataFrame:
    if check_hash and file_sha(path) != INPUT_SHA256[cohort]:
        raise ValueError(f"{cohort}: measured-source SHA256 mismatch")
    d = pd.read_csv(path, low_memory=False)
    needed = {"inat_taxon_id", "species", "photo_id", "latitude", "longitude",
              "morph", "global_classifiable"}
    if not needed.issubset(d):
        raise ValueError(f"{cohort}: missing columns {sorted(needed-set(d))}")
    d = d.loc[as_bool(d["global_classifiable"]) & d["morph"].isin(MORPHS),
              list(needed)].copy()
    for name in ("latitude", "longitude"):
        d[name] = pd.to_numeric(d[name], errors="coerce")
    d = d.loc[d.latitude.between(-90, 90) & d.longitude.between(-180, 180)].copy()
    n = d.groupby("inat_taxon_id").size()
    d = d.loc[d.inat_taxon_id.isin(n[n >= MIN_SPECIES_PHOTOS].index)].copy()
    if d["photo_id"].duplicated().any():
        raise ValueError(f"{cohort}: duplicated classifiable photo_id")
    d["abs_lat"] = d.latitude.abs()
    d["cohort"] = cohort
    d["genus"] = d["species"].astype(str).str.split().str[0]
    return d.sort_values(["inat_taxon_id", "photo_id"], kind="stable").reset_index(drop=True)


def sample_elevation(path: Path, d: pd.DataFrame) -> np.ndarray:
    pts = list(zip(d.longitude.to_numpy(float), d.latitude.to_numpy(float)))
    result = np.full(len(pts), np.nan)
    with rasterio.open(path) as raster:
        if raster.crs is None or raster.crs.to_epsg() != 4326:
            raise ValueError("Elevation raster must be geographic WGS84; no implicit reprojection")
        for start in range(0, len(pts), 2000):
            chunk = np.asarray(list(raster.sample(pts[start:start+2000])), dtype=float)[:, 0]
            if raster.nodata is not None:
                chunk[np.isclose(chunk, raster.nodata)] = np.nan
            result[start:start+len(chunk)] = chunk
    return result


def species_delta(x: np.ndarray, chromatic: np.ndarray, min_span: float,
                  *, rng: np.random.Generator, n_perm: int = N_PERM) -> tuple[float, np.ndarray] | None:
    """Equal-species standardized morph contrast; permute labels within species."""
    x = np.asarray(x, float)
    chromatic = np.asarray(chromatic, bool)
    good = np.isfinite(x)
    x, chromatic = x[good], chromatic[good]
    nw, nc = int((~chromatic).sum()), int(chromatic.sum())
    if nw < MIN_MORPH_PHOTOS or nc < MIN_MORPH_PHOTOS:
        return None
    if float(x.max() - x.min()) < min_span:
        return None
    sd = float(x.std(ddof=0))
    if sd <= 0:
        return None
    z = (x - float(x.mean())) / sd
    obs = float(z[chromatic].mean() - z[~chromatic].mean())
    null = np.empty(n_perm, float)
    # Exact within-species counts, locations and environmental distributions fixed.
    for b in range(n_perm):
        perm = rng.permutation(chromatic)
        null[b] = z[perm].mean() - z[~perm].mean()
    return obs, null


def residualize(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Within-species residual y after linear projection on x (sensitivity only)."""
    y = np.asarray(y, float)
    x = np.asarray(x, float)
    ans = np.full(len(y), np.nan)
    good = np.isfinite(y) & np.isfinite(x)
    if good.sum() < 2 or float(x[good].std()) <= 1e-10:
        return ans
    centered = x[good] - x[good].mean()
    beta = float(np.dot(centered, y[good]-y[good].mean()) / np.dot(centered, centered))
    ans[good] = y[good] - y[good].mean() - beta * centered
    return ans


def local_white_chromatic(g: pd.DataFrame) -> dict:
    """Observed local pair composition versus exact species-wide label opportunity."""
    n = len(g)
    if n < MIN_SPECIES_PHOTOS:
        return {"eligible": False}
    rad_lat = np.deg2rad(g.latitude.to_numpy(float))
    rad_lon = np.deg2rad(g.longitude.to_numpy(float))
    c = np.cos(rad_lat)
    xyz = np.column_stack([c*np.cos(rad_lon), c*np.sin(rad_lon), np.sin(rad_lat)])
    a, b = np.triu_indices(n, k=1)
    local = np.arccos(np.clip(np.sum(xyz[a]*xyz[b], axis=1), -1.0, 1.0))
    a, b = a[local*EARTH_RADIUS_KM <= LOCAL_RADIUS_KM], b[local*EARTH_RADIUS_KM <= LOCAL_RADIUS_KM]
    if len(a) < MIN_LOCAL_PAIRS:
        return {"eligible": False, "n_local_pairs": int(len(a))}
    code = pd.Categorical(g.morph, categories=MORPHS).codes
    counts = np.bincount(code, minlength=4)
    nw, nc = int(counts[0]), int(counts[1:].sum())
    p_wc = 2.0*nw*nc/(n*(n-1))
    p_hue = (nc*nc-int(np.sum(counts[1:]*counts[1:])))/(n*(n-1))
    white_a, white_b = code[a] == 0, code[b] == 0
    obs_wc = int(np.count_nonzero(white_a ^ white_b))
    obs_hue = int(np.count_nonzero(~white_a & ~white_b & (code[a] != code[b])))
    expected_wc, expected_hue = len(a)*p_wc, len(a)*p_hue
    return {
        "eligible": True, "n_local_pairs": int(len(a)),
        "white_chromatic_pairs": obs_wc,
        "white_chromatic_expected": float(expected_wc),
        "white_chromatic_ratio": float(obs_wc/expected_wc) if expected_wc >= 5 else None,
        "nonwhite_hue_pairs": obs_hue,
        "nonwhite_hue_expected": float(expected_hue),
        "nonwhite_hue_ratio": float(obs_hue/expected_hue) if expected_hue >= 5 else None,
        "has_local_white_chromatic_pair": bool(obs_wc > 0),
    }


def summarize_effect(entries: list[dict], predictor: str, cohort: str) -> dict:
    vals = np.array([r["delta"] for r in entries if r["predictor"] == predictor], float)
    nulls = [r["null"] for r in entries if r["predictor"] == predictor]
    if not len(vals):
        return {"estimable": False, "n_species": 0}
    mat = np.stack(nulls)
    null_means = mat.mean(axis=0)
    obs = float(vals.mean())
    p_upper = float((1+np.count_nonzero(null_means >= obs))/(N_PERM+1))
    p_two = float((1+np.count_nonzero(abs(null_means-null_means.mean()) >=
                                      abs(obs-null_means.mean())))/(N_PERM+1))
    rng = np.random.default_rng(seed_for("bootstrap", cohort, predictor))
    draws = vals[rng.integers(0, len(vals), (N_BOOT, len(vals)))].mean(axis=1)
    df = pd.DataFrame([(r["genus"], r["delta"]) for r in entries
                       if r["predictor"] == predictor], columns=["genus", "delta"])
    genus = df.groupby("genus").delta.mean()
    # Post hoc genus-level sign flip: conservatively change direction by genus,
    # rather than treating all photo locations/species as independent lineages.
    genus_values = genus.to_numpy(float)
    rng_genus = np.random.default_rng(seed_for("genus-signflip", cohort, predictor))
    signflips = rng_genus.choice((-1, 1), size=(9999, len(genus_values)))
    genus_null = (signflips @ genus_values)/len(genus_values)
    genus_obs = float(np.mean(genus_values))
    genus_p_positive = float((1+np.count_nonzero(genus_null>=genus_obs))/10000)
    return {
        "estimable": True, "n_species": int(len(vals)), "n_genera": int(len(genus)),
        "mean_delta_chromatic_minus_white_within_species_sd": obs,
        "median_species_delta_sd": float(np.median(vals)),
        "fraction_species_positive": float((vals > 0).mean()),
        "species_bootstrap_95CI_mean": [float(x) for x in np.quantile(draws, [0.025, 0.975])],
        "mean_genus_balanced_delta_sd": float(genus.mean()),
        "genus_signflip_p_positive": genus_p_positive,
        "genus_signflip_role": "posthoc genus-blocked sign-direction sensitivity, not full phylogenetic correction",
        "within_species_vertex_permutation_p_positive": p_upper,
        "within_species_vertex_permutation_p_two_sided": p_two,
        "null_mean": float(null_means.mean()),
    }


def holm_two(p1: float, p2: float) -> tuple[float, float]:
    # Primary simultaneous directional hypotheses: high |latitude| and high elevation.
    return (min(1.0, max(2*p1, p1) if p1 <= p2 else max(p1, 2*p2)),
            min(1.0, max(2*p2, p2) if p2 <= p1 else max(p2, 2*p1)))


def analyze(d: pd.DataFrame, cohort: str) -> tuple[dict, pd.DataFrame]:
    records = []
    effects = []
    for taxon_id, group in d.groupby("inat_taxon_id", sort=True):
        g = group.reset_index(drop=True)
        chrom = (g.morph.to_numpy() != "white")
        counts = g.morph.value_counts().reindex(MORPHS, fill_value=0)
        n_white, n_colour = int(counts["white"]), int(counts.iloc[1:].sum())
        coloured_hues = int((counts.iloc[1:] >= MIN_MORPH_PHOTOS).sum())
        row = {
            "inat_taxon_id": taxon_id, "species": str(g.species.iloc[0]),
            "genus": str(g.genus.iloc[0]), "cohort": cohort,
            "n_classifiable": len(g), "n_white": n_white, "n_colour": n_colour,
            "white_plus_colour": n_white >= MIN_MORPH_PHOTOS and n_colour >= MIN_MORPH_PHOTOS,
            "white_with_one_nonwhite_hue_at_least_5_each": (n_white >= MIN_MORPH_PHOTOS
                and coloured_hues >= 1),
            "dominant_top_two_white_nonwhite": (
                "white" in set(counts.sort_values(ascending=False, kind="stable").head(2).index)
                and int(np.sort(counts.to_numpy())[-2]) >= MIN_MORPH_PHOTOS
            ),
            "two_nonwhite_hues": coloured_hues >= 2,
            "abs_lat_span_degrees": float(g.abs_lat.max()-g.abs_lat.min()),
            "elevation_span_m": float(g.elevation_m.max()-g.elevation_m.min()) if g.elevation_m.notna().any() else None,
        }
        # The photo-derived "white" state is exposure-coupled; no functional pigment label implied.
        for predictor in ("abs_lat", "elevation_m"):
            x = g[predictor].to_numpy(float)
            fit = species_delta(x, chrom, MIN_AXIS_SPAN[predictor],
                                rng=np.random.default_rng(seed_for(cohort,taxon_id,predictor)))
            row[f"{predictor}_delta_sd"] = fit[0] if fit is not None else None
            if fit is not None:
                effects.append({"predictor": predictor, "delta": fit[0], "null": fit[1],
                                "genus": row["genus"]})
        # Sensitivity: omit yellow/orange because visible hues do not share one pigment chemistry.
        # Red/pink/blue/purple are still only photographic pigment-compatible proxies.
        rb = g.loc[g.morph.isin(("white", "red_pink", "blue_purple"))]
        rb_chrom = rb.morph.to_numpy() != "white"
        for predictor in ("abs_lat", "elevation_m"):
            fit_rb = species_delta(rb[predictor].to_numpy(float), rb_chrom,
                                   MIN_AXIS_SPAN[predictor],
                                   rng=np.random.default_rng(seed_for(cohort,taxon_id,"rb",predictor)))
            row[f"rb_{predictor}_delta_sd"] = None if fit_rb is None else fit_rb[0]
            if fit_rb is not None:
                effects.append({"predictor": f"rb_{predictor}", "delta": fit_rb[0],
                                "null": fit_rb[1], "genus": row["genus"]})
        # A noninferential collinearity check; do not promote to an independent primary test.
        alt_resid = residualize(g.elevation_m.to_numpy(float), g.abs_lat.to_numpy(float))
        sens = species_delta(alt_resid, chrom, 25.0,
                             rng=np.random.default_rng(seed_for(cohort,taxon_id,"elevation_resid")),
                             n_perm=0)
        row["elevation_resid_lat_delta_sd"] = None if sens is None else sens[0]
        row.update({"local_"+k: v for k,v in local_white_chromatic(g).items()})
        records.append(row)
    species = pd.DataFrame.from_records(records)
    lat = summarize_effect(effects, "abs_lat", cohort)
    alt = summarize_effect(effects, "elevation_m", cohort)
    if lat["estimable"] and alt["estimable"]:
        adj_lat, adj_alt = holm_two(
            lat["within_species_vertex_permutation_p_positive"],
            alt["within_species_vertex_permutation_p_positive"])
        lat["holm_two_axis_p_positive"] = adj_lat
        alt["holm_two_axis_p_positive"] = adj_alt
    local = species.loc[species.local_eligible.fillna(False)]
    wc_ratios = pd.to_numeric(local["local_white_chromatic_ratio"], errors="coerce").dropna()
    hue_ratios = pd.to_numeric(local["local_nonwhite_hue_ratio"], errors="coerce").dropna()
    # Paired within-species comparison guards against comparing two compositionally
    # distinct sets of species and local photographic opportunities.
    local_pair = local.loc[
        pd.to_numeric(local["local_white_chromatic_ratio"], errors="coerce").notna()
        & pd.to_numeric(local["local_nonwhite_hue_ratio"], errors="coerce").notna()
    ].copy()
    differences = (local_pair["local_white_chromatic_ratio"]
                   - local_pair["local_nonwhite_hue_ratio"]).to_numpy(float)
    if len(differences) >= 5:
        rng_p = np.random.default_rng(seed_for("pair-ratio", cohort))
        boot = differences[rng_p.integers(0,len(differences),
                                          size=(1999,len(differences)))].mean(axis=1)
        paired_ci = [float(v) for v in np.quantile(boot,[0.025,0.975])]
        paired_mean = float(differences.mean())
        # Only a secondary *exploratory* blocked-pair sign test, never adaptive proof.
        sign = rng_p.choice((-1,1), size=(9999,len(differences)))
        null = (sign @ differences)/len(differences)
        paired_p_two = float((1+np.count_nonzero(np.abs(null)>=abs(paired_mean)))/10000)
    else:
        paired_ci, paired_mean, paired_p_two = None, None, None
    result = {
        "n_species_eligible_global": int(len(species)),
        "n_species_white_plus_colour_at_least_5_each": int(species.white_plus_colour.sum()),
        "n_species_white_plus_one_nonwhite_hue_at_least_5_each": int(
            species.white_with_one_nonwhite_hue_at_least_5_each.sum()),
        "n_species_dominant_pair_white_nonwhite_at_least_5_each": int(
            species.dominant_top_two_white_nonwhite.sum()),
        "n_species_two_nonwhite_hues_at_least_5_each": int(species.two_nonwhite_hues.sum()),
        "n_species_both_types": int((species.white_plus_colour & species.two_nonwhite_hues).sum()),
        "n_species_white_colour_only": int((species.white_plus_colour & ~species.two_nonwhite_hues).sum()),
        "n_species_nonwhite_hue_only": int((~species.white_plus_colour & species.two_nonwhite_hues).sum()),
        "local_50km": {
            "n_evaluable_species_at_least_30_local_pairs": int(len(local)),
            "n_with_at_least_one_white_chromatic_local_pair": int(local.local_has_local_white_chromatic_pair.sum()),
            "n_with_white_chromatic_ratio_estimable": int(len(wc_ratios)),
            "mean_species_ratio_white_chromatic_observed_to_fixed_composition_expected": float(wc_ratios.mean()) if len(wc_ratios) else None,
            "median_species_ratio_white_chromatic": float(wc_ratios.median()) if len(wc_ratios) else None,
            "n_with_nonwhite_hue_ratio_estimable": int(len(hue_ratios)),
            "mean_species_ratio_nonwhite_hue_observed_to_fixed_composition_expected": float(hue_ratios.mean()) if len(hue_ratios) else None,
            "paired_species_white_colour_minus_nonwhite_hue": {
                "n_species_with_both_ratios_estimable": int(len(differences)),
                "mean_within_species_difference": paired_mean,
                "species_bootstrap_95CI": paired_ci,
                "paired_signflip_two_sided_p_exploratory": paired_p_two,
                "limitation": "coarse hues and varying expected discordances; no claim of morph-specific adaptation",
            },
        },
        "primary_high_absolute_latitude_chromatic": lat,
        "primary_high_elevation_chromatic": alt,
        "sensitivity_red_pink_blue_purple_vs_white": {
            "high_abs_latitude": summarize_effect(effects, "rb_abs_lat", cohort),
            "high_elevation": summarize_effect(effects, "rb_elevation_m", cohort),
            "role": "secondary descriptive red/blue hue comparison, not biochemical anthocyanin verification",
        },
        "sensitivity_latitude_adjusted_elevation": {
            "n_species": int(species.elevation_resid_lat_delta_sd.notna().sum()),
            "mean_species_delta_sd": float(species.elevation_resid_lat_delta_sd.mean())
            if species.elevation_resid_lat_delta_sd.notna().any() else None,
        },
    }
    return result, species


def main() -> None:
    ap = argparse.ArgumentParser()
    for cohort in INPUT_SHA256:
        ap.add_argument(f"--{cohort}", required=True, type=Path)
    ap.add_argument("--elevation", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = {
        "schema": "fcp_white_chromatic_clines_posthoc_v1",
        "date_jst": "2026-10-08",
        "role": "new_post_outcome_exploratory_ecological_screen",
        "status": "complete",
        "confirmatory_decisions_changed": False,
        "raw_data_sha256": INPUT_SHA256,
        "worldclim_elevation_10arcmin_raster_sha256": file_sha(args.elevation),
        "methods": {
            "latitude_axis": "absolute latitude (distance from equator), within species",
            "elevation_axis": "WorldClim2.1 10 arc-minute elevation, within species",
            "colour": "visible 4-class photo label; white vs all non-white labels",
            "species_inclusion": ">=40 classifiable photographs; >=5 per contrast arm; >=1 latitude-degree or >=100 elevation-metre sampled span",
            "primary_statistic": "equal-species mean of chromatic-minus-white environment difference / within-species SD",
            "null": "499 within-species fixed-count photo-label permutations; 2-axis Holm correction within cohort",
            "uncertainty": "1999 species-level bootstrap resamples; genus-balanced sensitivity",
            "local": "50 km local pair ratio versus species-wide exact fixed-colour composition; descriptive",
        },
        "hard_nonclaims": [
            "high-depth photo cohorts are not a random sample of all flowering species",
            "white detected in photos is exposure-coupled; not a validated anthocyanin deficiency",
            "visible flower colour is not UV pigment measurement",
            "the 50 km neighbourhood is not a demonstrated biological population",
            "WorldClim 10 arc-minute elevation is not field-measured altitude",
            "latitude and altitude correlations cannot distinguish climate, pollinator, history and sampling effects",
            "no plant fitness, local adaptation, adaptive maintenance or causal colour selection inferred",
            "all tests here are post-outcome and cannot upgrade prospective H2 evidence",
        ],
        "cohorts": {},
    }
    species_tables = []
    for cohort in INPUT_SHA256:
        data = read_cohort(getattr(args, cohort), cohort)
        data["elevation_m"] = sample_elevation(args.elevation, data)
        result, per_species = analyze(data, cohort)
        out["cohorts"][cohort] = result
        species_tables.append(per_species)
    pd.concat(species_tables, ignore_index=True).to_csv(args.outdir/"species_readout.csv", index=False)
    (args.outdir/"result.json").write_text(json.dumps(out, indent=2, ensure_ascii=False))
    print(json.dumps(out, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
