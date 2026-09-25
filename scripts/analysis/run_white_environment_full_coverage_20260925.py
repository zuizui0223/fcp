#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from scipy.stats import binomtest
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

MORPHS = ["white", "yellow_orange", "red_pink", "blue_purple"]
BIO_VARS = [1, 4, 5, 6, 7, 12, 14, 15, 17]
CONTRASTS = {
    "white_vs_anthocyanic": ["red_pink", "blue_purple"],
    "white_vs_all_nonwhite": ["yellow_orange", "red_pink", "blue_purple"],
    "white_vs_yellow_orange": ["yellow_orange"],
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def bool_series(s: pd.Series) -> pd.Series:
    if s.dtype == bool:
        return s.fillna(False)
    return s.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def holm_adjust(pvalues: list[float]) -> list[float]:
    out = [np.nan] * len(pvalues)
    valid = [(i, float(p)) for i, p in enumerate(pvalues) if np.isfinite(p)]
    valid.sort(key=lambda x: x[1])
    m = len(valid)
    prev = 0.0
    for rank, (i, p) in enumerate(valid):
        val = min(1.0, (m - rank) * p)
        val = max(prev, val)
        out[i] = val
        prev = val
    return out


def load_high_clip_photo_ids(high_clip_ids: Path, sealed_join_key: Path) -> set[int]:
    high = pd.read_csv(high_clip_ids, dtype={"measurement_id": str})
    join = pd.read_csv(sealed_join_key, dtype={"measurement_id": str})
    need_h = {"measurement_id"}
    need_j = {"measurement_id", "photo_id"}
    if not need_h <= set(high.columns):
        raise SystemExit(f"high-clip table missing {sorted(need_h - set(high.columns))}")
    if not need_j <= set(join.columns):
        raise SystemExit(f"sealed join key missing {sorted(need_j - set(join.columns))}")
    ids = set(high["measurement_id"].astype(str))
    photo = pd.to_numeric(join.loc[join["measurement_id"].isin(ids), "photo_id"], errors="coerce").dropna().astype("int64")
    if len(photo) != len(ids):
        raise SystemExit(f"high-clip mapping incomplete: {len(photo)} mapped for {len(ids)} measurement IDs")
    return set(photo.tolist())


def extract_worldclim(df: pd.DataFrame, worldclim_dir: Path) -> tuple[pd.DataFrame, dict[str, str]]:
    coords = list(zip(pd.to_numeric(df["longitude"], errors="coerce"), pd.to_numeric(df["latitude"], errors="coerce")))
    if any((not np.isfinite(x)) or (not np.isfinite(y)) for x, y in coords):
        raise SystemExit("nonfinite coordinates found in retained measured rows")

    hashes: dict[str, str] = {}
    out = df.copy()
    for bio in BIO_VARS:
        path = worldclim_dir / f"wc2.1_10m_bio_{bio}.tif"
        if not path.exists():
            raise SystemExit(f"missing WorldClim raster: {path}")
        hashes[path.name] = sha256_file(path)
        with rasterio.open(path) as src:
            vals = np.array([float(v[0]) for v in src.sample(coords)], dtype=float)
            if src.nodata is not None:
                vals[np.isclose(vals, float(src.nodata), equal_nan=False)] = np.nan
        out[f"bio{bio}"] = vals
    return out, hashes


def species_equal_screen(
    frame: pd.DataFrame,
    variables: list[str],
    contrast: str,
    comparison_states: list[str],
    min_per_state: int,
    n_bootstrap: int,
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    d = frame.loc[frame["morph"].eq("white") | frame["morph"].isin(comparison_states)].copy()
    d["is_white"] = d["morph"].eq("white")

    summary_rows = []
    species_rows = []

    for variable in variables:
        effects = []
        for species, g in d.groupby("species", sort=True):
            x = pd.to_numeric(g[variable], errors="coerce")
            ok = x.notna()
            g = g.loc[ok].copy()
            x = x.loc[ok].astype(float)
            n_white = int(g["is_white"].sum())
            n_comp = int((~g["is_white"]).sum())
            if n_white < min_per_state or n_comp < min_per_state:
                continue
            sd = float(x.std(ddof=0))
            if not np.isfinite(sd) or sd <= 0:
                continue
            z = (x - float(x.mean())) / sd
            effect = float(z.loc[g["is_white"]].mean() - z.loc[~g["is_white"]].mean())
            effects.append(effect)
            species_rows.append(
                {
                    "contrast": contrast,
                    "metric": variable,
                    "species": species,
                    "effect_white_minus_comparison_sd": effect,
                    "n_white": n_white,
                    "n_comparison": n_comp,
                    "n_rows": int(len(g)),
                }
            )

        effects_arr = np.asarray(effects, dtype=float)
        if len(effects_arr) == 0:
            continue
        boots = np.empty(n_bootstrap, dtype=float)
        for b in range(n_bootstrap):
            boots[b] = float(rng.choice(effects_arr, size=len(effects_arr), replace=True).mean())
        pos = int((effects_arr > 0).sum())
        neg = int((effects_arr < 0).sum())
        n_nonzero = pos + neg
        sign_p = float(binomtest(pos, n_nonzero, 0.5).pvalue) if n_nonzero else np.nan
        summary_rows.append(
            {
                "contrast": contrast,
                "metric": variable,
                "n_species": int(len(effects_arr)),
                "mean_effect_sd": float(effects_arr.mean()),
                "median_effect_sd": float(np.median(effects_arr)),
                "bootstrap_ci_low": float(np.quantile(boots, 0.025)),
                "bootstrap_ci_high": float(np.quantile(boots, 0.975)),
                "positive_species": pos,
                "negative_species": neg,
                "sign_p": sign_p,
            }
        )

    summary = pd.DataFrame(summary_rows)
    if not summary.empty:
        summary["sign_p_holm"] = holm_adjust(summary["sign_p"].tolist())
    return summary, pd.DataFrame(species_rows)


def fit_outcome_blind_pca(frame: pd.DataFrame, variables: list[str]) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    x = frame[variables].apply(pd.to_numeric, errors="coerce")
    ok = x.notna().all(axis=1)
    if int(ok.sum()) < 100:
        raise SystemExit("too few complete climate rows for PCA")
    scaler = StandardScaler().fit(x.loc[ok].to_numpy(float))
    pca = PCA(n_components=3).fit(scaler.transform(x.loc[ok].to_numpy(float)))

    projected = frame.copy()
    projected[["PC1", "PC2", "PC3"]] = np.nan
    projected.loc[ok, ["PC1", "PC2", "PC3"]] = pca.transform(scaler.transform(x.loc[ok].to_numpy(float)))

    loadings = pd.DataFrame(
        pca.components_.T,
        index=variables,
        columns=["PC1", "PC2", "PC3"],
    ).reset_index(names="metric")

    meta = {
        "explained_variance_ratio": [float(v) for v in pca.explained_variance_ratio_],
        "scaler_mean": {v: float(m) for v, m in zip(variables, scaler.mean_)},
        "scaler_scale": {v: float(s) for v, s in zip(variables, scaler.scale_)},
    }
    return projected, loadings, meta


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--measured", required=True)
    p.add_argument("--high-clip-ids", required=True)
    p.add_argument("--sealed-join-key", required=True)
    p.add_argument("--worldclim-dir", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--min-per-state", type=int, default=5)
    p.add_argument("--bootstraps", type=int, default=9999)
    p.add_argument("--seed", type=int, default=20260925)
    args = p.parse_args()

    measured_path = Path(args.measured)
    high_path = Path(args.high_clip_ids)
    join_path = Path(args.sealed_join_key)
    wc_dir = Path(args.worldclim_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    measured = pd.read_csv(measured_path, low_memory=False)
    required = {
        "species", "photo_id", "latitude", "longitude", "morph", "global_classifiable"
    }
    missing = sorted(required - set(measured.columns))
    if missing:
        raise SystemExit(f"measured table missing columns: {missing}")

    high_photo_ids = load_high_clip_photo_ids(high_path, join_path)

    keep = bool_series(measured["global_classifiable"]) & measured["morph"].astype(str).isin(MORPHS)
    frame = measured.loc[keep].copy()
    frame["photo_id"] = pd.to_numeric(frame["photo_id"], errors="raise").astype("int64")
    pre_highclip_rows = int(len(frame))
    frame = frame.loc[~frame["photo_id"].isin(high_photo_ids)].copy()
    post_highclip_rows = int(len(frame))

    frame, wc_hashes = extract_worldclim(frame, wc_dir)
    variables = [f"bio{i}" for i in BIO_VARS]
    climate_complete = frame[variables].notna().all(axis=1)
    climate_frame = frame.loc[climate_complete].copy()

    rng = np.random.default_rng(args.seed)
    all_summary = []
    all_species = []
    for contrast, comparison in CONTRASTS.items():
        summary, species = species_equal_screen(
            climate_frame,
            variables,
            contrast,
            comparison,
            args.min_per_state,
            args.bootstraps,
            rng,
        )
        all_summary.append(summary)
        all_species.append(species)

    metric_summary = pd.concat(all_summary, ignore_index=True)
    species_effects = pd.concat(all_species, ignore_index=True)
    metric_summary.to_csv(outdir / "metric_summary.csv", index=False)
    species_effects.to_csv(outdir / "species_effects.csv", index=False)

    pca_frame, pca_loadings, pca_meta = fit_outcome_blind_pca(climate_frame, variables)
    pca_loadings.to_csv(outdir / "pca_loadings.csv", index=False)

    pca_summary_parts = []
    pca_species_parts = []
    for contrast, comparison in CONTRASTS.items():
        summary, species = species_equal_screen(
            pca_frame,
            ["PC1", "PC2", "PC3"],
            contrast,
            comparison,
            args.min_per_state,
            args.bootstraps,
            rng,
        )
        pca_summary_parts.append(summary)
        pca_species_parts.append(species)
    pca_summary = pd.concat(pca_summary_parts, ignore_index=True)
    pca_species = pd.concat(pca_species_parts, ignore_index=True)
    pca_summary.to_csv(outdir / "pca_summary.csv", index=False)
    pca_species.to_csv(outdir / "pca_species_effects.csv", index=False)

    result = {
        "schema": "white_environment_full_coverage_screen_v1",
        "status": "complete",
        "date_jst": "2026-09-25",
        "inference_status": "exploratory_full_coverage_replication_after_partial_screen_opened",
        "primary_contrast": "white_vs_anthocyanic",
        "specificity_contrast": "white_vs_yellow_orange",
        "worldclim": {
            "version": "2.1",
            "resolution": "10 arc-minute",
            "variables": variables,
            "raster_sha256": wc_hashes,
        },
        "inputs": {
            "measured_sha256": sha256_file(measured_path),
            "high_clip_ids_sha256": sha256_file(high_path),
            "sealed_join_key_sha256": sha256_file(join_path),
        },
        "counts": {
            "measured_rows": int(len(measured)),
            "classifiable_four_state_rows_before_high_clip_exclusion": pre_highclip_rows,
            "high_clip_photo_ids": int(len(high_photo_ids)),
            "classifiable_four_state_rows_after_high_clip_exclusion": post_highclip_rows,
            "climate_complete_rows": int(len(climate_frame)),
            "climate_complete_species": int(climate_frame["species"].nunique()),
        },
        "analysis": {
            "min_rows_per_state_per_species": int(args.min_per_state),
            "species_bootstraps": int(args.bootstraps),
            "seed": int(args.seed),
            "species_weighting": "equal weight across species after within-species standardization",
            "multiplicity": "Holm within each contrast separately",
            "pca": pca_meta,
        },
        "hard_nonclaims": [
            "does not identify pigmented-to-white evolutionary direction",
            "does not establish anthocyanin loss as a molecular mechanism",
            "does not establish causal climatic adaptation",
            "does not test pollinator causation",
            "does not modify the frozen New Phytologist manuscript claim",
        ],
        "files": {
            "metric_summary": "metric_summary.csv",
            "species_effects": "species_effects.csv",
            "pca_summary": "pca_summary.csv",
            "pca_species_effects": "pca_species_effects.csv",
            "pca_loadings": "pca_loadings.csv",
        },
    }
    (outdir / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")

    print(metric_summary.to_string(index=False))
    print("\nPCA loadings")
    print(pca_loadings.to_string(index=False))
    print("\nPCA summary")
    print(pca_summary.to_string(index=False))
    print("\n" + json.dumps(result["counts"], indent=2))


if __name__ == "__main__":
    main()
