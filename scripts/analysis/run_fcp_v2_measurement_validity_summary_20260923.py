#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

BIO_MORPHS = ("white", "yellow_orange", "red_pink", "blue_purple")
PALETTE = (
    "white",
    "yellow",
    "orange",
    "red",
    "pink",
    "magenta",
    "purple",
    "blue",
    "bronze",
)
FRACTIONS = [f"flower_fraction_{x}" for x in PALETTE]
N_CLASSIFIABLE_MIN = 40
PRIMARY_THRESHOLD = 0.10
EPS = 1e-12


def _bool_series(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.strip().str.casefold().isin(
        {"1", "true", "yes"}
    )


def _is_classifiable(frame: pd.DataFrame) -> pd.Series:
    return (
        frame["morph"].astype(str).isin(BIO_MORPHS)
        & frame["measurement_status"].astype(str).eq(
            "classified_four_state_morph"
        )
    )


def _normalize_palette(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    mass = x.sum(axis=1)
    out = np.zeros_like(x, dtype=float)
    keep = np.isfinite(x).all(axis=1) & (mass > 0)
    out[keep] = x[keep] / mass[keep, None]
    return out


def _hellinger_rows(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    return np.linalg.norm(np.sqrt(p) - np.sqrt(q), axis=1) / math.sqrt(2.0)


def _qwhite() -> np.ndarray:
    q = np.asarray([1.0] + [-1.0 / 8.0] * 8, dtype=float)
    return q / np.linalg.norm(q)


def _four_colour_matrix(frame: pd.DataFrame) -> np.ndarray:
    """Collapse the frozen nine biological palette fractions to the original four-group vector."""
    p = _normalize_palette(frame[FRACTIONS].to_numpy(float))
    return np.column_stack(
        [
            p[:, 0],
            p[:, 1] + p[:, 2] + p[:, 8],
            p[:, 3] + p[:, 4] + p[:, 5],
            p[:, 6] + p[:, 7],
        ]
    )


def _great_circle_pairwise_km(latitude: np.ndarray, longitude: np.ndarray) -> np.ndarray:
    lat = np.deg2rad(np.asarray(latitude, dtype=float))
    lon = np.deg2rad(np.asarray(longitude, dtype=float))
    if lat.ndim != 1 or lon.ndim != 1 or lat.shape != lon.shape or len(lat) < 2:
        raise ValueError("latitude/longitude vectors are invalid")
    if np.any(~np.isfinite(lat)) or np.any(~np.isfinite(lon)):
        raise ValueError("coordinates must be finite")
    c = np.cos(lat)
    xyz = np.column_stack([c * np.cos(lon), c * np.sin(lon), np.sin(lat)])
    dot = np.clip(xyz @ xyz.T, -1.0, 1.0)
    d = np.arccos(dot) * 6371.0088
    return d[np.triu_indices(len(lat), k=1)]


def _jensen_shannon_pairwise(traits: np.ndarray) -> np.ndarray:
    p = np.asarray(traits, dtype=float)
    if p.ndim != 2 or p.shape[0] < 2 or p.shape[1] != 4:
        raise ValueError("four-colour traits must have shape (n>=2, 4)")
    if np.any(~np.isfinite(p)) or np.any(p < 0):
        raise ValueError("four-colour traits must be finite and nonnegative")
    mass = p.sum(axis=1)
    if np.any(mass <= 0):
        raise ValueError("four-colour rows must have positive mass")
    p = p / mass[:, None]
    a = p[:, None, :]
    b = p[None, :, :]
    m = 0.5 * (a + b)
    with np.errstate(divide="ignore", invalid="ignore"):
        ka = np.where(a > 0, a * np.log2(a / m), 0.0).sum(axis=2)
        kb = np.where(b > 0, b * np.log2(b / m), 0.0).sum(axis=2)
    jsd = np.clip(0.5 * (ka + kb), 0.0, 1.0)
    return jsd[np.triu_indices(len(p), k=1)]


def _spatial_rho(latitude: np.ndarray, longitude: np.ndarray, traits: np.ndarray) -> float:
    geo = _great_circle_pairwise_km(latitude, longitude)
    colour = _jensen_shannon_pairwise(traits)
    if np.ptp(geo) <= 1e-12:
        return float("nan")
    if np.ptp(colour) <= 1e-15:
        return 0.0
    return float(spearmanr(geo, colour).statistic)


def _species_spatial(frame: pd.DataFrame, condition_id: str) -> pd.DataFrame:
    x = frame.loc[
        frame["condition_id"].astype(str).eq(condition_id)
        & _is_classifiable(frame)
    ].copy()
    rows = []
    for (panel, species), g in x.groupby(["panel", "species"], sort=True):
        n = int(len(g))
        if n < N_CLASSIFIABLE_MIN:
            continue
        lat = pd.to_numeric(g["latitude"], errors="coerce").to_numpy(float)
        lon = pd.to_numeric(g["longitude"], errors="coerce").to_numpy(float)
        if not np.isfinite(lat).all() or not np.isfinite(lon).all():
            continue
        traits = _four_colour_matrix(g)
        if np.any(traits.sum(axis=1) <= 0):
            continue
        try:
            rho = _spatial_rho(lat, lon, traits)
        except ValueError:
            continue
        if not np.isfinite(rho):
            continue
        rows.append(
            {
                "condition_id": condition_id,
                "panel": str(panel),
                "species": str(species),
                "n_classifiable": n,
                "spatial_rho": float(rho),
            }
        )
    return pd.DataFrame(rows)


def _spatial_pair_summary(
    baseline: pd.DataFrame,
    other: pd.DataFrame,
    *,
    scope: str,
    condition_id: str,
) -> dict:
    a = baseline if scope == "all" else baseline.loc[baseline["panel"].eq(scope)]
    b = other if scope == "all" else other.loc[other["panel"].eq(scope)]
    pair = a[["species", "spatial_rho"]].merge(
        b[["species", "spatial_rho"]],
        on="species",
        how="inner",
        suffixes=("_baseline", "_condition"),
        validate="one_to_one",
    )
    if len(pair) == 0:
        return {"scope": scope, "condition_id": condition_id, "species": 0}
    x = pair["spatial_rho_baseline"].to_numpy(float)
    y = pair["spatial_rho_condition"].to_numpy(float)
    intercept, slope = _calibration(x, y)
    diff = y - x
    return {
        "scope": scope,
        "condition_id": condition_id,
        "species": int(len(pair)),
        "spearman_rho": (
            float(spearmanr(x, y).statistic) if len(pair) >= 2 else float("nan")
        ),
        "lin_ccc": _ccc(x, y),
        "calibration_intercept": intercept,
        "calibration_slope": slope,
        "signed_spatial_rho_change": _summary_numeric(diff),
        "absolute_spatial_rho_change": _summary_numeric(np.abs(diff)),
    }


def _D_spatial_descriptive(
    D_table: pd.DataFrame,
    spatial_table: pd.DataFrame,
    *,
    scope: str,
    condition_id: str,
) -> dict:
    d = D_table.loc[D_table["condition_id"].eq(condition_id)].copy()
    s = spatial_table.loc[spatial_table["condition_id"].eq(condition_id)].copy()
    if scope != "all":
        d = d.loc[d["panel"].eq(scope)]
        s = s.loc[s["panel"].eq(scope)]
    pair = d[["species", "D"]].merge(
        s[["species", "spatial_rho"]],
        on="species",
        how="inner",
        validate="one_to_one",
    )
    if len(pair) < 2:
        return {
            "scope": scope,
            "condition_id": condition_id,
            "species": int(len(pair)),
            "spearman_D_spatial_rho": float("nan"),
        }
    return {
        "scope": scope,
        "condition_id": condition_id,
        "species": int(len(pair)),
        "spearman_D_spatial_rho": float(
            spearmanr(
                pair["D"].to_numpy(float),
                pair["spatial_rho"].to_numpy(float),
            ).statistic
        ),
        "p_value_generated": False,
    }


def _ccc(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 2:
        return float("nan")
    vx = float(np.var(x, ddof=0))
    vy = float(np.var(y, ddof=0))
    mx = float(np.mean(x))
    my = float(np.mean(y))
    cov = float(np.mean((x - mx) * (y - my)))
    denom = vx + vy + (mx - my) ** 2
    return float(2.0 * cov / denom) if denom > 0 else float("nan")


def _calibration(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 2 or float(np.var(x)) <= 0:
        return float("nan"), float("nan")
    slope, intercept = np.polyfit(x, y, deg=1)
    return float(intercept), float(slope)


def _summary_numeric(values: pd.Series | np.ndarray) -> dict[str, float | int]:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return {
            "n": 0,
            "mean": float("nan"),
            "median": float("nan"),
            "q025": float("nan"),
            "q975": float("nan"),
        }
    return {
        "n": int(len(x)),
        "mean": float(np.mean(x)),
        "median": float(np.median(x)),
        "q025": float(np.quantile(x, 0.025)),
        "q975": float(np.quantile(x, 0.975)),
    }


def _condition_pairs(
    biological: pd.DataFrame,
    join_key: pd.DataFrame,
    condition_id: str,
) -> pd.DataFrame:
    base = biological.loc[
        biological["condition_id"].astype(str).eq("baseline")
    ].copy()
    cf = biological.loc[
        biological["condition_id"].astype(str).eq(condition_id)
    ].copy()
    if base["measurement_id"].duplicated().any() or cf["measurement_id"].duplicated().any():
        raise RuntimeError("condition table contains duplicate measurement IDs")

    base = base.merge(
        join_key[["measurement_id", "panel", "species"]],
        on="measurement_id",
        how="inner",
        validate="one_to_one",
    )
    cf_keep = [
        "measurement_id",
        "morph",
        "measurement_status",
        "counterfactual_status",
        *FRACTIONS,
    ]
    cf = cf[cf_keep].rename(
        columns={
            "morph": "cf_morph",
            "measurement_status": "cf_measurement_status",
            "counterfactual_status": "cf_counterfactual_status",
            **{c: f"cf_{c}" for c in FRACTIONS},
        }
    )
    pairs = base.merge(
        cf,
        on="measurement_id",
        how="inner",
        validate="one_to_one",
    )
    return pairs


def _mv1_for_scope(pairs: pd.DataFrame, scope: str) -> tuple[dict, pd.DataFrame]:
    x = pairs if scope == "all" else pairs.loc[pairs["panel"].eq(scope)]
    if len(x) == 0:
        return {"scope": scope, "pairs": 0}, pd.DataFrame()

    base_class = (
        x["morph"].astype(str).isin(BIO_MORPHS)
        & x["measurement_status"].astype(str).eq(
            "classified_four_state_morph"
        )
    )
    cf_class = (
        x["cf_morph"].astype(str).isin(BIO_MORPHS)
        & x["cf_measurement_status"].astype(str).eq(
            "classified_four_state_morph"
        )
    )
    disagreement = x["morph"].astype(str).ne(x["cf_morph"].astype(str))
    both = base_class & cf_class
    flip = disagreement & both

    p0 = _normalize_palette(x[FRACTIONS].to_numpy(float))
    p1 = _normalize_palette(
        x[[f"cf_{c}" for c in FRACTIONS]].to_numpy(float)
    )
    valid_palette = (p0.sum(axis=1) > 0) & (p1.sum(axis=1) > 0)
    hellinger = np.full(len(x), np.nan)
    if valid_palette.any():
        hellinger[valid_palette] = _hellinger_rows(
            p0[valid_palette], p1[valid_palette]
        )

    per_species = pd.DataFrame(
        {
            "species": x["species"].astype(str).to_numpy(),
            "disagreement": disagreement.astype(float).to_numpy(),
            "both_classifiable": both.astype(float).to_numpy(),
            "flip": flip.astype(float).to_numpy(),
            "hellinger": hellinger,
            "base_classifiable": base_class.astype(float).to_numpy(),
            "cf_classifiable": cf_class.astype(float).to_numpy(),
        }
    ).groupby("species", sort=True).agg(
        rows=("disagreement", "size"),
        complete_state_disagreement=("disagreement", "mean"),
        both_classifiable_rows=("both_classifiable", "sum"),
        morph_flip_joint=("flip", "sum"),
        mean_hellinger=("hellinger", "mean"),
        baseline_classifiable_fraction=("base_classifiable", "mean"),
        counterfactual_classifiable_fraction=("cf_classifiable", "mean"),
    ).reset_index()
    per_species["morph_flip_rate_both_classifiable"] = np.where(
        per_species["both_classifiable_rows"] > 0,
        per_species["morph_flip_joint"]
        / per_species["both_classifiable_rows"],
        np.nan,
    )
    per_species["classifiability_change"] = (
        per_species["counterfactual_classifiable_fraction"]
        - per_species["baseline_classifiable_fraction"]
    )

    summary = {
        "scope": scope,
        "pairs": int(len(x)),
        "species": int(x["species"].nunique()),
        "complete_state_disagreement_image": float(disagreement.mean()),
        "both_classifiable_rows": int(both.sum()),
        "morph_flip_rate_both_classifiable_image": (
            float(flip.sum() / both.sum()) if both.sum() else float("nan")
        ),
        "baseline_classifiable_fraction_image": float(base_class.mean()),
        "counterfactual_classifiable_fraction_image": float(cf_class.mean()),
        "classifiability_change_image": float(cf_class.mean() - base_class.mean()),
        "hellinger_image": _summary_numeric(hellinger),
        "equal_species_complete_state_disagreement": _summary_numeric(
            per_species["complete_state_disagreement"]
        ),
        "equal_species_morph_flip_rate": _summary_numeric(
            per_species["morph_flip_rate_both_classifiable"]
        ),
        "equal_species_hellinger": _summary_numeric(
            per_species["mean_hellinger"]
        ),
        "equal_species_classifiability_change": _summary_numeric(
            per_species["classifiability_change"]
        ),
    }
    per_species.insert(0, "scope", scope)
    return summary, per_species


def _transition_matrix(pairs: pd.DataFrame, condition_id: str) -> pd.DataFrame:
    labels = (*BIO_MORPHS, "mixed_uncertain")
    base = pairs["morph"].astype(str).where(
        pairs["morph"].astype(str).isin(labels),
        "mixed_uncertain",
    )
    cf = pairs["cf_morph"].astype(str).where(
        pairs["cf_morph"].astype(str).isin(labels),
        "mixed_uncertain",
    )
    table = (
        pd.crosstab(base, cf, dropna=False)
        .reindex(index=labels, columns=labels, fill_value=0)
        .stack()
        .rename("rows")
        .reset_index()
        .rename(columns={"morph": "baseline_state", "cf_morph": "counterfactual_state"})
    )
    table.insert(0, "condition_id", condition_id)
    return table


def _mv3a_for_scope(
    pairs: pd.DataFrame,
    scope: str,
) -> tuple[dict, pd.DataFrame]:
    x = pairs if scope == "all" else pairs.loc[pairs["panel"].eq(scope)]
    if len(x) == 0:
        return {"scope": scope, "pairs": 0}, pd.DataFrame()
    p0 = _normalize_palette(x[FRACTIONS].to_numpy(float))
    p1 = _normalize_palette(
        x[[f"cf_{c}" for c in FRACTIONS]].to_numpy(float)
    )
    d = p1 - p0
    norms = np.linalg.norm(d, axis=1)
    valid = (
        (p0.sum(axis=1) > 0)
        & (p1.sum(axis=1) > 0)
        & np.isfinite(norms)
        & (norms > EPS)
    )
    q = _qwhite()
    t_white = np.full(len(x), np.nan)
    if valid.any():
        u = d[valid] / norms[valid, None]
        t_white[valid] = (u @ q) ** 2
    white_change = p1[:, 0] - p0[:, 0]
    white_change[~((p0.sum(axis=1) > 0) & (p1.sum(axis=1) > 0))] = np.nan

    per = pd.DataFrame(
        {
            "species": x["species"].astype(str).to_numpy(),
            "technical_displacement_nonzero": valid.astype(float),
            "T_white": t_white,
            "white_fraction_change": white_change,
        }
    ).groupby("species", sort=True).agg(
        rows=("technical_displacement_nonzero", "size"),
        nonzero_fraction=("technical_displacement_nonzero", "mean"),
        mean_T_white=("T_white", "mean"),
        mean_white_fraction_change=("white_fraction_change", "mean"),
    ).reset_index()

    summary = {
        "scope": scope,
        "pairs": int(len(x)),
        "species": int(x["species"].nunique()),
        "nonzero_technical_displacements": int(valid.sum()),
        "nonzero_fraction_image": float(valid.mean()),
        "T_white_image": _summary_numeric(t_white),
        "white_fraction_change_image": _summary_numeric(white_change),
        "equal_species_nonzero_fraction": _summary_numeric(
            per["nonzero_fraction"]
        ),
        "equal_species_T_white": _summary_numeric(per["mean_T_white"]),
        "equal_species_white_fraction_change": _summary_numeric(
            per["mean_white_fraction_change"]
        ),
    }
    per.insert(0, "scope", scope)
    return summary, per


def _species_D(frame: pd.DataFrame, condition_id: str) -> pd.DataFrame:
    x = frame.loc[
        frame["condition_id"].astype(str).eq(condition_id)
        & _is_classifiable(frame)
    ].copy()
    rows = []
    for (panel, species), g in x.groupby(["panel", "species"], sort=True):
        n = int(len(g))
        if n < N_CLASSIFIABLE_MIN:
            continue
        counts = (
            g["morph"].astype(str).value_counts()
            .reindex(BIO_MORPHS, fill_value=0)
            .to_numpy(float)
        )
        p = counts / counts.sum()
        D = float(1.0 - np.sum(p * p))
        rows.append(
            {
                "condition_id": condition_id,
                "panel": str(panel),
                "species": str(species),
                "n_classifiable": n,
                "D": D,
            }
        )
    return pd.DataFrame(rows)


def _D_pair_summary(
    baseline: pd.DataFrame,
    other: pd.DataFrame,
    *,
    scope: str,
    condition_id: str,
) -> dict:
    a = baseline if scope == "all" else baseline.loc[baseline["panel"].eq(scope)]
    b = other if scope == "all" else other.loc[other["panel"].eq(scope)]
    pair = a[["species", "D"]].merge(
        b[["species", "D"]],
        on="species",
        how="inner",
        suffixes=("_baseline", "_counterfactual"),
        validate="one_to_one",
    )
    if len(pair) == 0:
        return {"scope": scope, "condition_id": condition_id, "species": 0}
    xb = pair["D_baseline"].to_numpy(float)
    yc = pair["D_counterfactual"].to_numpy(float)
    rho = float(spearmanr(xb, yc).statistic) if len(pair) >= 2 else float("nan")
    intercept, slope = _calibration(xb, yc)
    diff = yc - xb
    return {
        "scope": scope,
        "condition_id": condition_id,
        "species": int(len(pair)),
        "spearman_rho": rho,
        "lin_ccc": _ccc(xb, yc),
        "calibration_intercept": intercept,
        "calibration_slope": slope,
        "signed_D_change": _summary_numeric(diff),
        "absolute_D_change": _summary_numeric(np.abs(diff)),
    }


def _deterministic_two_means(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(x, dtype=float)
    grand = x.mean(axis=0)
    i0 = int(np.argmax(np.sum((x - grand) ** 2, axis=1)))
    d0 = np.sum((x - x[i0]) ** 2, axis=1)
    i1 = int(np.argmax(d0))
    if float(d0[i1]) <= EPS:
        labels = np.zeros(len(x), dtype=int)
        labels[len(x) // 2 :] = 1
    else:
        centers = np.vstack([x[i0], x[i1]])
        labels = np.full(len(x), -1, dtype=int)
        for _ in range(200):
            dist = np.sum(
                (x[:, None, :] - centers[None, :, :]) ** 2,
                axis=2,
            )
            new = np.argmin(dist, axis=1).astype(int)
            if np.all(new == new[0]):
                only = int(new[0])
                new[int(np.argmax(dist[:, only]))] = 1 - only
            if np.array_equal(new, labels):
                break
            labels = new
            for k in (0, 1):
                members = x[labels == k]
                if len(members) == 0:
                    raise RuntimeError("empty two-means cluster")
                centers[k] = members.mean(axis=0)
        else:
            raise RuntimeError("two-means did not converge")
    return labels, np.bincount(labels, minlength=2)


def _coarse_second_fraction(g: pd.DataFrame) -> float:
    c = Counter(g["morph"].astype(str))
    order = sorted(BIO_MORPHS, key=lambda m: (-c.get(m, 0), m))
    return float(c.get(order[1], 0) / len(g))


def _vector_table(
    frame: pd.DataFrame,
    condition_id: str,
    scope: str,
) -> pd.DataFrame:
    x = frame.loc[
        frame["condition_id"].astype(str).eq(condition_id)
        & _is_classifiable(frame)
    ].copy()
    if scope != "all":
        x = x.loc[x["panel"].eq(scope)]
    rows = []
    q = _qwhite()
    for species, g in x.groupby("species", sort=True):
        n = int(len(g))
        if n < N_CLASSIFIABLE_MIN:
            continue
        if _coarse_second_fraction(g) < PRIMARY_THRESHOLD:
            continue
        p = _normalize_palette(g[FRACTIONS].to_numpy(float))
        if np.any(p.sum(axis=1) <= 0):
            continue
        labels, counts = _deterministic_two_means(np.sqrt(p))
        if counts.min() / n < PRIMARY_THRESHOLD:
            continue
        if counts[0] > counts[1]:
            major, minor = 0, 1
        elif counts[1] > counts[0]:
            major, minor = 1, 0
        else:
            c0 = p[labels == 0].mean(axis=0)
            c1 = p[labels == 1].mean(axis=0)
            major, minor = (
                (0, 1)
                if tuple(c0.tolist()) <= tuple(c1.tolist())
                else (1, 0)
            )
        delta = (
            p[labels == minor].mean(axis=0)
            - p[labels == major].mean(axis=0)
        )
        norm = float(np.linalg.norm(delta))
        if not np.isfinite(norm) or norm <= EPS:
            continue
        u = delta / norm
        rows.append(
            {
                "condition_id": condition_id,
                "scope": scope,
                "species": str(species),
                "n_classifiable": n,
                "minor_cluster_fraction": float(counts.min() / n),
                "projection_sq": float((u @ q) ** 2),
                **{
                    f"u_{PALETTE[j]}": float(u[j])
                    for j in range(len(PALETTE))
                },
            }
        )
    return pd.DataFrame(rows)


def _W_pair_summary(
    base: pd.DataFrame,
    other: pd.DataFrame,
    *,
    scope: str,
    condition_id: str,
) -> dict:
    if len(base) == 0 or len(other) == 0:
        return {
            "scope": scope,
            "condition_id": condition_id,
            "baseline_vector_species": int(len(base)),
            "condition_vector_species": int(len(other)),
            "common_vector_species": 0,
        }
    W_base = float(base["projection_sq"].mean())
    W_other = float(other["projection_sq"].mean())
    unit_cols = [f"u_{x}" for x in PALETTE]
    pair = base[["species", "projection_sq", *unit_cols]].merge(
        other[["species", "projection_sq", *unit_cols]],
        on="species",
        how="inner",
        suffixes=("_baseline", "_counterfactual"),
        validate="one_to_one",
    )
    if len(pair) == 0:
        return {
            "scope": scope,
            "condition_id": condition_id,
            "baseline_vector_species": int(len(base)),
            "condition_vector_species": int(len(other)),
            "baseline_W_own_set": W_base,
            "condition_W_own_set": W_other,
            "common_vector_species": 0,
        }
    base_proj = pair["projection_sq_baseline"].to_numpy(float)
    cf_proj = pair["projection_sq_counterfactual"].to_numpy(float)
    ub = pair[[f"{c}_baseline" for c in unit_cols]].to_numpy(float)
    uc = pair[[f"{c}_counterfactual" for c in unit_cols]].to_numpy(float)
    cosine = np.sum(ub * uc, axis=1)
    return {
        "scope": scope,
        "condition_id": condition_id,
        "baseline_vector_species": int(len(base)),
        "condition_vector_species": int(len(other)),
        "baseline_W_own_set": W_base,
        "condition_W_own_set": W_other,
        "W_own_set_change": float(W_other - W_base),
        "common_vector_species": int(len(pair)),
        "baseline_W_common": float(np.mean(base_proj)),
        "condition_W_common": float(np.mean(cf_proj)),
        "W_common_change": float(np.mean(cf_proj) - np.mean(base_proj)),
        "projection_contribution_change": _summary_numeric(
            cf_proj - base_proj
        ),
        "delta_vector_cosine": _summary_numeric(cosine),
    }


def _make_baseline_stratum_conditions(
    biological: pd.DataFrame,
    strata: pd.DataFrame,
) -> pd.DataFrame:
    base = biological.loc[
        biological["condition_id"].astype(str).eq("baseline")
    ].copy()
    keep_cols = [
        "measurement_id",
        "stratum_high_flower_near_clip",
        "stratum_high_background_near_clip",
        "stratum_low_flower_background_lab_separation",
        "stratum_roi_unstable_fixed",
    ]
    missing = [c for c in keep_cols if c not in strata.columns]
    if missing:
        raise RuntimeError(f"technical strata missing columns: {missing}")
    x = base.merge(
        strata[keep_cols],
        on="measurement_id",
        how="inner",
        validate="one_to_one",
    )
    for c in keep_cols[1:]:
        x[c] = _bool_series(x[c])

    specs = {
        "baseline_excl_high_flower_near_clip": ~x[
            "stratum_high_flower_near_clip"
        ],
        "baseline_excl_high_background_near_clip": ~x[
            "stratum_high_background_near_clip"
        ],
        "baseline_excl_low_lab_separation": ~x[
            "stratum_low_flower_background_lab_separation"
        ],
        "baseline_excl_roi_unstable": ~x[
            "stratum_roi_unstable_fixed"
        ],
        "baseline_excl_union_primary_technical": ~(
            x["stratum_high_flower_near_clip"]
            | x["stratum_high_background_near_clip"]
            | x["stratum_low_flower_background_lab_separation"]
            | x["stratum_roi_unstable_fixed"]
        ),
    }
    outputs = []
    base_columns = biological.columns
    for name, mask in specs.items():
        y = x.loc[mask, base_columns].copy()
        y["condition_id"] = name
        y["condition_family"] = "technical_stratum_exclusion"
        y["condition_value"] = name
        outputs.append(y)
    return pd.concat(outputs, ignore_index=True) if outputs else biological.iloc[0:0].copy()


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--biological-long", type=Path, required=True)
    p.add_argument("--terminal-receipt", type=Path, required=True)
    p.add_argument("--sealed-join-key", type=Path, required=True)
    p.add_argument("--technical-strata", type=Path, required=True)
    p.add_argument("--prior-d", type=Path)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()

    biological = pd.read_csv(
        args.biological_long,
        dtype={"measurement_id": str},
        low_memory=False,
    )
    terminal = pd.read_csv(
        args.terminal_receipt,
        dtype={"measurement_id": str},
    )
    join_key = pd.read_csv(
        args.sealed_join_key,
        dtype={"measurement_id": str},
    )
    strata = pd.read_csv(
        args.technical_strata,
        dtype={"measurement_id": str},
    )

    if len(terminal) != 40_000 or terminal["measurement_id"].nunique() != 40_000:
        raise RuntimeError("Pass-B terminal denominator is not 40,000 unique IDs")
    if len(join_key) != 40_000 or join_key["measurement_id"].nunique() != 40_000:
        raise RuntimeError("sealed join denominator is not 40,000 unique IDs")
    if len(strata) != 40_000 or strata["measurement_id"].nunique() != 40_000:
        raise RuntimeError("technical strata denominator is not 40,000 unique IDs")
    if set(terminal["measurement_id"]) != set(join_key["measurement_id"]):
        raise RuntimeError("terminal receipt and join key IDs differ")
    if set(strata["measurement_id"]) != set(join_key["measurement_id"]):
        raise RuntimeError("technical strata and join key IDs differ")

    required_join = {"measurement_id", "panel", "species", "latitude", "longitude"}
    if not required_join.issubset(join_key.columns):
        raise RuntimeError("sealed join key lacks species/panel")
    if join_key["species"].nunique() != 400:
        raise RuntimeError("species denominator drift after Pass-B opening")
    if set(join_key["panel"].astype(str)) != {"P", "N"}:
        raise RuntimeError("panel labels drifted")

    source_counts = terminal["pass_b_terminal_status"].astype(str).value_counts()
    source_match_ids = set(
        terminal.loc[
            terminal["source_byte_match"].pipe(_bool_series),
            "measurement_id",
        ].astype(str)
    )
    if not set(biological["measurement_id"].astype(str)) <= source_match_ids:
        raise RuntimeError("biological output contains non-byte-matched rows")

    biological = biological.merge(
        join_key[["measurement_id", "panel", "species", "latitude", "longitude"]],
        on="measurement_id",
        how="inner",
        validate="many_to_one",
    )
    stratum_conditions = _make_baseline_stratum_conditions(
        biological, strata
    )
    analysis_frame = pd.concat(
        [biological, stratum_conditions],
        ignore_index=True,
    )

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)

    # MV1 and MV3a same-image counterfactuals.
    condition_ids = sorted(
        c
        for c in biological["condition_id"].astype(str).unique()
        if c != "baseline"
    )
    mv1_rows = []
    mv1_species = []
    transitions = []
    mv3a_rows = []
    mv3a_species = []
    for condition_id in condition_ids:
        pairs = _condition_pairs(
            biological.drop(columns=["panel", "species"]),
            join_key,
            condition_id,
        )
        transitions.append(_transition_matrix(pairs, condition_id))
        for scope in ("all", "P", "N"):
            s1, ps1 = _mv1_for_scope(pairs, scope)
            s1["condition_id"] = condition_id
            mv1_rows.append(s1)
            if len(ps1):
                ps1.insert(1, "condition_id", condition_id)
                mv1_species.append(ps1)
            s3, ps3 = _mv3a_for_scope(pairs, scope)
            s3["condition_id"] = condition_id
            mv3a_rows.append(s3)
            if len(ps3):
                ps3.insert(1, "condition_id", condition_id)
                mv3a_species.append(ps3)

    pd.DataFrame(mv1_rows).to_json(
        out / "mv1_image_invariance_summary.json",
        orient="records",
        indent=2,
    )
    if mv1_species:
        pd.concat(mv1_species, ignore_index=True).to_csv(
            out / "mv1_equal_species.csv",
            index=False,
            lineterminator="\n",
        )
    if transitions:
        pd.concat(transitions, ignore_index=True).to_csv(
            out / "mv1_transition_matrices.csv",
            index=False,
            lineterminator="\n",
        )
    pd.DataFrame(mv3a_rows).to_json(
        out / "mv3a_technical_qwhite_summary.json",
        orient="records",
        indent=2,
    )
    if mv3a_species:
        pd.concat(mv3a_species, ignore_index=True).to_csv(
            out / "mv3a_equal_species.csv",
            index=False,
            lineterminator="\n",
        )

    # MV2 D invariance for baseline, all-row fixed-mask conditions and frozen
    # technical-stratum exclusions.
    allrow_conditions = ["baseline"] + sorted(
        biological.loc[
            biological["condition_family"].astype(str).eq(
                "fixed_mask_exposure"
            ),
            "condition_id",
        ].astype(str).unique().tolist()
    ) + sorted(
        stratum_conditions["condition_id"].astype(str).unique().tolist()
    )
    D_tables = []
    for condition_id in allrow_conditions:
        table = _species_D(analysis_frame, condition_id)
        if len(table):
            D_tables.append(table)
    D_all = (
        pd.concat(D_tables, ignore_index=True)
        if D_tables
        else pd.DataFrame()
    )
    D_all.to_csv(
        out / "mv2_species_D_by_condition.csv",
        index=False,
        lineterminator="\n",
    )
    D_summaries = []
    baseline_D = D_all.loc[D_all["condition_id"].eq("baseline")]
    for condition_id in allrow_conditions:
        if condition_id == "baseline":
            continue
        other = D_all.loc[D_all["condition_id"].eq(condition_id)]
        for scope in ("all", "P", "N"):
            D_summaries.append(
                _D_pair_summary(
                    baseline_D,
                    other,
                    scope=scope,
                    condition_id=condition_id,
                )
            )
    pd.DataFrame(D_summaries).to_json(
        out / "mv2_D_invariance_summary.json",
        orient="records",
        indent=2,
    )

    prior_transport = None
    if args.prior_d is not None:
        prior = pd.read_csv(args.prior_d)
        if not {"species", "D"}.issubset(prior.columns):
            raise RuntimeError("prior D table schema drift")
        pbase = baseline_D.loc[
            baseline_D["panel"].eq("P"),
            ["species", "D"],
        ]
        pair = pbase.merge(
            prior[["species", "D"]].rename(columns={"D": "D_prior"}),
            on="species",
            how="inner",
            validate="one_to_one",
        ).rename(columns={"D": "D_fresh"})
        if len(pair):
            x = pair["D_prior"].to_numpy(float)
            y = pair["D_fresh"].to_numpy(float)
            intercept, slope = _calibration(x, y)
            prior_transport = {
                "species": int(len(pair)),
                "spearman_rho": float(spearmanr(x, y).statistic)
                if len(pair) >= 2
                else float("nan"),
                "lin_ccc": _ccc(x, y),
                "calibration_intercept": intercept,
                "calibration_slope": slope,
                "signed_D_change_fresh_minus_prior": _summary_numeric(y - x),
                "absolute_D_change": _summary_numeric(np.abs(y - x)),
            }
            pair.to_csv(
                out / "mv2_panel_P_prior_D_pairs.csv",
                index=False,
                lineterminator="\n",
            )

    # MV3b q_white/W invariance under all-row fixed-mask conditions and frozen
    # technical-stratum exclusions. No new null distribution is generated.
    vector_tables = []
    W_summaries = []
    for scope in ("all", "P", "N"):
        base_vec = _vector_table(
            analysis_frame,
            "baseline",
            scope,
        )
        if len(base_vec):
            vector_tables.append(base_vec)
        for condition_id in allrow_conditions:
            if condition_id == "baseline":
                continue
            other_vec = _vector_table(
                analysis_frame,
                condition_id,
                scope,
            )
            if len(other_vec):
                vector_tables.append(other_vec)
            W_summaries.append(
                _W_pair_summary(
                    base_vec,
                    other_vec,
                    scope=scope,
                    condition_id=condition_id,
                )
            )
    vector_all = (
        pd.concat(vector_tables, ignore_index=True)
        if vector_tables
        else pd.DataFrame()
    )
    vector_all.to_csv(
        out / "mv3b_species_vectors_by_condition.csv",
        index=False,
        lineterminator="\n",
    )
    pd.DataFrame(W_summaries).to_json(
        out / "mv3b_W_invariance_summary.json",
        orient="records",
        indent=2,
    )

    # MV4 spatial measurement sensitivity. This reuses the frozen FCP/RGFCA
    # species statistic: all unordered photo pairs, great-circle distance,
    # Jensen-Shannon divergence of the continuous four-group colour vectors,
    # classifiable rows only, >=40 rows/species. No new permutation null.
    spatial_tables = []
    spatial_summaries = []
    for condition_id in allrow_conditions:
        table = _species_spatial(analysis_frame, condition_id)
        if len(table):
            spatial_tables.append(table)
    spatial_all = (
        pd.concat(spatial_tables, ignore_index=True)
        if spatial_tables
        else pd.DataFrame()
    )
    spatial_all.to_csv(
        out / "mv4_species_spatial_by_condition.csv",
        index=False,
        lineterminator="\n",
    )
    baseline_spatial = (
        spatial_all.loc[spatial_all["condition_id"].eq("baseline")]
        if len(spatial_all)
        else pd.DataFrame()
    )
    for condition_id in allrow_conditions:
        if condition_id == "baseline":
            continue
        other = (
            spatial_all.loc[spatial_all["condition_id"].eq(condition_id)]
            if len(spatial_all)
            else pd.DataFrame()
        )
        for scope in ("all", "P", "N"):
            spatial_summaries.append(
                _spatial_pair_summary(
                    baseline_spatial,
                    other,
                    scope=scope,
                    condition_id=condition_id,
                )
            )
    pd.DataFrame(spatial_summaries).to_json(
        out / "mv4_spatial_invariance_summary.json",
        orient="records",
        indent=2,
    )

    D_spatial_rows = []
    if len(spatial_all) and len(D_all):
        for condition_id in allrow_conditions:
            for scope in ("all", "P", "N"):
                D_spatial_rows.append(
                    _D_spatial_descriptive(
                        D_all,
                        spatial_all,
                        scope=scope,
                        condition_id=condition_id,
                    )
                )
    pd.DataFrame(D_spatial_rows).to_json(
        out / "mv4_D_spatial_descriptive.json",
        orient="records",
        indent=2,
    )

    result = {
        "schema": "fcp_v2_measurement_validity_summary_v1",
        "status": "COMPLETE_MV1_MV2_MV3_MV4_WITHOUT_NEW_NULL",
        "source_identity": {
            "terminal_rows": int(len(terminal)),
            "source_byte_matched_rows": int(
                terminal["source_byte_match"].pipe(_bool_series).sum()
            ),
            "terminal_status_counts": {
                str(k): int(v) for k, v in source_counts.to_dict().items()
            },
        },
        "biological_rows": int(len(biological)),
        "species": int(join_key["species"].nunique()),
        "panels": {
            str(k): int(v)
            for k, v in join_key.groupby("panel")["species"].nunique().to_dict().items()
        },
        "mv1_conditions": condition_ids,
        "mv2_conditions": allrow_conditions,
        "mv3_q_white": _qwhite().tolist(),
        "mv4_conditions": allrow_conditions,
        "mv4_spatial_definition": "Spearman(all-pair great-circle distance, JSD of frozen four-group continuous colour vector), classifiable rows only, >=40/species",
        "new_null_distribution_generated": False,
        "prior_D_transport": prior_transport,
        "hard_boundaries": [
            "No q_white refit.",
            "No new structured null family.",
            "Heavy 20-row/species counterfactuals are not promoted to the original H2 W test.",
            "Technical strata were frozen before biological opening.",
            "MV4 generates no new spatial permutation null or p-value.",
            "Heavy 20-row/species perturbations are not used for primary MV4 spatial rho.",
        ],
    }
    (out / "result.json").write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=True))


if __name__ == "__main__":
    main()
