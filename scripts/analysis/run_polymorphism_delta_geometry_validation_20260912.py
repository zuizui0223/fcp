#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import linalg

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "polymorphism_delta_geometry_validation_20260912"
DISCOVERY = ROOT / "data" / "derived" / "global_monte_carlo_measured_photos_v1.csv"
RESERVE = ROOT / "data" / "derived" / "rgfca_reserve_replication_measured_photos_v1.csv"

BIO = ["white", "yellow", "orange", "red", "pink", "magenta", "purple", "blue", "bronze"]
FRACTIONS = [f"flower_fraction_{c}" for c in BIO]
MORPHS = ["white", "yellow_orange", "red_pink", "blue_purple"]
MIN_CLASSIFIABLE = 40
THRESHOLDS = {"primary_0_10": 0.10, "strict_0_20": 0.20}
N_NULL = 10_000
SEED = 20260912
EPS = 1e-12
MAX_KMEANS_ITER = 200


def bool_series(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.strip().str.lower().isin(["true", "1", "yes"])


def normalize_rows(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(x, dtype=float)
    finite = np.isfinite(x).all(axis=1)
    sums = np.nansum(x, axis=1)
    good = finite & (sums > 0)
    out = np.full_like(x, np.nan, dtype=float)
    out[good] = x[good] / sums[good, None]
    return out, good


def load_classifiable(path: Path) -> pd.DataFrame:
    usecols = ["species", "morph", "global_classifiable"] + FRACTIONS
    df = pd.read_csv(path, usecols=usecols)
    df["species"] = df["species"].astype(str)
    df["morph"] = df["morph"].fillna("").astype(str)
    keep = bool_series(df["global_classifiable"]) & df["morph"].isin(MORPHS)
    work = df.loc[keep].copy()
    x, good = normalize_rows(work[FRACTIONS].to_numpy(float))
    if not good.all():
        raise RuntimeError(
            f"classifiable rows with invalid biological palette in {path}: {int((~good).sum())}"
        )
    work.loc[:, FRACTIONS] = x
    return work


def nclass40_count(work: pd.DataFrame) -> int:
    sizes = work.groupby("species", sort=False).size()
    return int((sizes >= MIN_CLASSIFIABLE).sum())


def h1_coarse_second_fraction(g: pd.DataFrame) -> tuple[float, str, str, dict[str, int]]:
    n = int(len(g))
    cnt = Counter(g["morph"].astype(str))
    order = sorted(MORPHS, key=lambda m: (-cnt.get(m, 0), m))
    primary, secondary = order[:2]
    second_fraction = float(cnt.get(secondary, 0) / n)
    return second_fraction, primary, secondary, {m: int(cnt.get(m, 0)) for m in MORPHS}


def deterministic_two_means(x: np.ndarray) -> tuple[np.ndarray, np.ndarray, float | None]:
    """Two-means on continuous Hellinger coordinates without using morph labels."""
    x = np.asarray(x, dtype=float)
    if x.ndim != 2 or len(x) < 2:
        raise ValueError("two-means requires at least two rows")
    if not np.isfinite(x).all():
        raise ValueError("two-means received non-finite values")

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
        for _ in range(MAX_KMEANS_ITER):
            dist = np.sum((x[:, None, :] - centers[None, :, :]) ** 2, axis=2)
            new_labels = np.argmin(dist, axis=1).astype(int)
            if np.all(new_labels == new_labels[0]):
                only = int(new_labels[0])
                other = 1 - only
                far = int(np.argmax(dist[:, only]))
                new_labels[far] = other
            if np.array_equal(new_labels, labels):
                break
            labels = new_labels
            for k in (0, 1):
                members = x[labels == k]
                if len(members) == 0:
                    raise RuntimeError("deterministic two-means produced an empty cluster")
                centers[k] = members.mean(axis=0)
        else:
            raise RuntimeError("deterministic two-means did not converge")

    counts = np.bincount(labels, minlength=2)
    centers_final = np.vstack([x[labels == k].mean(axis=0) for k in (0, 1)])
    between = float(np.linalg.norm(centers_final[0] - centers_final[1]))
    within_ss = float(
        sum(np.sum((x[labels == k] - centers_final[k]) ** 2) for k in (0, 1))
    )
    within_rms = math.sqrt(within_ss / len(x)) if within_ss > EPS else 0.0
    separation_ratio = None if within_rms <= EPS else float(between / within_rms)
    return labels, counts, separation_ratio


def normalized_mutual_information(labels: np.ndarray, morphs: np.ndarray) -> float | None:
    labels = np.asarray(labels, dtype=int)
    morphs = np.asarray(morphs, dtype=str)
    n = len(labels)
    if n == 0:
        return None

    unique_morphs = sorted(set(morphs.tolist()))
    morph_to_j = {m: j for j, m in enumerate(unique_morphs)}
    table = np.zeros((2, len(unique_morphs)), dtype=float)
    for a, m in zip(labels, morphs):
        table[int(a), morph_to_j[str(m)]] += 1.0
    p = table / n
    pa = p.sum(axis=1)
    pb = p.sum(axis=0)

    mi = 0.0
    for i in range(p.shape[0]):
        for j in range(p.shape[1]):
            if p[i, j] > 0 and pa[i] > 0 and pb[j] > 0:
                mi += float(p[i, j] * math.log(p[i, j] / (pa[i] * pb[j])))
    ha = -float(sum(q * math.log(q) for q in pa if q > 0))
    hb = -float(sum(q * math.log(q) for q in pb if q > 0))
    denom = math.sqrt(ha * hb) if ha > 0 and hb > 0 else 0.0
    return None if denom <= EPS else float(mi / denom)


def species_delta_vectors(
    work: pd.DataFrame, second_min: float
) -> tuple[pd.DataFrame, np.ndarray, dict[str, int]]:
    rows: list[dict[str, object]] = []
    vectors: list[np.ndarray] = []
    gate_counts = {
        "nclass_ge40": 0,
        "H1_coarse_second_pass": 0,
        "continuous_minor_fraction_pass": 0,
        "nonzero_delta_pass": 0,
    }

    for sp, g in work.groupby("species", sort=True):
        n = int(len(g))
        if n < MIN_CLASSIFIABLE:
            continue
        gate_counts["nclass_ge40"] += 1

        coarse_second, coarse_primary, coarse_secondary, coarse_counts = h1_coarse_second_fraction(g)
        if coarse_second < second_min:
            continue
        gate_counts["H1_coarse_second_pass"] += 1

        compositions = g[FRACTIONS].to_numpy(float)
        hellinger = np.sqrt(compositions)
        labels, counts, separation_ratio = deterministic_two_means(hellinger)

        minor_fraction = float(counts.min() / n)
        if minor_fraction < second_min:
            continue
        gate_counts["continuous_minor_fraction_pass"] += 1

        # Order is used only to make exported tables deterministic. The primary
        # test is axial/sign-invariant, so swapping the two modes changes no H2 statistic.
        if counts[0] > counts[1]:
            major_label, minor_label = 0, 1
        elif counts[1] > counts[0]:
            major_label, minor_label = 1, 0
        else:
            c0 = compositions[labels == 0].mean(axis=0)
            c1 = compositions[labels == 1].mean(axis=0)
            major_label, minor_label = (0, 1) if tuple(c0.tolist()) <= tuple(c1.tolist()) else (1, 0)

        mu1 = compositions[labels == major_label].mean(axis=0)
        mu2 = compositions[labels == minor_label].mean(axis=0)
        delta = mu2 - mu1
        if not np.isclose(float(delta.sum()), 0.0, atol=1e-10):
            raise RuntimeError(f"Delta does not lie in zero-sum subspace for {sp}")
        norm = float(np.linalg.norm(delta))
        if not np.isfinite(norm) or norm <= EPS:
            continue
        gate_counts["nonzero_delta_pass"] += 1

        nmi = normalized_mutual_information(labels, g["morph"].to_numpy(str))
        row: dict[str, object] = {
            "species": sp,
            "n_classifiable": n,
            "H1_coarse_primary_morph": coarse_primary,
            "H1_coarse_secondary_morph": coarse_secondary,
            "H1_coarse_second_fraction": coarse_second,
            "continuous_major_n": int(counts[major_label]),
            "continuous_minor_n": int(counts[minor_label]),
            "continuous_minor_fraction": minor_fraction,
            "continuous_hellinger_separation_ratio": separation_ratio,
            "continuous_vs_coarse_nmi": nmi,
            "delta_norm": norm,
        }
        for m in MORPHS:
            row[f"coarse_n_{m}"] = int(coarse_counts[m])
        for j, colour in enumerate(BIO):
            row[f"mu1_{colour}"] = float(mu1[j])
            row[f"mu2_{colour}"] = float(mu2[j])
            row[f"delta_{colour}"] = float(delta[j])
        rows.append(row)
        vectors.append(delta)

    table = pd.DataFrame(rows)
    if len(table):
        table = table.sort_values("species").reset_index(drop=True)
    if not vectors:
        return table, np.empty((0, len(BIO)), dtype=float), gate_counts
    return table, np.vstack(vectors), gate_counts


def zero_sum_basis() -> np.ndarray:
    basis = linalg.null_space(np.ones((1, len(BIO)), dtype=float))
    expected = (len(BIO), len(BIO) - 1)
    if basis.shape != expected:
        raise RuntimeError(f"unexpected zero-sum basis shape: {basis.shape} != {expected}")
    if not np.allclose(basis.T @ basis, np.eye(len(BIO) - 1), atol=1e-10):
        raise RuntimeError("zero-sum basis is not orthonormal")
    return basis


def unit_tangent(vectors9: np.ndarray, basis: np.ndarray) -> np.ndarray:
    z = np.asarray(vectors9, float) @ basis
    norms = np.linalg.norm(z, axis=1)
    if np.any(~np.isfinite(norms)) or np.any(norms <= EPS):
        raise RuntimeError("non-finite or zero tangent vector")
    return z / norms[:, None]


def scatter_axis(u: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    if len(u) < 2:
        raise RuntimeError("too few vectors for scatter geometry")
    matrix = np.einsum("ni,nj->ij", u, u, optimize=True) / len(u)
    vals, vecs = np.linalg.eigh(matrix)
    order = np.argsort(vals)[::-1]
    vals = vals[order]
    vecs = vecs[:, order]
    return float(vals[0]), vecs[:, 0], vals


def isotropic_lambda_null(n_species: int, dim: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    out = np.empty(N_NULL, dtype=float)
    chunk = 250
    for start in range(0, N_NULL, chunk):
        end = min(N_NULL, start + chunk)
        z = rng.normal(size=(end - start, n_species, dim))
        z /= np.linalg.norm(z, axis=2, keepdims=True)
        matrix = np.einsum("bni,bnj->bij", z, z, optimize=True) / n_species
        out[start:end] = np.linalg.eigvalsh(matrix)[:, -1]
    return out


def fixed_axis_projection_null(n_species: int, dim: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    out = np.empty(N_NULL, dtype=float)
    chunk = 500
    for start in range(0, N_NULL, chunk):
        end = min(N_NULL, start + chunk)
        z = rng.normal(size=(end - start, n_species, dim))
        z /= np.linalg.norm(z, axis=2, keepdims=True)
        out[start:end] = np.mean(z[:, :, 0] ** 2, axis=1)
    return out


def mc_upper(observed: float, null: np.ndarray) -> float:
    return float((1 + np.sum(null >= observed)) / (len(null) + 1))


def signed_resultant(u: np.ndarray) -> float:
    return float(np.linalg.norm(np.mean(u, axis=0)))


def finite_summary(series: pd.Series) -> dict[str, float | int | None]:
    x = pd.to_numeric(series, errors="coerce").to_numpy(float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return {"n": 0, "median": None, "q25": None, "q75": None}
    return {
        "n": int(len(x)),
        "median": float(np.median(x)),
        "q25": float(np.quantile(x, 0.25)),
        "q75": float(np.quantile(x, 0.75)),
    }


def run_threshold(
    discovery: pd.DataFrame,
    reserve: pd.DataFrame,
    basis: np.ndarray,
    label: str,
    threshold: float,
    seed_offset: int,
) -> tuple[dict, pd.DataFrame, pd.DataFrame, np.ndarray, np.ndarray]:
    dtab, dv, dgates = species_delta_vectors(discovery, threshold)
    rtab, rv, rgates = species_delta_vectors(reserve, threshold)
    if len(dtab) < 20 or len(rtab) < 20:
        raise RuntimeError(
            f"insufficient label-free H2 vectors at {label}: discovery={len(dtab)}, reserve={len(rtab)}"
        )

    du = unit_tangent(dv, basis)
    ru = unit_tangent(rv, basis)
    dim = du.shape[1]

    d_lambda, d_axis, d_eigs = scatter_axis(du)
    d_null = isotropic_lambda_null(len(du), dim, SEED + seed_offset)
    d_p = mc_upper(d_lambda, d_null)

    r_lambda, r_axis, r_eigs = scatter_axis(ru)
    r_null = isotropic_lambda_null(len(ru), dim, SEED + seed_offset + 1)
    r_p = mc_upper(r_lambda, r_null)

    transport = float(np.mean((ru @ d_axis) ** 2))
    transport_null = fixed_axis_projection_null(len(ru), dim, SEED + seed_offset + 2)
    transport_p = mc_upper(transport, transport_null)
    axis_alignment = float(abs(np.dot(d_axis, r_axis)))

    palette_axis = basis @ d_axis
    if not np.isclose(float(palette_axis.sum()), 0.0, atol=1e-10):
        raise RuntimeError("mapped discovery axis is not zero-sum")

    discovery_pass = bool(d_p < 0.05)
    reserve_independent_pass = bool(r_p < 0.05)
    transport_pass = bool(transport_p < 0.05)
    claim_pass = bool(discovery_pass and transport_pass)

    result = {
        "coarse_H1_second_fraction_min": threshold,
        "continuous_H2_minor_cluster_fraction_min": threshold,
        "H2_mode_construction": "deterministic_two_means_on_Hellinger_transformed_9_colour_palette_without_morph_labels",
        "discovery": {
            "gate_counts": dgates,
            "H2_species": int(len(dtab)),
            "lambda1": d_lambda,
            "isotropic_upper_p": d_p,
            "eigenvalues": [float(x) for x in d_eigs],
            "signed_mean_resultant": signed_resultant(du),
            "minor_fraction_summary": finite_summary(dtab["continuous_minor_fraction"]),
            "separation_summary": finite_summary(dtab["continuous_hellinger_separation_ratio"]),
            "continuous_vs_coarse_nmi_summary": finite_summary(dtab["continuous_vs_coarse_nmi"]),
        },
        "reserve": {
            "gate_counts": rgates,
            "H2_species": int(len(rtab)),
            "lambda1": r_lambda,
            "isotropic_upper_p": r_p,
            "eigenvalues": [float(x) for x in r_eigs],
            "signed_mean_resultant": signed_resultant(ru),
            "minor_fraction_summary": finite_summary(rtab["continuous_minor_fraction"]),
            "separation_summary": finite_summary(rtab["continuous_hellinger_separation_ratio"]),
            "continuous_vs_coarse_nmi_summary": finite_summary(rtab["continuous_vs_coarse_nmi"]),
        },
        "transport": {
            "reserve_mean_squared_projection_on_frozen_discovery_axis": transport,
            "isotropic_upper_p": transport_p,
            "absolute_discovery_reserve_leading_axis_dot_diagnostic": axis_alignment,
            "discovery_axis_palette_loadings": {BIO[i]: float(palette_axis[i]) for i in range(len(BIO))},
        },
        "decision": {
            "discovery_anisotropy_pass": discovery_pass,
            "reserve_independent_anisotropy_diagnostic_pass": reserve_independent_pass,
            "reserve_transport_on_frozen_discovery_axis_pass": transport_pass,
            "H2_transportable_direction_pass": claim_pass,
        },
    }
    return result, dtab, rtab, d_null, transport_null


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    discovery = load_classifiable(DISCOVERY)
    reserve = load_classifiable(RESERVE)

    d_n40 = nclass40_count(discovery)
    r_n40 = nclass40_count(reserve)
    if d_n40 != 369:
        raise RuntimeError(f"discovery nclass>=40 fingerprint drift: {d_n40} != 369")
    if r_n40 != 363:
        raise RuntimeError(f"reserve nclass>=40 fingerprint drift: {r_n40} != 363")

    basis = zero_sum_basis()
    threshold_results: dict[str, dict] = {}
    for j, (label, threshold) in enumerate(THRESHOLDS.items()):
        result, dtab, rtab, d_null, transport_null = run_threshold(
            discovery, reserve, basis, label, threshold, 10 * j
        )
        threshold_results[label] = result
        dtab.to_csv(OUT / f"{label}_discovery_delta_vectors.csv", index=False)
        rtab.to_csv(OUT / f"{label}_reserve_delta_vectors.csv", index=False)
        pd.DataFrame({"lambda1": d_null}).to_csv(
            OUT / f"{label}_discovery_isotropic_lambda_null.csv", index=False
        )
        pd.DataFrame({"mean_squared_projection": transport_null}).to_csv(
            OUT / f"{label}_reserve_transport_null.csv", index=False
        )

    primary = threshold_results["primary_0_10"]
    strict = threshold_results["strict_0_20"]
    final_pass = bool(
        primary["decision"]["H2_transportable_direction_pass"]
        and strict["decision"]["H2_transportable_direction_pass"]
    )
    if final_pass:
        verdict = "H2_LABEL_FREE_DIRECTION_TRANSPORTS_AND_SURVIVES_STRICT_THRESHOLD"
    elif primary["decision"]["H2_transportable_direction_pass"]:
        verdict = "H2_LABEL_FREE_PRIMARY_ONLY_STRICT_SENSITIVITY_FAILS"
    elif primary["decision"]["discovery_anisotropy_pass"]:
        verdict = "H2_LABEL_FREE_DISCOVERY_ONLY_NO_RESERVE_TRANSPORT"
    else:
        verdict = "H2_LABEL_FREE_NO_RECURRENT_DELTA_DIRECTION"

    result = {
        "analysis": "polymorphism_delta_geometry_validation",
        "date_jst": "2026-09-12",
        "protocol": "docs/POLYMORPHISM_42111_H1_H2_ELIGIBILITY_PROTOCOL_20260912.md",
        "new_image_acquisition": False,
        "mixed_uncertain_used": False,
        "four_state_labels_used_for_H2_geometry": False,
        "four_state_labels_used_for_H1_admission_only": True,
        "biological_palette_dimensions": BIO,
        "H2_continuous_space": "Hellinger_transform_of_normalized_9_colour_flower_palette",
        "zero_sum_dimension_for_delta": int(basis.shape[1]),
        "minimum_classifiable": MIN_CLASSIFIABLE,
        "null_replicates": N_NULL,
        "discovery_species_nclass_ge40": d_n40,
        "reserve_species_nclass_ge40": r_n40,
        "thresholds": threshold_results,
        "decision": {
            "primary_and_strict_transport_pass": final_pass,
            "verdict": verdict,
        },
        "claim_boundary": (
            "Validation-cohort test of recurrent species-specific continuous-palette colour-axis geometry, "
            "conditional on the frozen coarse H1 admission gate. H2 modes are estimated without four-state "
            "labels. The result does not estimate 42,111-species polymorphism prevalence and does not open "
            "H3 ecology/phylogeny. A surviving H2 requires a stronger structure-preserving null before "
            "ecological interpretation."
        ),
    }
    (OUT / "result.json").write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8"
    )

    lines = [
        "# Label-free species-specific Delta geometry validation",
        "",
        f"**Verdict: `{verdict}`**",
        "",
        f"- discovery nclass>=40 frame: **{d_n40} species**",
        f"- reserve nclass>=40 frame: **{r_n40} species**",
        "- four-state labels used to construct H2 modes: **false**",
        "- mixed/unresolved images used as biological modes: **false**",
        "",
    ]
    for label in THRESHOLDS:
        x = threshold_results[label]
        lines += [
            f"## {label}",
            "",
            f"- discovery H2 N: **{x['discovery']['H2_species']}**, lambda1 **{x['discovery']['lambda1']:.6f}**, p **{x['discovery']['isotropic_upper_p']:.6g}**",
            f"- reserve H2 N: **{x['reserve']['H2_species']}**, lambda1 **{x['reserve']['lambda1']:.6f}**, p **{x['reserve']['isotropic_upper_p']:.6g}**",
            f"- reserve projection on frozen discovery axis: **{x['transport']['reserve_mean_squared_projection_on_frozen_discovery_axis']:.6f}**, p **{x['transport']['isotropic_upper_p']:.6g}**",
            f"- absolute discovery/reserve axis alignment (diagnostic): **{x['transport']['absolute_discovery_reserve_leading_axis_dot_diagnostic']:.6f}**",
            "",
        ]
    (OUT / "RESULT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
