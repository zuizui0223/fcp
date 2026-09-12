#!/usr/bin/env python3
from __future__ import annotations

import json
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
        raise RuntimeError(f"classifiable rows with invalid biological palette in {path}: {int((~good).sum())}")
    work.loc[:, FRACTIONS] = x
    return work


def nclass40_count(work: pd.DataFrame) -> int:
    sizes = work.groupby("species", sort=False).size()
    return int((sizes >= MIN_CLASSIFIABLE).sum())


def species_delta_vectors(work: pd.DataFrame, second_min: float) -> tuple[pd.DataFrame, np.ndarray]:
    rows: list[dict[str, object]] = []
    vectors: list[np.ndarray] = []

    for sp, g in work.groupby("species", sort=True):
        n = int(len(g))
        if n < MIN_CLASSIFIABLE:
            continue

        cnt = Counter(g["morph"].astype(str))
        order = sorted(MORPHS, key=lambda m: (-cnt.get(m, 0), m))
        primary, secondary = order[:2]
        primary_n = int(cnt.get(primary, 0))
        secondary_n = int(cnt.get(secondary, 0))
        third_n = int(cnt.get(order[2], 0))
        second_fraction = float(secondary_n / n)

        if second_fraction < second_min:
            continue

        # A tie between second and third states makes the top-two pair itself
        # underidentified. A tie between first and second only flips Delta's sign,
        # which does not affect the sign-invariant primary statistic.
        second_third_tie = bool(secondary_n == third_n)
        if second_third_tie:
            continue

        x1 = g.loc[g["morph"].eq(primary), FRACTIONS].to_numpy(float)
        x2 = g.loc[g["morph"].eq(secondary), FRACTIONS].to_numpy(float)
        if len(x1) != primary_n or len(x2) != secondary_n:
            raise RuntimeError("state count mismatch")
        mu1 = x1.mean(axis=0)
        mu2 = x2.mean(axis=0)
        delta = mu2 - mu1
        if not np.isclose(float(delta.sum()), 0.0, atol=1e-10):
            raise RuntimeError(f"Delta does not lie in zero-sum subspace for {sp}")
        norm = float(np.linalg.norm(delta))
        if not np.isfinite(norm) or norm <= EPS:
            continue

        row: dict[str, object] = {
            "species": sp,
            "n_classifiable": n,
            "primary_morph": primary,
            "secondary_morph": secondary,
            "primary_n": primary_n,
            "secondary_n": secondary_n,
            "second_fraction": second_fraction,
            "delta_norm": norm,
        }
        for j, colour in enumerate(BIO):
            row[f"mu1_{colour}"] = float(mu1[j])
            row[f"mu2_{colour}"] = float(mu2[j])
            row[f"delta_{colour}"] = float(delta[j])
        rows.append(row)
        vectors.append(delta)

    table = pd.DataFrame(rows).sort_values("species").reset_index(drop=True)
    if not vectors:
        return table, np.empty((0, len(BIO)), dtype=float)
    return table, np.vstack(vectors)


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
    M = np.einsum("ni,nj->ij", u, u, optimize=True) / len(u)
    vals, vecs = np.linalg.eigh(M)
    order = np.argsort(vals)[::-1]
    vals = vals[order]
    vecs = vecs[:, order]
    axis = vecs[:, 0]
    return float(vals[0]), axis, vals


def isotropic_lambda_null(n_species: int, dim: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    out = np.empty(N_NULL, dtype=float)
    chunk = 250
    for start in range(0, N_NULL, chunk):
        end = min(N_NULL, start + chunk)
        z = rng.normal(size=(end - start, n_species, dim))
        z /= np.linalg.norm(z, axis=2, keepdims=True)
        M = np.einsum("bni,bnj->bij", z, z, optimize=True) / n_species
        out[start:end] = np.linalg.eigvalsh(M)[:, -1]
    return out


def fixed_axis_projection_null(n_species: int, dim: int, seed: int) -> np.ndarray:
    # For a fixed unit axis and independent isotropic unit vectors, the axis can
    # be represented as coordinate 0 without loss of generality.
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


def pair_counts(table: pd.DataFrame) -> dict[str, int]:
    if table.empty:
        return {}
    pair = table.apply(
        lambda r: "__".join(sorted([str(r["primary_morph"]), str(r["secondary_morph"])])),
        axis=1,
    )
    vc = pair.value_counts().sort_index()
    return {str(k): int(v) for k, v in vc.items()}


def signed_resultant(u: np.ndarray) -> float:
    return float(np.linalg.norm(np.mean(u, axis=0)))


def run_threshold(
    discovery: pd.DataFrame,
    reserve: pd.DataFrame,
    basis: np.ndarray,
    label: str,
    threshold: float,
    seed_offset: int,
) -> tuple[dict, pd.DataFrame, pd.DataFrame, np.ndarray, np.ndarray]:
    dtab, dv = species_delta_vectors(discovery, threshold)
    rtab, rv = species_delta_vectors(reserve, threshold)
    if len(dtab) < 20 or len(rtab) < 20:
        raise RuntimeError(f"insufficient H2 vectors at {label}: discovery={len(dtab)}, reserve={len(rtab)}")

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
        "second_fraction_min": threshold,
        "discovery": {
            "H2_species": int(len(dtab)),
            "lambda1": d_lambda,
            "isotropic_upper_p": d_p,
            "eigenvalues": [float(x) for x in d_eigs],
            "signed_mean_resultant": signed_resultant(du),
            "pair_counts": pair_counts(dtab),
        },
        "reserve": {
            "H2_species": int(len(rtab)),
            "lambda1": r_lambda,
            "isotropic_upper_p": r_p,
            "eigenvalues": [float(x) for x in r_eigs],
            "signed_mean_resultant": signed_resultant(ru),
            "pair_counts": pair_counts(rtab),
        },
        "transport": {
            "reserve_mean_squared_projection_on_discovery_axis": transport,
            "isotropic_upper_p": transport_p,
            "absolute_discovery_reserve_leading_axis_dot": axis_alignment,
            "discovery_axis_palette_loadings": {BIO[i]: float(palette_axis[i]) for i in range(len(BIO))},
        },
        "decision": {
            "discovery_anisotropy_pass": discovery_pass,
            "reserve_independent_anisotropy_pass": reserve_independent_pass,
            "reserve_transport_on_discovery_axis_pass": transport_pass,
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
        pd.DataFrame({"lambda1": d_null}).to_csv(OUT / f"{label}_discovery_isotropic_lambda_null.csv", index=False)
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
        verdict = "H2_DIRECTION_TRANSPORTS_AND_SURVIVES_STRICT_THRESHOLD"
    elif primary["decision"]["H2_transportable_direction_pass"]:
        verdict = "H2_PRIMARY_ONLY_STRICT_SENSITIVITY_FAILS"
    elif primary["decision"]["discovery_anisotropy_pass"]:
        verdict = "H2_DISCOVERY_ONLY_NO_RESERVE_TRANSPORT"
    else:
        verdict = "H2_NO_RECURRENT_DELTA_DIRECTION"

    result = {
        "analysis": "polymorphism_delta_geometry_validation",
        "date_jst": "2026-09-12",
        "protocol": "docs/POLYMORPHISM_42111_H1_H2_ELIGIBILITY_PROTOCOL_20260912.md",
        "new_image_acquisition": False,
        "mixed_uncertain_used": False,
        "biological_palette_dimensions": BIO,
        "zero_sum_dimension": int(basis.shape[1]),
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
            "Validation-cohort test of recurrent species-specific primary/secondary colour-axis geometry only. "
            "It does not estimate 42,111-species polymorphism prevalence and does not open H3 ecology/phylogeny."
        ),
    }
    (OUT / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    lines = [
        "# Species-specific Delta geometry validation",
        "",
        f"**Verdict: `{verdict}`**",
        "",
        f"- discovery nclass>=40 frame: **{d_n40} species**",
        f"- reserve nclass>=40 frame: **{r_n40} species**",
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
            f"- reserve projection on discovery axis: **{x['transport']['reserve_mean_squared_projection_on_discovery_axis']:.6f}**, p **{x['transport']['isotropic_upper_p']:.6g}**",
            f"- absolute discovery/reserve axis alignment: **{x['transport']['absolute_discovery_reserve_leading_axis_dot']:.6f}**",
            "",
        ]
    (OUT / "RESULT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
