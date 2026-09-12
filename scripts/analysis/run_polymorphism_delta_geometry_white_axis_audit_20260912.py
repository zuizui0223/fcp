#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "analysis"))

from run_polymorphism_delta_geometry_validation_20260912 import (  # noqa: E402
    BIO,
    THRESHOLDS,
    load_classifiable,
)
from run_polymorphism_delta_geometry_structured_null_20260912 import (  # noqa: E402
    prepare_selected,
    permute_within_coarse_morph,
    delta_units_from_compositions,
)

DISCOVERY = ROOT / "data" / "derived" / "global_monte_carlo_measured_photos_v1.csv"
RESERVE = ROOT / "data" / "derived" / "rgfca_reserve_replication_measured_photos_v1.csv"
H2_OUT = ROOT / "results" / "polymorphism_delta_geometry_validation_20260912"
OUT = ROOT / "results" / "polymorphism_delta_geometry_white_axis_audit_20260912"
N_ISO = 10_000
N_STRUCTURED = 999
SEED = 20260912
EPS = 1e-12


def unit_from_table(table: pd.DataFrame) -> np.ndarray:
    x = table[[f"delta_{c}" for c in BIO]].to_numpy(float)
    n = np.linalg.norm(x, axis=1)
    if np.any(~np.isfinite(n)) or np.any(n <= EPS):
        raise RuntimeError("non-finite or zero observed Delta")
    return x / n[:, None]


def white_contrast() -> np.ndarray:
    # Canonical zero-sum achromatic contrast: white versus the equal mean of all
    # eight non-white palette coordinates. Defined a priori from coordinate labels.
    q = np.array([1.0] + [-1.0 / 8.0] * 8, dtype=float)
    q /= np.linalg.norm(q)
    if not np.isclose(float(q.sum()), 0.0, atol=1e-12):
        raise RuntimeError("white contrast is not zero-sum")
    return q


def residualize_white(units9: np.ndarray, q: np.ndarray) -> np.ndarray:
    r = units9 - (units9 @ q)[:, None] * q[None, :]
    n = np.linalg.norm(r, axis=1)
    if np.any(~np.isfinite(n)) or np.any(n <= EPS):
        raise RuntimeError("white-axis residualization produced a degenerate vector")
    return r / n[:, None]


def leading_axis(units9: np.ndarray) -> tuple[float, np.ndarray]:
    matrix = units9.T @ units9 / len(units9)
    vals, vecs = np.linalg.eigh(matrix)
    j = int(np.argmax(vals))
    return float(vals[j]), vecs[:, j]


def mean_squared_projection(units9: np.ndarray, axis9: np.ndarray) -> float:
    axis9 = axis9 / np.linalg.norm(axis9)
    return float(np.mean((units9 @ axis9) ** 2))


def isotropic_lambda_null(n_species: int, dim: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    out = np.empty(N_ISO, dtype=float)
    chunk = 250
    for start in range(0, N_ISO, chunk):
        end = min(N_ISO, start + chunk)
        z = rng.normal(size=(end - start, n_species, dim))
        z /= np.linalg.norm(z, axis=2, keepdims=True)
        matrix = np.einsum("bni,bnj->bij", z, z, optimize=True) / n_species
        out[start:end] = np.linalg.eigvalsh(matrix)[:, -1]
    return out


def fixed_axis_projection_null(n_species: int, dim: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    out = np.empty(N_ISO, dtype=float)
    chunk = 500
    for start in range(0, N_ISO, chunk):
        end = min(N_ISO, start + chunk)
        z = rng.normal(size=(end - start, n_species, dim))
        z /= np.linalg.norm(z, axis=2, keepdims=True)
        out[start:end] = np.mean(z[:, :, 0] ** 2, axis=1)
    return out


def mc_upper(observed: float, null: np.ndarray) -> float:
    return float((1 + np.sum(null >= observed)) / (len(null) + 1))


def qsummary(x: np.ndarray) -> dict[str, float]:
    return {
        "mean": float(np.mean(x)),
        "sd": float(np.std(x, ddof=1)),
        "q025": float(np.quantile(x, 0.025)),
        "q50": float(np.quantile(x, 0.5)),
        "q975": float(np.quantile(x, 0.975)),
    }


def orient_axis(axis: np.ndarray) -> np.ndarray:
    axis = axis.copy()
    j = int(np.argmax(np.abs(axis)))
    if axis[j] < 0:
        axis *= -1
    return axis


def white_pair_mask(table: pd.DataFrame) -> np.ndarray:
    return (
        table["H1_coarse_primary_morph"].eq("white")
        | table["H1_coarse_secondary_morph"].eq("white")
    ).to_numpy(bool)


def structured_residual_null(
    dprep: dict,
    rprep: dict,
    q: np.ndarray,
    observed_axis9: np.ndarray,
    seed_offset: int,
) -> tuple[np.ndarray, np.ndarray]:
    rng_d = np.random.default_rng(SEED + seed_offset)
    rng_r = np.random.default_rng(SEED + seed_offset + 1)
    null_lambda = np.empty(N_STRUCTURED, dtype=float)
    null_transport = np.empty(N_STRUCTURED, dtype=float)
    for b in range(N_STRUCTURED):
        dperm = permute_within_coarse_morph(dprep, rng_d)
        du = delta_units_from_compositions(dperm, dprep["species_indices"])
        dur = residualize_white(du, q)
        null_lambda[b] = leading_axis(dur)[0]

        rperm = permute_within_coarse_morph(rprep, rng_r)
        ru = delta_units_from_compositions(rperm, rprep["species_indices"])
        rur = residualize_white(ru, q)
        null_transport[b] = mean_squared_projection(rur, observed_axis9)
    return null_lambda, null_transport


def structured_plain_null(
    dprep: dict,
    rprep: dict,
    observed_axis9: np.ndarray,
    seed_offset: int,
) -> tuple[np.ndarray, np.ndarray]:
    rng_d = np.random.default_rng(SEED + seed_offset)
    rng_r = np.random.default_rng(SEED + seed_offset + 1)
    null_lambda = np.empty(N_STRUCTURED, dtype=float)
    null_transport = np.empty(N_STRUCTURED, dtype=float)
    for b in range(N_STRUCTURED):
        dperm = permute_within_coarse_morph(dprep, rng_d)
        du = delta_units_from_compositions(dperm, dprep["species_indices"])
        null_lambda[b] = leading_axis(du)[0]

        rperm = permute_within_coarse_morph(rprep, rng_r)
        ru = delta_units_from_compositions(rperm, rprep["species_indices"])
        null_transport[b] = mean_squared_projection(ru, observed_axis9)
    return null_lambda, null_transport


def analyse_threshold(
    label: str,
    discovery: pd.DataFrame,
    reserve: pd.DataFrame,
    q: np.ndarray,
    seed_offset: int,
) -> tuple[dict, dict[str, np.ndarray]]:
    dtab = pd.read_csv(H2_OUT / f"{label}_discovery_delta_vectors.csv")
    rtab = pd.read_csv(H2_OUT / f"{label}_reserve_delta_vectors.csv")
    du = unit_from_table(dtab)
    ru = unit_from_table(rtab)

    # Test 1: remove the canonical white-vs-equal-nonwhite contrast from every
    # observed species axis, then ask whether recurrent geometry remains.
    dur = residualize_white(du, q)
    rur = residualize_white(ru, q)
    resid_lambda, resid_axis = leading_axis(dur)
    resid_axis = orient_axis(resid_axis)
    resid_transport = mean_squared_projection(rur, resid_axis)
    iso_lambda = isotropic_lambda_null(len(dur), 7, SEED + seed_offset)
    iso_transport = fixed_axis_projection_null(len(rur), 7, SEED + seed_offset + 1)

    dprep = prepare_selected(discovery, dtab["species"].astype(str).tolist())
    rprep = prepare_selected(reserve, rtab["species"].astype(str).tolist())
    structured_lambda, structured_transport = structured_residual_null(
        dprep, rprep, q, resid_axis, seed_offset + 10
    )

    residual = {
        "definition": (
            "Project each unit Delta_i onto the orthogonal complement of the a-priori "
            "white-versus-equal-mean-of-eight-nonwhite zero-sum contrast, renormalize, "
            "then repeat discovery concentration and frozen-axis reserve transport."
        ),
        "residual_dimension": 7,
        "discovery_species": int(len(dtab)),
        "reserve_species": int(len(rtab)),
        "discovery_lambda1": resid_lambda,
        "isotropic_upper_p": mc_upper(resid_lambda, iso_lambda),
        "structured_null_upper_p": mc_upper(resid_lambda, structured_lambda),
        "reserve_mean_squared_projection_on_frozen_residual_axis": resid_transport,
        "reserve_isotropic_upper_p": mc_upper(resid_transport, iso_transport),
        "reserve_structured_null_upper_p": mc_upper(resid_transport, structured_transport),
        "discovery_residual_axis_palette_loadings": {
            BIO[i]: float(resid_axis[i]) for i in range(len(BIO))
        },
        "structured_lambda_null_summary": qsummary(structured_lambda),
        "structured_transport_null_summary": qsummary(structured_transport),
    }
    residual["pass"] = bool(
        residual["structured_null_upper_p"] < 0.05
        and residual["reserve_structured_null_upper_p"] < 0.05
    )

    # Test 2: restrict to species whose two leading coarse H1 states do not
    # include white. This is a diagnostic because strict-threshold N is small.
    d_nonwhite = ~white_pair_mask(dtab)
    r_nonwhite = ~white_pair_mask(rtab)
    dnw_tab = dtab.loc[d_nonwhite].reset_index(drop=True)
    rnw_tab = rtab.loc[r_nonwhite].reset_index(drop=True)
    dnw = du[d_nonwhite]
    rnw = ru[r_nonwhite]
    nonwhite: dict[str, object] = {
        "definition": "Species whose H1 primary and secondary coarse states both exclude white.",
        "discovery_species": int(len(dnw)),
        "reserve_species": int(len(rnw)),
        "role": "diagnostic_not_primary_gate",
    }
    nulls: dict[str, np.ndarray] = {
        "residual_iso_lambda": iso_lambda,
        "residual_iso_transport": iso_transport,
        "residual_structured_lambda": structured_lambda,
        "residual_structured_transport": structured_transport,
    }

    if len(dnw) >= 8 and len(rnw) >= 8:
        nw_lambda, nw_axis = leading_axis(dnw)
        nw_axis = orient_axis(nw_axis)
        nw_transport = mean_squared_projection(rnw, nw_axis)
        nw_iso_lambda = isotropic_lambda_null(len(dnw), 8, SEED + seed_offset + 20)
        nw_iso_transport = fixed_axis_projection_null(len(rnw), 8, SEED + seed_offset + 21)

        dnprep = prepare_selected(discovery, dnw_tab["species"].astype(str).tolist())
        rnprep = prepare_selected(reserve, rnw_tab["species"].astype(str).tolist())
        nw_struct_lambda, nw_struct_transport = structured_plain_null(
            dnprep, rnprep, nw_axis, seed_offset + 30
        )
        nonwhite.update(
            {
                "discovery_lambda1": nw_lambda,
                "isotropic_upper_p": mc_upper(nw_lambda, nw_iso_lambda),
                "structured_null_upper_p": mc_upper(nw_lambda, nw_struct_lambda),
                "reserve_mean_squared_projection_on_frozen_axis": nw_transport,
                "reserve_isotropic_upper_p": mc_upper(nw_transport, nw_iso_transport),
                "reserve_structured_null_upper_p": mc_upper(nw_transport, nw_struct_transport),
                "discovery_axis_palette_loadings": {
                    BIO[i]: float(nw_axis[i]) for i in range(len(BIO))
                },
                "structured_lambda_null_summary": qsummary(nw_struct_lambda),
                "structured_transport_null_summary": qsummary(nw_struct_transport),
                "pass": bool(
                    mc_upper(nw_lambda, nw_struct_lambda) < 0.05
                    and mc_upper(nw_transport, nw_struct_transport) < 0.05
                ),
            }
        )
        nulls.update(
            {
                "nonwhite_iso_lambda": nw_iso_lambda,
                "nonwhite_iso_transport": nw_iso_transport,
                "nonwhite_structured_lambda": nw_struct_lambda,
                "nonwhite_structured_transport": nw_struct_transport,
            }
        )
    else:
        nonwhite["pass"] = None
        nonwhite["reason"] = "too_few_species_for_prespecified_diagnostic_minimum_8_per_cohort"

    result = {
        "white_pair_fraction_discovery": float(white_pair_mask(dtab).mean()),
        "white_pair_fraction_reserve": float(white_pair_mask(rtab).mean()),
        "white_contrast_palette_loadings": {BIO[i]: float(q[i]) for i in range(len(BIO))},
        "white_contrast_removed": residual,
        "nonwhite_top_two_subset": nonwhite,
    }
    return result, nulls


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    discovery = load_classifiable(DISCOVERY)
    reserve = load_classifiable(RESERVE)
    q = white_contrast()

    results: dict[str, dict] = {}
    for j, label in enumerate(THRESHOLDS):
        result, nulls = analyse_threshold(label, discovery, reserve, q, 100 * j)
        results[label] = result
        for name, values in nulls.items():
            pd.DataFrame({name: values}).to_csv(OUT / f"{label}_{name}.csv", index=False)

    primary_resid = bool(results["primary_0_10"]["white_contrast_removed"]["pass"])
    strict_resid = bool(results["strict_0_20"]["white_contrast_removed"]["pass"])
    primary_nonwhite = results["primary_0_10"]["nonwhite_top_two_subset"].get("pass")

    if primary_resid and strict_resid and primary_nonwhite is True:
        verdict = "H2_EXTENDS_BEYOND_WHITE_AXIS_AND_REPLICATES_WITHOUT_WHITE_TOP_PAIR"
    elif primary_resid and strict_resid:
        verdict = "H2_EXTENDS_BEYOND_WHITE_AXIS_NONWHITE_SUBSET_UNRESOLVED"
    elif primary_resid:
        verdict = "H2_WHITE_RESIDUAL_PRIMARY_ONLY_STRICT_FAILS"
    else:
        verdict = "H2_DOMINATED_BY_WHITE_VERSUS_NONWHITE_AXIS"

    output = {
        "analysis": "polymorphism_delta_geometry_white_axis_audit",
        "date_jst": "2026-09-12",
        "isotropic_null_replicates": N_ISO,
        "structured_null_replicates": N_STRUCTURED,
        "thresholds": results,
        "decision": {
            "primary_white_residual_pass": primary_resid,
            "strict_white_residual_pass": strict_resid,
            "primary_nonwhite_subset_pass": primary_nonwhite,
            "verdict": verdict,
        },
        "claim_boundary": (
            "This audit asks whether the validated H2 recurrent geometry is exhausted by a canonical "
            "white-versus-nonwhite compositional contrast. Passing the residual and nonwhite diagnostics "
            "supports additional recurrent hue geometry, but does not identify ecological or genetic mechanism."
        ),
    }
    (OUT / "result.json").write_text(
        json.dumps(output, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8"
    )

    lines = [
        "# H2 white-axis audit",
        "",
        f"**Verdict: `{verdict}`**",
        "",
    ]
    for label in THRESHOLDS:
        x = results[label]
        r = x["white_contrast_removed"]
        nw = x["nonwhite_top_two_subset"]
        lines += [
            f"## {label}",
            "",
            f"- white-pair fraction discovery / reserve: **{x['white_pair_fraction_discovery']:.3f} / {x['white_pair_fraction_reserve']:.3f}**",
            f"- after removing canonical white contrast, discovery lambda1: **{r['discovery_lambda1']:.6f}**",
            f"- residual structured-null p: **{r['structured_null_upper_p']:.6g}**",
            f"- residual reserve projection: **{r['reserve_mean_squared_projection_on_frozen_residual_axis']:.6f}**",
            f"- residual structured-null transport p: **{r['reserve_structured_null_upper_p']:.6g}**",
            f"- nonwhite top-two subset N discovery / reserve: **{nw['discovery_species']} / {nw['reserve_species']}**",
            f"- nonwhite subset pass: **{nw.get('pass')}**",
            "",
        ]
    (OUT / "RESULT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
