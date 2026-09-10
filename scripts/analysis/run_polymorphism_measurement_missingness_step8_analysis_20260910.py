#!/usr/bin/env python3
"""Separate pure technical ROI/flip failures from ambiguity-related missingness and retest FCP main claims."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "polymorphism_measurement_missingness_step8_analysis_20260910"
OUT.mkdir(parents=True, exist_ok=True)

DISC_RAW = ROOT / "data" / "derived" / "global_monte_carlo_measured_photos_v1.csv"
RES_RAW = ROOT / "data" / "derived" / "rgfca_reserve_replication_measured_photos_v1.csv"
DISC_ATTR = ROOT / "results" / "polymorphism_species_attributes_step4_association_20260910" / "species_analysis_table.csv"
RES_ATTR = ROOT / "results" / "polymorphism_spatial_span_adjusted_step6_20260910" / "reserve_species_span_adjusted_input.csv"
DISC_SPATIAL = ROOT / "results" / "polymorphism_spatial_subset_step5_20260909" / "species_membership_and_observed_rho.csv"
RES_SPATIAL = ROOT / "results" / "polymorphism_spatial_reserve_step5b_20260909" / "reserve_species_membership_and_observed_metrics.csv"
PREFLIGHT = ROOT / "results" / "polymorphism_measurement_missingness_step8_preflight_20260910" / "result.json"
STEP7B = ROOT / "results" / "polymorphism_genus_adjustment_decomposition_step7b_20260910" / "result.json"

STATUSES = [
    "classified_four_state_morph",
    "not_evaluable_roi_or_flip_gate",
    "not_evaluable_ambiguous_palette_composition",
    "not_evaluable_no_biological_palette_mass",
]
STATUS_TO_RATE = {
    "classified_four_state_morph": "classifiable_rate",
    "not_evaluable_roi_or_flip_gate": "technical_failure_rate",
    "not_evaluable_ambiguous_palette_composition": "ambiguous_palette_rate",
    "not_evaluable_no_biological_palette_mass": "no_biological_palette_rate",
}
N_PERM = 20_000
SPATIAL_NULL = 999
SEED = 20260920
EXPECTED_DISC = 369
EXPECTED_RES = 363
RES_METRICS = ["primary", "observer_pair_exclusion", "calendar_quarter_stratification", "matched_background_differential"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--discovery-artifact-root", type=Path, required=True)
    p.add_argument("--reserve-artifact-root", type=Path, required=True)
    return p.parse_args()


def write_json(path: Path, obj: Any) -> None:
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def build_status_panel(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, usecols=["species", "inat_taxon_id", "measurement_status"], low_memory=False)
    if df["species"].isna().any() or df["measurement_status"].isna().any():
        raise RuntimeError(f"missing species/status in {path}")
    observed_statuses = sorted(df["measurement_status"].astype(str).unique().tolist())
    if observed_statuses != sorted(STATUSES):
        raise RuntimeError(f"measurement-status vocabulary drift in {path}: {observed_statuses}")
    size = df.groupby("species", sort=False).size()
    if not size.eq(100).all():
        raise RuntimeError(f"fixed 100-photo denominator failed in {path}: min={size.min()} max={size.max()}")
    rows: list[dict[str, Any]] = []
    for species, g in df.groupby("species", sort=True):
        counts = g["measurement_status"].astype(str).value_counts()
        row: dict[str, Any] = {
            "species": str(species),
            "inat_taxon_id_status": int(g["inat_taxon_id"].iloc[0]),
            "raw_photo_rows": int(len(g)),
        }
        total = len(g)
        for status, col in STATUS_TO_RATE.items():
            row[col] = float(counts.get(status, 0) / total)
            row[col.replace("_rate", "_count")] = int(counts.get(status, 0))
        rows.append(row)
    out = pd.DataFrame(rows)
    rates = out[list(STATUS_TO_RATE.values())].sum(axis=1)
    if not np.allclose(rates.to_numpy(float), 1.0, atol=1e-12):
        raise RuntimeError("measurement process classes do not sum to one")
    return out


def add_genus(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["genus"] = out["species"].astype(str).str.strip().str.split().str[0]
    return out


def repeated_group_pairs(labels: pd.Series) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    lab = labels.fillna("").astype(str).to_numpy()
    groups: dict[str, np.ndarray] = {}
    for name in sorted(set(lab)):
        if not name:
            continue
        idx = np.flatnonzero(lab == name)
        if len(idx) >= 2:
            groups[name] = idx
    species_idx = sorted({int(i) for idx in groups.values() for i in idx})
    remap = {old: new for new, old in enumerate(species_idx)}
    pi: list[int] = []
    pj: list[int] = []
    ww: list[float] = []
    ng = len(groups)
    if ng:
        for _name, idx in groups.items():
            npairs = len(idx) * (len(idx) - 1) // 2
            w = 1.0 / (ng * npairs)
            for a in range(len(idx) - 1):
                for b in range(a + 1, len(idx)):
                    pi.append(remap[int(idx[a])])
                    pj.append(remap[int(idx[b])])
                    ww.append(w)
    meta = {
        "repeated_groups": int(ng),
        "species_in_repeated_groups": int(len(species_idx)),
        "gate_pass": bool(ng >= 10 and len(species_idx) >= 30),
    }
    return np.asarray(species_idx, int), np.asarray(pi, int), np.asarray(pj, int), np.asarray(ww, float), meta


def clustering_test(values: np.ndarray, labels: pd.Series, *, seed: int) -> dict[str, Any]:
    idx, pi, pj, weights, meta = repeated_group_pairs(labels)
    res = dict(meta)
    if not res["gate_pass"]:
        res.update({"tested": False, "clustering_gain": None, "p_lower": 1.0})
        return res
    z = np.asarray(values, float)[idx]
    if not np.isfinite(z).all() or not np.isclose(weights.sum(), 1.0):
        raise RuntimeError("invalid genus clustering input")
    observed = float(np.sum(np.abs(z[pi] - z[pj]) * weights))
    rng = np.random.default_rng(seed)
    null = np.empty(N_PERM, float)
    for start in range(0, N_PERM, 500):
        stop = min(N_PERM, start + 500)
        perms = np.stack([rng.permutation(len(z)) for _ in range(stop - start)])
        zz = z[perms]
        null[start:stop] = np.sum(np.abs(zz[:, pi] - zz[:, pj]) * weights[None, :], axis=1)
    nm = float(null.mean())
    gain = float(1.0 - observed / nm) if nm > 0 else math.nan
    p = float((1 + np.sum(null <= observed + 1e-15)) / (N_PERM + 1))
    res.update({
        "tested": True,
        "W_observed": observed,
        "null_mean_W": nm,
        "null_q025": float(np.quantile(null, 0.025)),
        "null_q975": float(np.quantile(null, 0.975)),
        "clustering_gain": gain,
        "p_lower": p,
        "supported_at_0_05": bool(gain > 0 and p < 0.05),
    })
    return res


def rank_residual(values: np.ndarray, controls: list[np.ndarray]) -> np.ndarray:
    y = stats.rankdata(np.asarray(values, float)).astype(float)
    cols = [np.ones(len(y), float)]
    for v in controls:
        v = np.asarray(v, float)
        if not np.isfinite(v).all():
            raise RuntimeError("nonfinite residualization control")
        cols.append(stats.rankdata(v).astype(float))
    X = np.column_stack(cols)
    return y - X @ np.linalg.lstsq(X, y, rcond=None)[0]


def partial_rank(values_x: np.ndarray, values_y: np.ndarray, controls: list[np.ndarray]) -> float:
    ex = rank_residual(values_x, controls)
    ey = rank_residual(values_y, controls)
    return float(np.corrcoef(ex, ey)[0, 1])


def load_discovery_spatial_null(root: Path) -> tuple[pd.DataFrame, np.ndarray]:
    files = sorted(root.glob("global-rgfca-within-species-omnibus-shard-*/species_rho_permutations.csv"))
    if len(files) != 20:
        raise RuntimeError(f"expected 20 discovery spatial shard files, found {len(files)}")
    long = pd.concat([pd.read_csv(p, usecols=["species", "permutation_index", "rho"]) for p in files], ignore_index=True)
    if long.duplicated(["species", "permutation_index"]).any():
        raise RuntimeError("duplicate discovery spatial null row")
    wide = long.pivot(index="species", columns="permutation_index", values="rho").sort_index()
    if len(wide) != EXPECTED_DISC or -1 not in wide.columns or not all(i in wide.columns for i in range(SPATIAL_NULL)):
        raise RuntimeError(f"discovery spatial null shape drift: {wide.shape}")
    obs = wide[-1].rename("spatial_observed_rho").reset_index()
    return obs, wide[list(range(SPATIAL_NULL))].to_numpy(float)


def load_reserve_spatial_null(root: Path) -> tuple[pd.DataFrame, np.ndarray]:
    dirs = sorted([p for p in root.glob("reserve-inference-shard-*") if p.is_dir()])
    if len(dirs) != 20:
        raise RuntimeError(f"expected 20 reserve spatial shard dirs, found {len(dirs)}")
    frames = []
    null_by_taxon: dict[int, np.ndarray] = {}
    for d in dirs:
        receipt = json.loads((d / "receipt.json").read_text(encoding="utf-8"))
        if receipt.get("metrics") != RES_METRICS or int(receipt.get("permutations", -1)) != SPATIAL_NULL:
            raise RuntimeError(f"reserve spatial receipt drift: {d}")
        sp = pd.read_csv(d / "species.csv")
        z = np.load(d / "null.npz")
        tids = z["taxon_ids"].astype(int)
        arr = z["null"].astype(float)
        if arr.shape != (len(tids), len(RES_METRICS), SPATIAL_NULL):
            raise RuntimeError(f"reserve spatial null shape drift: {arr.shape}")
        frames.append(sp)
        for i, tid in enumerate(tids):
            if int(tid) in null_by_taxon:
                raise RuntimeError(f"duplicate reserve taxon {tid}")
            null_by_taxon[int(tid)] = arr[i]
    species = pd.concat(frames, ignore_index=True).sort_values("inat_taxon_id").reset_index(drop=True)
    if len(species) != EXPECTED_RES or species["inat_taxon_id"].duplicated().any():
        raise RuntimeError("reserve spatial species drift")
    tids = species["inat_taxon_id"].astype(int).to_numpy()
    if set(tids) != set(null_by_taxon):
        raise RuntimeError("reserve spatial null identity mismatch")
    null = np.stack([null_by_taxon[int(t)] for t in tids], axis=0)
    return species, null


def spatial_null_test(d: np.ndarray, observed_y: np.ndarray, null_y: np.ndarray, span: np.ndarray, technical: np.ndarray) -> dict[str, Any]:
    observed = partial_rank(d, observed_y, [span, technical])
    vals = np.array([partial_rank(d, null_y[:, k], [span, technical]) for k in range(SPATIAL_NULL)], float)
    p = float((1 + np.sum(vals >= observed - 1e-15)) / (SPATIAL_NULL + 1))
    return {
        "observed_partial_spearman": observed,
        "p_upper_geometry_preserving_spatial_null": p,
        "null_mean": float(vals.mean()),
        "null_q025": float(np.quantile(vals, 0.025)),
        "null_q975": float(np.quantile(vals, 0.975)),
        "positive_supported_at_0_05": bool(observed > 0 and p < 0.05),
    }


def diagnostic_block(df: pd.DataFrame, *, span_col: str, seed_offset: int) -> dict[str, Any]:
    D = df["D"].to_numpy(float)
    genus = df["genus"]
    return {
        "rho_D_technical_failure_rate": float(stats.spearmanr(D, df["technical_failure_rate"]).statistic),
        "rho_D_ambiguous_palette_rate": float(stats.spearmanr(D, df["ambiguous_palette_rate"]).statistic),
        "rho_D_no_biological_palette_rate": float(stats.spearmanr(D, df["no_biological_palette_rate"]).statistic),
        "rho_D_classifiable_rate": float(stats.spearmanr(D, df["classifiable_rate"]).statistic),
        "rho_D_sampled_span": float(stats.spearmanr(D, df[span_col]).statistic),
        "technical_failure_genus_clustering": clustering_test(df["technical_failure_rate"].to_numpy(float), genus, seed=SEED + seed_offset + 1),
        "ambiguous_palette_genus_clustering": clustering_test(df["ambiguous_palette_rate"].to_numpy(float), genus, seed=SEED + seed_offset + 2),
    }


def genus_adjusted_block(df: pd.DataFrame, *, span_col: str, seed_offset: int) -> dict[str, Any]:
    D = df["D"].to_numpy(float)
    n = df["n_classifiable"].to_numpy(float)
    Dunb = D * n / (n - 1.0)
    span = df[span_col].to_numpy(float)
    tech = df["technical_failure_rate"].to_numpy(float)
    rD = rank_residual(D, [span, tech])
    rDu = rank_residual(Dunb, [span, tech])
    return {
        "span_plus_technical_adjusted_D": clustering_test(rD, df["genus"], seed=SEED + seed_offset + 1),
        "span_plus_technical_adjusted_D_unbiased": clustering_test(rDu, df["genus"], seed=SEED + seed_offset + 2),
    }


def main() -> None:
    args = parse_args()
    pre = json.loads(PREFLIGHT.read_text(encoding="utf-8"))
    step7b = json.loads(STEP7B.read_text(encoding="utf-8"))
    if pre.get("D_association_computed") is not False or pre.get("every_species_exactly_100_rows") is not True:
        raise RuntimeError("Step8 preflight firewall/fixed-denominator receipt invalid")

    disc_status = build_status_panel(DISC_RAW)
    res_status = build_status_panel(RES_RAW)
    disc_attr = add_genus(pd.read_csv(DISC_ATTR)).sort_values("species").reset_index(drop=True)
    res_attr = add_genus(pd.read_csv(RES_ATTR)).sort_values("species").reset_index(drop=True)
    if len(disc_attr) != EXPECTED_DISC or len(res_attr) != EXPECTED_RES:
        raise RuntimeError("eligible D frame count drift")

    discovery = disc_attr.merge(disc_status, on="species", how="inner", validate="one_to_one")
    reserve = res_attr.merge(res_status, on="species", how="inner", validate="one_to_one")
    if len(discovery) != EXPECTED_DISC or len(reserve) != EXPECTED_RES:
        raise RuntimeError("status/D join lost eligible species")
    if not np.array_equal(discovery["n_classifiable"].astype(int).to_numpy(), discovery["classifiable_count"].astype(int).to_numpy()):
        raise RuntimeError("discovery n_classifiable != exact classified_four_state_morph count")
    if not np.array_equal(reserve["n_classifiable"].astype(int).to_numpy(), reserve["classifiable_count"].astype(int).to_numpy()):
        raise RuntimeError("reserve n_classifiable != exact classified_four_state_morph count")

    disc_diag = diagnostic_block(discovery, span_col="log1p_span_primary", seed_offset=100)
    res_diag = diagnostic_block(reserve, span_col="reserve_log1p_span", seed_offset=200)
    disc_genus = genus_adjusted_block(discovery, span_col="log1p_span_primary", seed_offset=300)
    res_genus = genus_adjusted_block(reserve, span_col="reserve_log1p_span", seed_offset=400)

    disc_obs_null, disc_null = load_discovery_spatial_null(args.discovery_artifact_root)
    disc_sp = pd.read_csv(DISC_SPATIAL)[["species", "spatial_observed_rho"]]
    disc_check = disc_obs_null.merge(disc_sp, on="species", validate="one_to_one", suffixes=("_null", "_saved")).sort_values("species")
    if np.max(np.abs(disc_check["spatial_observed_rho_null"] - disc_check["spatial_observed_rho_saved"])) > 1e-14:
        raise RuntimeError("discovery observed spatial anchor mismatch")
    dspace = discovery.sort_values("species").reset_index(drop=True)
    if not np.array_equal(dspace["species"].astype(str).to_numpy(), disc_obs_null.sort_values("species")["species"].astype(str).to_numpy()):
        raise RuntimeError("discovery spatial null order mismatch")
    dd = dspace["D"].to_numpy(float)
    nd = dspace["n_classifiable"].to_numpy(float)
    dud = dd * nd / (nd - 1.0)
    sd = dspace["log1p_span_primary"].to_numpy(float)
    td = dspace["technical_failure_rate"].to_numpy(float)
    yd = disc_obs_null.sort_values("species")["spatial_observed_rho"].to_numpy(float)
    disc_spatial = spatial_null_test(dd, yd, disc_null, sd, td)
    disc_spatial_unb = spatial_null_test(dud, yd, disc_null, sd, td)

    res_obs, res_null = load_reserve_spatial_null(args.reserve_artifact_root)
    res_saved = pd.read_csv(RES_SPATIAL).sort_values("inat_taxon_id").reset_index(drop=True)
    rspace = reserve.sort_values("inat_taxon_id_status").reset_index(drop=True)
    # Eligible reserve status taxon IDs must match the frozen spatial inference frame.
    if not np.array_equal(rspace["inat_taxon_id_status"].astype(int).to_numpy(), res_obs["inat_taxon_id"].astype(int).to_numpy()):
        raise RuntimeError("reserve spatial null/status taxon order mismatch")
    chk = res_obs[["inat_taxon_id", "rho_primary", "rho_matched_background_differential"]].merge(
        res_saved[["inat_taxon_id", "rho_primary", "rho_matched_background_differential"]],
        on="inat_taxon_id", validate="one_to_one", suffixes=("_artifact", "_saved")
    )
    for col in ["rho_primary", "rho_matched_background_differential"]:
        if np.max(np.abs(chk[f"{col}_artifact"] - chk[f"{col}_saved"])) > 1e-14:
            raise RuntimeError(f"reserve observed spatial anchor mismatch: {col}")
    dr = rspace["D"].to_numpy(float)
    nr = rspace["n_classifiable"].to_numpy(float)
    dur = dr * nr / (nr - 1.0)
    sr = rspace["reserve_log1p_span"].to_numpy(float)
    tr = rspace["technical_failure_rate"].to_numpy(float)
    yr = res_obs["rho_primary"].to_numpy(float)
    ybg = res_obs["rho_matched_background_differential"].to_numpy(float)
    res_spatial_primary = spatial_null_test(dr, yr, res_null[:, 0, :], sr, tr)
    res_spatial_primary_unb = spatial_null_test(dur, yr, res_null[:, 0, :], sr, tr)
    res_spatial_bg = spatial_null_test(dr, ybg, res_null[:, 3, :], sr, tr)
    res_spatial_bg_unb = spatial_null_test(dur, ybg, res_null[:, 3, :], sr, tr)

    genus_survives = bool(res_genus["span_plus_technical_adjusted_D"]["supported_at_0_05"])
    genus_unb_survives = bool(res_genus["span_plus_technical_adjusted_D_unbiased"]["supported_at_0_05"])
    spatial_reserve_survives = bool(res_spatial_primary["positive_supported_at_0_05"])
    spatial_bg_survives = bool(res_spatial_bg["positive_supported_at_0_05"])

    result = {
        "analysis": "polymorphism_measurement_missingness_step8_analysis",
        "date_jst": "2026-09-10",
        "protocol": "docs/POLYMORPHISM_MEASUREMENT_MISSINGNESS_STEP8_ANALYSIS_PROTOCOL_20260910.md",
        "measurement_class_definition": STATUS_TO_RATE,
        "fixed_raw_photo_denominator": {"discovery": 100, "reserve": 100},
        "n_permutations_genus": N_PERM,
        "spatial_null_permutations_reused": SPATIAL_NULL,
        "new_spatial_permutations_generated": False,
        "discovery_diagnostics": disc_diag,
        "reserve_diagnostics": res_diag,
        "discovery_genus": disc_genus,
        "reserve_genus": res_genus,
        "discovery_spatial": {"D": disc_spatial, "D_unbiased": disc_spatial_unb},
        "reserve_spatial_primary": {"D": res_spatial_primary, "D_unbiased": res_spatial_primary_unb},
        "reserve_spatial_matched_background_differential": {"D": res_spatial_bg, "D_unbiased": res_spatial_bg_unb},
        "comparison_to_step7b": {
            "reserve_span_only": step7b["reserve"]["span_only_adjusted_D"],
            "reserve_n_classifiable_only": step7b["reserve"]["n_classifiable_only_adjusted_D"],
        },
        "decision": {
            "reserve_genus_clustering_survives_span_plus_pure_technical_failure": genus_survives,
            "reserve_genus_clustering_D_unbiased_survives_span_plus_pure_technical_failure": genus_unb_survives,
            "reserve_primary_D_spatial_survives_span_plus_pure_technical_failure": spatial_reserve_survives,
            "reserve_flower_minus_background_D_spatial_survives_span_plus_pure_technical_failure": spatial_bg_survives,
            "pure_technical_roi_flip_failure_sufficient_explanation_for_genus_signal": False if genus_survives else None,
            "pure_technical_roi_flip_failure_sufficient_explanation_for_flower_specific_spatial_signal": False if spatial_bg_survives else None,
            "claim_boundary": "Survival rules out sampled span plus clearly technical ROI/flip failure as a sufficient explanation. Ambiguous-palette and no-biological-palette missingness remain unresolved and are not recoded as hidden morphs.",
        },
    }

    discovery.to_csv(OUT / "discovery_measurement_process_panel.csv", index=False)
    reserve.to_csv(OUT / "reserve_measurement_process_panel.csv", index=False)
    write_json(OUT / "result.json", result)

    def c_line(label: str, x: dict[str, Any]) -> str:
        return f"- {label}: gain **{x['clustering_gain']:.6f}**, p_lower **{x['p_lower']:.6g}**"

    lines = [
        "# Polymorphism measurement missingness — Step 8 analysis result",
        "",
        f"**Reserve genus survives span + pure technical-failure control: `{genus_survives}`.**",
        f"**Reserve genus D_unbiased survives span + pure technical-failure control: `{genus_unb_survives}`.**",
        f"**Reserve primary D-spatial survives span + pure technical-failure control: `{spatial_reserve_survives}`.**",
        f"**Reserve flower-minus-background D-spatial survives span + pure technical-failure control: `{spatial_bg_survives}`.**",
        "",
        "## What `technical failure` means",
        "",
        "Only `not_evaluable_roi_or_flip_gate` is counted as technical failure. Ambiguous palette and no-biological-palette rows are kept separate.",
        "",
        "## Measurement-process diagnostics",
        "",
        f"- discovery rho(D, technical failure) = **{disc_diag['rho_D_technical_failure_rate']:.6f}**",
        f"- discovery rho(D, ambiguous palette) = **{disc_diag['rho_D_ambiguous_palette_rate']:.6f}**",
        f"- reserve rho(D, technical failure) = **{res_diag['rho_D_technical_failure_rate']:.6f}**",
        f"- reserve rho(D, ambiguous palette) = **{res_diag['rho_D_ambiguous_palette_rate']:.6f}**",
        "",
        "## Genus clustering after sampled span + technical failure",
        "",
        c_line("discovery D", disc_genus["span_plus_technical_adjusted_D"]),
        c_line("discovery D_unbiased", disc_genus["span_plus_technical_adjusted_D_unbiased"]),
        c_line("reserve D", res_genus["span_plus_technical_adjusted_D"]),
        c_line("reserve D_unbiased", res_genus["span_plus_technical_adjusted_D_unbiased"]),
        "",
        "## Geometry-preserving D-spatial tests after sampled span + technical failure",
        "",
        f"- discovery primary: partial rho **{disc_spatial['observed_partial_spearman']:.6f}**, p_upper **{disc_spatial['p_upper_geometry_preserving_spatial_null']:.6g}**",
        f"- reserve primary: partial rho **{res_spatial_primary['observed_partial_spearman']:.6f}**, p_upper **{res_spatial_primary['p_upper_geometry_preserving_spatial_null']:.6g}**",
        f"- reserve flower-minus-background: partial rho **{res_spatial_bg['observed_partial_spearman']:.6f}**, p_upper **{res_spatial_bg['p_upper_geometry_preserving_spatial_null']:.6g}**",
        "",
        "## Interpretation boundary",
        "",
        "These tests separate clearly operational ROI/flip failures from ambiguity-related missingness. A surviving signal is not explained by sampled geographic span plus the rate of clear technical failures alone. Ambiguous palette rows remain unresolved measurement/biological ambiguity and are not treated as hidden morphs.",
    ]
    (OUT / "RESULT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(result["decision"], indent=2))
    print((OUT / "RESULT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
