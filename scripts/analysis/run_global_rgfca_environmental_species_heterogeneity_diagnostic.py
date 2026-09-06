#!/usr/bin/env python3
"""Omnibus post-outcome species-level environmental response heterogeneity diagnostic."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/global_rgfca_environmental_species_heterogeneity_diagnostic_contract_v1.json"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--artifact-root", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def find_species_file(root: Path, artifact_name: str, interaction: bool = False) -> Path:
    base = root / artifact_name
    if not base.exists():
        raise RuntimeError(f"missing artifact directory: {base}")
    pattern = "*interaction_family_species_v1.csv" if interaction else "*block_species_v1.csv"
    files = sorted(base.rglob(pattern))
    if len(files) != 1:
        raise RuntimeError(f"expected exactly one species file in {artifact_name}, found {len(files)}: {files}")
    return files[0]


def load_effect_matrix(contract: dict[str, object], artifact_root: Path) -> tuple[pd.DataFrame, dict[str, str]]:
    parts = []
    input_sha = {}
    for entry in contract["inputs"]["main_effect_artifacts"]:
        block = str(entry["block"])
        path = find_species_file(artifact_root, str(entry["artifact_name"]))
        digest = sha256_file(path)
        if digest != str(entry["species_sha256"]):
            raise RuntimeError(f"species SHA mismatch for {block}: {digest}")
        df = pd.read_csv(path)
        effect_cols = [c for c in df.columns if c.startswith("partial_rho_colour_vs_") and c.endswith("_given_distance")]
        if len(effect_cols) != 1:
            raise RuntimeError(f"cannot identify one effect column for {block}: {effect_cols}")
        x = df[["species", effect_cols[0]]].copy()
        x["species"] = x["species"].astype(str)
        x = x.rename(columns={effect_cols[0]: f"main__{block}"})
        if x["species"].duplicated().any():
            raise RuntimeError(f"duplicate species in {block}")
        parts.append(x.set_index("species"))
        input_sha[block] = digest

    inter = contract["inputs"]["interaction_artifact"]
    ipath = find_species_file(artifact_root, str(inter["artifact_name"]), interaction=True)
    idigest = sha256_file(ipath)
    if idigest != str(inter["species_sha256"]):
        raise RuntimeError(f"interaction species SHA mismatch: {idigest}")
    raw = pd.read_csv(ipath)
    required = {"interaction", "species", "beta_interaction"}
    missing = sorted(required - set(raw.columns))
    if missing:
        raise RuntimeError(f"interaction species file missing {missing}")
    raw["species"] = raw["species"].astype(str)
    pivot = raw.pivot(index="species", columns="interaction", values="beta_interaction")
    if pivot.shape[1] != 10:
        raise RuntimeError(f"expected 10 interaction features, found {pivot.shape[1]}")
    pivot.columns = [f"interaction__{c}" for c in pivot.columns]
    parts.append(pivot)
    input_sha["ten_interaction_family"] = idigest

    matrix = pd.concat(parts, axis=1, join="inner").sort_index()
    expected_species = int(contract["inputs"]["expected_species"])
    expected_features = int(contract["inputs"]["expected_features"])
    if matrix.shape != (expected_species, expected_features):
        raise RuntimeError(f"effect matrix shape {matrix.shape} != {(expected_species, expected_features)}")
    if not np.isfinite(matrix.to_numpy(float)).all():
        raise RuntimeError("effect matrix contains non-finite values")
    return matrix, input_sha


def standardize(matrix: pd.DataFrame) -> np.ndarray:
    x = matrix.to_numpy(float)
    mean = x.mean(axis=0)
    sd = x.std(axis=0, ddof=0)
    if np.any(sd <= 0):
        raise RuntimeError("constant effect feature prevents standardization")
    return (x - mean) / sd


def pca(z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    u, s, vt = np.linalg.svd(z, full_matrices=False)
    explained = np.square(s) / np.sum(np.square(s))
    scores = z @ vt.T
    return vt, explained, scores


def match_components(reference: np.ndarray, boot: np.ndarray, n_components: int = 3) -> list[tuple[int, float]]:
    available = set(range(boot.shape[0]))
    matches = []
    for r in range(n_components):
        candidates = sorted(available)
        dots = np.asarray([float(np.dot(reference[r], boot[j])) for j in candidates])
        best_local = int(np.argmax(np.abs(dots)))
        best = candidates[best_local]
        sign = 1.0 if dots[best_local] >= 0 else -1.0
        matches.append((best, sign))
        available.remove(best)
    return matches


def bootstrap_pca(z: np.ndarray, reference: np.ndarray, reps: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    n, p = z.shape
    loading = np.empty((reps, 3, p), dtype=np.float32)
    explained = np.empty((reps, 3), dtype=np.float32)
    for b in range(reps):
        idx = rng.integers(0, n, size=n)
        sample = z[idx]
        sample = (sample - sample.mean(axis=0)) / sample.std(axis=0, ddof=0)
        vt, ev, _ = pca(sample)
        for r, (j, sign) in enumerate(match_components(reference, vt, 3)):
            loading[b, r] = (vt[j] * sign).astype(np.float32)
            explained[b, r] = float(ev[j])
    return loading, explained


def main() -> int:
    args = parse_args()
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    if contract.get("status") != "postoutcome_exploratory_omnibus_frozen_before_species_effect_matrix_is_reopened":
        raise RuntimeError("heterogeneity diagnostic contract drift")
    matrix, input_sha = load_effect_matrix(contract, args.artifact_root)
    z = standardize(matrix)
    features = list(matrix.columns)

    corr = pd.DataFrame(z, index=matrix.index, columns=features).corr(method="spearman")
    vt, explained, scores = pca(z)
    loading_df = pd.DataFrame(vt[:5].T, index=features, columns=[f"PC{i}" for i in range(1, 6)])
    loading_df.index.name = "feature"

    boot_cfg = contract["analysis"]["bootstrap_stability"]
    boot_loading, boot_explained = bootstrap_pca(
        z, vt, int(boot_cfg["replicates"]), int(boot_cfg["seed"])
    )
    boot_rows = []
    for pc in range(3):
        for j, feature in enumerate(features):
            vals = boot_loading[:, pc, j]
            boot_rows.append({
                "component": f"PC{pc+1}",
                "feature": feature,
                "reference_loading": float(vt[pc, j]),
                "loading_q10": float(np.quantile(vals, 0.10)),
                "loading_q50": float(np.quantile(vals, 0.50)),
                "loading_q90": float(np.quantile(vals, 0.90)),
            })
    boot_df = pd.DataFrame(boot_rows)

    species_out = matrix.copy()
    for pc in range(5):
        species_out[f"PC{pc+1}_score"] = scores[:, pc]
    species_out = species_out.reset_index()

    corr_long = corr.rename_axis(index="feature_a", columns="feature_b").stack().rename("spearman_rho").reset_index()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    species_out.to_csv(args.output_dir / "global_rgfca_environmental_species_effect_matrix_v1.csv", index=False)
    corr_long.to_csv(args.output_dir / "global_rgfca_environmental_species_effect_correlations_v1.csv", index=False)
    loading_df.reset_index().to_csv(args.output_dir / "global_rgfca_environmental_species_effect_pca_loadings_v1.csv", index=False)
    boot_df.to_csv(args.output_dir / "global_rgfca_environmental_species_effect_pca_bootstrap_v1.csv", index=False)

    thermal = "main__thermal_regime"
    dryness = "main__atmospheric_energy_dryness"
    td_corr = float(corr.loc[thermal, dryness])
    ev_boot = {
        f"PC{i+1}": {
            "reference": float(explained[i]),
            "q10": float(np.quantile(boot_explained[:, i], 0.10)),
            "q50": float(np.quantile(boot_explained[:, i], 0.50)),
            "q90": float(np.quantile(boot_explained[:, i], 0.90)),
        }
        for i in range(3)
    }
    result = {
        "protocol": contract["protocol"],
        "status": "complete_omnibus_species_environmental_heterogeneity_diagnostic",
        "n_species": int(matrix.shape[0]),
        "n_features": int(matrix.shape[1]),
        "features": features,
        "explained_variance_fraction_PC1_to_PC5": {f"PC{i+1}": float(explained[i]) for i in range(5)},
        "bootstrap_explained_variance_PC1_to_PC3": ev_boot,
        "thermal_vs_atmospheric_dryness_species_effect_spearman": td_corr,
        "thermal_dryness_interpretation_guard": "Descriptive member of the complete 15-feature matrix; not a confirmatory combined heat-dryness test and does not open variable decomposition.",
        "input_species_sha256": input_sha,
        "decision_guard": {
            "five_block_main_effect_panel_reclassified": False,
            "ten_interaction_family_reclassified": False,
            "individual_variable_decomposition_opened": False,
            "causal_or_local_adaptation_language_allowed": False,
        },
        "files": {
            "species_matrix": "data/derived/global_rgfca_environmental_species_effect_matrix_v1.csv",
            "correlations": "data/derived/global_rgfca_environmental_species_effect_correlations_v1.csv",
            "pca_loadings": "data/derived/global_rgfca_environmental_species_effect_pca_loadings_v1.csv",
            "pca_bootstrap": "data/derived/global_rgfca_environmental_species_effect_pca_bootstrap_v1.csv"
        }
    }
    (args.output_dir / "global_rgfca_environmental_species_heterogeneity_diagnostic_result_v1.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
