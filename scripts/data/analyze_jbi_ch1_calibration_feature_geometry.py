#!/usr/bin/env python3
"""Diagnose within-species calibration feature geometry without assigning final morphs.

The input is the quota-independent 480-photo Florence feature table. This analysis is
calibration-only and deliberately ignores geography, observer, date, environment and
held-out evaluation rows. Gaussian-mixture component support is summarized with BIC
and bootstrap selection frequencies, but components are not converted into biological
colour states here.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import warnings

import numpy as np
from sklearn.exceptions import ConvergenceWarning
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler

PROTOCOL = "jbi-ch1-calibration-feature-geometry-v1"
RANDOM_SEED = 20260829
BOOTSTRAPS = 200

# Maximum component counts come from the predeclared literature-constrained candidate
# state structure. They are caps for diagnostic exploration, not forced component counts.
SPECIES_SPECS = {
    "Ipomoea purpurea": {
        "source": "candidate_scores",
        "features": ["white", "pink", "blue_purple"],
        "max_components": 3,
    },
    "Raphanus sativus": {
        "source": "visual_colour_axes",
        "features": ["anthocyanin_like_signal", "carotenoid_like_signal"],
        "max_components": 4,
    },
    "Gentiana lutea": {
        "source": "candidate_scores",
        "features": ["yellow", "orange"],
        "max_components": 2,
    },
    "Dactylorhiza sambucina": {
        "source": "candidate_scores",
        "features": ["yellow", "purple"],
        "max_components": 2,
    },
    "Antirrhinum majus": {
        "source": "candidate_scores",
        "features": ["magenta_pseudomajus_like", "yellow_striatum_like"],
        "max_components": 3,
    },
    "Lysimachia arvensis": {
        "source": "candidate_scores",
        "features": ["blue", "red"],
        "max_components": 2,
    },
}


def load_rows(path: Path) -> list[dict]:
    rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(rows) != 480:
        raise ValueError(f"expected 480 calibration feature rows, found {len(rows)}")
    if len({str(r["photo_id"]) for r in rows}) != 480:
        raise ValueError("duplicate photo IDs in calibration features")
    if any(r.get("evaluation_row") is not False or r.get("final_label") is not False for r in rows):
        raise ValueError("evaluation/final-label firewall violation")
    return rows


def feature_vector(row: dict, species: str) -> list[float]:
    spec = SPECIES_SPECS[species]
    payload = row.get(spec["source"], {})
    if not isinstance(payload, dict):
        raise ValueError(f"{species}: missing {spec['source']} payload")
    try:
        values = [float(payload[name]) for name in spec["features"]]
    except KeyError as exc:
        raise ValueError(f"{species}: missing feature {exc.args[0]}") from exc
    if not np.isfinite(values).all():
        raise ValueError(f"{species}: non-finite feature values")
    return values


def fit_bic_grid(x: np.ndarray, max_components: int, seed: int = RANDOM_SEED) -> list[dict]:
    out = []
    for k in range(1, max_components + 1):
        model = GaussianMixture(
            n_components=k,
            covariance_type="full",
            reg_covar=1e-5,
            n_init=20,
            random_state=seed + k,
            max_iter=500,
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            model.fit(x)
        out.append({
            "components": k,
            "bic": float(model.bic(x)),
            "aic": float(model.aic(x)),
            "converged": bool(model.converged_),
        })
    return out


def best_bic_components(grid: list[dict]) -> int:
    return int(min(grid, key=lambda row: row["bic"])["components"])


def bootstrap_component_support(
    x: np.ndarray,
    max_components: int,
    n_bootstraps: int = BOOTSTRAPS,
    seed: int = RANDOM_SEED,
) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    counts: Counter[int] = Counter()
    n = len(x)
    for b in range(n_bootstraps):
        idx = rng.integers(0, n, size=n)
        sample = x[idx]
        grid = fit_bic_grid(sample, max_components, seed=seed + 1000 + b * 17)
        counts[best_bic_components(grid)] += 1
    return {
        str(k): counts[k] / n_bootstraps
        for k in range(1, max_components + 1)
    }


def quantiles(x: np.ndarray, feature_names: list[str]) -> dict[str, dict[str, float]]:
    probs = [0.05, 0.25, 0.5, 0.75, 0.95]
    q = np.quantile(x, probs, axis=0)
    result = {}
    for j, name in enumerate(feature_names):
        result[name] = {
            "p05": float(q[0, j]),
            "p25": float(q[1, j]),
            "p50": float(q[2, j]),
            "p75": float(q[3, j]),
            "p95": float(q[4, j]),
        }
    return result


def analyze_species(rows: list[dict], species: str) -> dict:
    spec = SPECIES_SPECS[species]
    group = [r for r in rows if r.get("species") == species]
    if len(group) != 80:
        raise ValueError(f"{species}: expected 80 rows, found {len(group)}")
    ok = [r for r in group if r.get("feature_status") == "ok"]
    if len(ok) < 10:
        return {
            "species": species,
            "status": "insufficient_feature_rows",
            "n_total": len(group),
            "n_feature_ok": len(ok),
            "n_localization_failed": len(group) - len(ok),
            "feature_source": spec["source"],
            "features": spec["features"],
            "literature_candidate_component_cap": spec["max_components"],
        }

    raw = np.asarray([feature_vector(r, species) for r in ok], dtype=float)
    scaled = StandardScaler().fit_transform(raw)
    max_components = min(int(spec["max_components"]), max(1, len(ok) // 10))
    grid = fit_bic_grid(scaled, max_components)
    selected = best_bic_components(grid)
    bootstrap = bootstrap_component_support(scaled, max_components)
    return {
        "species": species,
        "status": "feature_geometry_measured_not_interpreted_as_morph_labels",
        "n_total": len(group),
        "n_feature_ok": len(ok),
        "n_localization_failed": len(group) - len(ok),
        "feature_source": spec["source"],
        "features": spec["features"],
        "literature_candidate_component_cap": spec["max_components"],
        "components_evaluated": list(range(1, max_components + 1)),
        "bic_grid": grid,
        "bic_selected_components": selected,
        "bootstrap_bic_selected_component_frequency": bootstrap,
        "raw_feature_quantiles": quantiles(raw, list(spec["features"])),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/calibration/jbi_ch1_florence_calibration_features_v1.jsonl"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("docs/supporting/jbi_ch1_calibration_feature_geometry_v1.json"),
    )
    args = parser.parse_args()

    rows = load_rows(args.input)
    species = sorted({str(r["species"]) for r in rows})
    if species != sorted(SPECIES_SPECS):
        raise ValueError(f"unexpected species set: {species}")

    result = {
        "protocol": PROTOCOL,
        "status": "pre_condition_feature_geometry_complete_not_final_labels",
        "calibration_only": True,
        "evaluation_rows_opened": False,
        "final_label": False,
        "geography_used": False,
        "observer_used": False,
        "date_used": False,
        "environment_used": False,
        "candidate_argmax_used_as_training_label": False,
        "condition_filter_applied": False,
        "gmm_components_are_biological_morph_labels": False,
        "bootstrap_replicates": BOOTSTRAPS,
        "species": [analyze_species(rows, sp) for sp in species],
        "next_gate": "repeat feature-geometry diagnosis after an independently validated fresh/senescent calibration condition filter; only then freeze species-specific measurement rules",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
