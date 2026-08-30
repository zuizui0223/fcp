#!/usr/bin/env python3
"""Run the preregistered Chapter 1 Stage-A spatial test.

Primary analysis
----------------
* one colour-blind spherical kNN graph per species (k=5),
* complete continuous colour vectors kept intact,
* within-species vector permutations (9,999),
* species-specific graph discontinuities combined with equal species weight.

Sensitivity analyses repeat the identical procedure at k=3 and k=8.  This script does
not tune a threshold, choose clusters, create colour classes, or inspect biological
outcomes.  A representation is selected only by a fixed, name-based priority list.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np

from fcp_pipeline.spatial_graph import (
    species_conditioned_graph_permutation_null,
    species_conditioned_knn_graph,
)

EXPECTED_EVALUATION_ROWS = 720
PRIMARY_K = 5
SENSITIVITY_K = (3, 8)
DEFAULT_PERMUTATIONS = 9_999
DEFAULT_SEED = 20260830

# Fixed priority.  Candidate order is part of the analysis contract and never depends on
# observed values, geography, species effects, or a result from the permutation test.
STANDARDIZED_VECTOR_PATHS: tuple[str, ...] = (
    "standardized_colour_vector",
    "standardized_color_vector",
    "colour_z_vector",
    "color_z_vector",
    "z_colour_vector",
    "z_color_vector",
    "features.standardized_colour_vector",
    "features.standardized_color_vector",
    "colour_features.standardized_vector",
    "color_features.standardized_vector",
    "measurement.standardized_colour_vector",
    "measurement.standardized_color_vector",
)
RAW_VECTOR_PATHS: tuple[str, ...] = (
    "continuous_colour_vector",
    "continuous_color_vector",
    "colour_vector",
    "color_vector",
    "cielab",
    "lab",
    "mean_cielab",
    "mean_lab",
    "flower_cielab",
    "flower_lab",
    "roi_cielab",
    "roi_lab",
    "features.continuous_colour_vector",
    "features.continuous_color_vector",
    "features.colour_vector",
    "features.color_vector",
    "features.cielab",
    "features.lab",
    "colour_features.vector",
    "color_features.vector",
    "measurement.colour_vector",
    "measurement.color_vector",
)
STANDARDIZED_SCALAR_SETS: tuple[tuple[str, ...], ...] = (
    ("z_L", "z_a", "z_b"),
    ("L_z", "a_z", "b_z"),
    ("lab_L_z", "lab_a_z", "lab_b_z"),
    ("cielab_L_z", "cielab_a_z", "cielab_b_z"),
    ("features.z_L", "features.z_a", "features.z_b"),
    ("features.L_z", "features.a_z", "features.b_z"),
    ("z_lightness", "z_chroma", "z_hue_cos", "z_hue_sin"),
    ("features.z_lightness", "features.z_chroma", "features.z_hue_cos", "features.z_hue_sin"),
)
RAW_SCALAR_SETS: tuple[tuple[str, ...], ...] = (
    ("L", "a", "b"),
    ("lab_L", "lab_a", "lab_b"),
    ("cielab_L", "cielab_a", "cielab_b"),
    ("mean_L", "mean_a", "mean_b"),
    ("mean_lab_L", "mean_lab_a", "mean_lab_b"),
    ("features.L", "features.a", "features.b"),
    ("features.lab_L", "features.lab_a", "features.lab_b"),
    ("features.cielab_L", "features.cielab_a", "features.cielab_b"),
    ("lightness", "chroma", "hue_cos", "hue_sin"),
    ("features.lightness", "features.chroma", "features.hue_cos", "features.hue_sin"),
)


@dataclass(frozen=True)
class Representation:
    kind: str
    paths: tuple[str, ...]
    standardized: bool
    dimension: int

    @property
    def label(self) -> str:
        return f"{self.kind}:" + ",".join(self.paths)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def scalar_id(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text.endswith(".0"):
        try:
            return str(int(float(text)))
        except ValueError:
            pass
    return text


def get_path(row: dict[str, Any], path: str) -> Any:
    value: Any = row
    for part in path.split("."):
        if not isinstance(value, dict) or part not in value:
            raise KeyError(path)
        value = value[part]
    return value


def parse_vector(value: Any) -> np.ndarray:
    if isinstance(value, str):
        text = value.strip()
        if not text:
            raise ValueError("empty vector string")
        try:
            value = json.loads(text)
        except json.JSONDecodeError:
            value = [part.strip() for part in text.replace(";", ",").split(",") if part.strip()]
    if isinstance(value, dict):
        # Fixed conventional orders only; never alphabetical selection of arbitrary fields.
        for keys in (("L", "a", "b"), ("l", "a", "b"), ("lightness", "chroma", "hue_cos", "hue_sin")):
            if all(key in value for key in keys):
                value = [value[key] for key in keys]
                break
        else:
            raise ValueError("unrecognized vector object")
    array = np.asarray(value, dtype=float)
    if array.ndim != 1 or len(array) < 2:
        raise ValueError("colour vector must be one-dimensional with at least two components")
    if not np.isfinite(array).all():
        raise ValueError("colour vector contains non-finite values")
    return array


def extract_with_representation(row: dict[str, Any], rep: Representation) -> np.ndarray:
    if rep.kind == "vector":
        vector = parse_vector(get_path(row, rep.paths[0]))
    elif rep.kind == "scalars":
        vector = np.asarray([float(get_path(row, path)) for path in rep.paths], dtype=float)
    else:
        raise RuntimeError(f"unknown representation kind {rep.kind}")
    if vector.shape != (rep.dimension,) or not np.isfinite(vector).all():
        raise ValueError(f"row does not match frozen representation {rep.label}")
    return vector


def discover_representation(rows: Sequence[dict[str, Any]]) -> Representation:
    if not rows:
        raise ValueError("cannot discover a representation from zero rows")

    for standardized, vector_paths, scalar_sets in (
        (True, STANDARDIZED_VECTOR_PATHS, STANDARDIZED_SCALAR_SETS),
        (False, RAW_VECTOR_PATHS, RAW_SCALAR_SETS),
    ):
        for path in vector_paths:
            try:
                vectors = [parse_vector(get_path(row, path)) for row in rows]
            except (KeyError, TypeError, ValueError):
                continue
            dimensions = {len(vector) for vector in vectors}
            if len(dimensions) == 1:
                return Representation("vector", (path,), standardized, dimensions.pop())
        for paths in scalar_sets:
            try:
                matrix = np.asarray(
                    [[float(get_path(row, path)) for path in paths] for row in rows],
                    dtype=float,
                )
            except (KeyError, TypeError, ValueError):
                continue
            if matrix.shape == (len(rows), len(paths)) and np.isfinite(matrix).all():
                return Representation("scalars", tuple(paths), standardized, len(paths))

    available = sorted({key for row in rows[:50] for key in row.keys()})
    raise ValueError(
        "no preregistered continuous-colour representation was found; top-level fields were: "
        + ", ".join(available)
    )


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            value = json.loads(line)
            if not isinstance(value, dict):
                raise ValueError(f"{path}:{line_number}: expected a JSON object")
            rows.append(value)
    return rows


def read_csv(path: Path) -> list[dict[str, Any]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def candidate_calibration_files(root: Path, evaluation_path: Path) -> list[Path]:
    patterns = ("*.jsonl", "*.csv")
    candidates: list[Path] = []
    for pattern in patterns:
        for path in root.rglob(pattern):
            if path.resolve() == evaluation_path.resolve() or not path.is_file():
                continue
            lower = str(path).lower()
            if any(token in lower for token in ("calibration", "colour_feature", "color_feature", "florence")):
                if not any(token in lower for token in ("evaluation_720", "stage_a_null", "source_shards")):
                    candidates.append(path)
    return sorted(set(candidates))


def load_calibration_rows(
    root: Path,
    evaluation_path: Path,
    rep: Representation,
    expected_species: set[str],
) -> tuple[list[dict[str, Any]], list[str]]:
    selected: dict[tuple[str, str], dict[str, Any]] = {}
    used_files: list[str] = []
    for path in candidate_calibration_files(root, evaluation_path):
        try:
            rows = read_jsonl(path) if path.suffix.lower() == ".jsonl" else read_csv(path)
        except Exception:
            continue
        accepted = 0
        for index, row in enumerate(rows):
            species = str(row.get("species", "")).strip()
            if species not in expected_species:
                continue
            split = str(row.get("split", "")).strip().lower()
            is_calibration = (
                split == "calibration"
                or row.get("calibration_only") is True
                or str(row.get("calibration_only", "")).strip().lower() == "true"
            )
            if not is_calibration:
                continue
            try:
                extract_with_representation(row, rep)
            except Exception:
                continue
            photo_id = scalar_id(row.get("photo_id")) or f"{path}:{index}"
            key = (species, photo_id)
            if key not in selected:
                selected[key] = row
                accepted += 1
        if accepted:
            used_files.append(str(path.relative_to(root)))
    rows = list(selected.values())
    counts = Counter(str(row.get("species", "")).strip() for row in rows)
    missing = sorted(species for species in expected_species if counts[species] < 2)
    if missing:
        raise ValueError(
            "raw evaluation vectors require frozen calibration measurements; fewer than two "
            f"compatible calibration rows were found for: {', '.join(missing)}"
        )
    return rows, used_files


def standardize_from_calibration(
    evaluation: np.ndarray,
    eval_species: np.ndarray,
    calibration_rows: Sequence[dict[str, Any]],
    rep: Representation,
) -> tuple[np.ndarray, dict[str, dict[str, list[float] | int]]]:
    result = np.empty_like(evaluation, dtype=float)
    scalers: dict[str, dict[str, list[float] | int]] = {}
    for species in np.unique(eval_species):
        calibration = np.vstack(
            [
                extract_with_representation(row, rep)
                for row in calibration_rows
                if str(row.get("species", "")).strip() == species
            ]
        )
        mean = calibration.mean(axis=0)
        sd = calibration.std(axis=0, ddof=1)
        if not np.isfinite(mean).all() or not np.isfinite(sd).all() or np.any(sd <= 0):
            raise ValueError(f"{species}: invalid calibration mean/SD for representation {rep.label}")
        idx = eval_species == species
        result[idx] = (evaluation[idx] - mean) / sd
        scalers[str(species)] = {
            "n_calibration": int(len(calibration)),
            "mean": mean.tolist(),
            "sd": sd.tolist(),
        }
    if not np.isfinite(result).all():
        raise ValueError("standardization produced non-finite values")
    return result, scalers


def monte_carlo_upper_p(observed: float, null: np.ndarray) -> float:
    return float((1 + np.count_nonzero(null >= observed)) / (len(null) + 1))


def summarize_null(observed: float, null: np.ndarray) -> dict[str, float | str]:
    mean = float(np.mean(null))
    sd = float(np.std(null, ddof=1))
    z = float((observed - mean) / sd) if sd > 0 else math.nan
    return {
        "observed": float(observed),
        "null_mean": mean,
        "null_sd": sd,
        "null_q025": float(np.quantile(null, 0.025)),
        "null_q500": float(np.quantile(null, 0.5)),
        "null_q975": float(np.quantile(null, 0.975)),
        "effect_z": z,
        "direction": "more_discontinuous" if observed > mean else "more_spatially_clustered",
        "p_upper": monte_carlo_upper_p(observed, null),
    }


def write_csv(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    materialized = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not materialized:
        raise ValueError(f"cannot write empty CSV {path}")
    fields: list[str] = []
    for row in materialized:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(materialized)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evaluation-jsonl", type=Path, required=True)
    parser.add_argument("--calibration-root", type=Path, default=Path("."))
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--permutations", type=int, default=DEFAULT_PERMUTATIONS)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = parser.parse_args()

    if args.permutations != DEFAULT_PERMUTATIONS:
        raise ValueError(
            f"Stage A contract requires exactly {DEFAULT_PERMUTATIONS} permutations, "
            f"received {args.permutations}"
        )
    rows = read_jsonl(args.evaluation_jsonl)
    if len(rows) != EXPECTED_EVALUATION_ROWS:
        raise ValueError(
            f"Stage A requires {EXPECTED_EVALUATION_ROWS} evaluation rows, found {len(rows)}"
        )
    species = np.asarray([str(row.get("species", "")).strip() for row in rows], dtype=object)
    latitude = np.asarray([float(row["latitude"]) for row in rows], dtype=float)
    longitude = np.asarray([float(row["longitude"]) for row in rows], dtype=float)
    if any(not label for label in species) or not np.isfinite(latitude).all() or not np.isfinite(longitude).all():
        raise ValueError("species and finite coordinates are required for every evaluation row")
    counts = Counter(species.tolist())
    if len(counts) != 6 or set(counts.values()) != {120}:
        raise ValueError(f"expected six species with 120 rows each, found {dict(counts)}")

    representation = discover_representation(rows)
    raw_values = np.vstack([extract_with_representation(row, representation) for row in rows])
    calibration_files: list[str] = []
    scalers: dict[str, Any] | None = None
    if representation.standardized:
        standardized = raw_values
        scaler_provenance = "row_level_frozen_standardized_features"
    else:
        calibration_rows, calibration_files = load_calibration_rows(
            args.calibration_root,
            args.evaluation_jsonl,
            representation,
            set(counts),
        )
        standardized, scalers = standardize_from_calibration(
            raw_values,
            species,
            calibration_rows,
            representation,
        )
        scaler_provenance = "frozen_calibration_rows_species_specific_mean_sd"

    args.output_dir.mkdir(parents=True, exist_ok=True)
    species_rows: list[dict[str, Any]] = []
    global_rows: list[dict[str, Any]] = []
    null_arrays: dict[str, np.ndarray] = {}
    seeds = np.random.SeedSequence(args.seed).spawn(1 + len(SENSITIVITY_K))

    for k, seed_sequence in zip((PRIMARY_K, *SENSITIVITY_K), seeds, strict=True):
        graph = species_conditioned_knn_graph(latitude, longitude, species, k=k)
        rng = np.random.default_rng(seed_sequence)
        result = species_conditioned_graph_permutation_null(
            standardized,
            species,
            graph,
            n_permutations=args.permutations,
            rng=rng,
        )
        observed_global = float(result["observed_global_equal_species_mean"])
        null_global = np.asarray(result["null_global_equal_species_mean"], dtype=float)
        summary = summarize_null(observed_global, null_global)
        global_rows.append(
            {
                "analysis": "primary" if k == PRIMARY_K else "sensitivity",
                "k": k,
                "n_observations": len(rows),
                "n_species": len(counts),
                "n_edges": int(len(graph.edges)),
                "n_permutations": args.permutations,
                **summary,
            }
        )
        null_arrays[f"k{k}_global"] = null_global

        order = np.asarray(result["species"], dtype=object)
        observed_species = np.asarray(result["observed_species_q"], dtype=float)
        null_species = np.asarray(result["null_species_q"], dtype=float)
        null_arrays[f"k{k}_species"] = null_species
        for column, label in enumerate(order):
            species_summary = summarize_null(observed_species[column], null_species[:, column])
            species_rows.append(
                {
                    "analysis": "primary" if k == PRIMARY_K else "sensitivity",
                    "k": k,
                    "species": str(label),
                    "n_observations": counts[str(label)],
                    "n_edges": int(np.count_nonzero(graph.edge_species == label)),
                    "n_permutations": args.permutations,
                    **species_summary,
                }
            )

    global_csv = args.output_dir / "stage_a_global_summary.csv"
    species_csv = args.output_dir / "stage_a_species_summary.csv"
    null_npz = args.output_dir / "stage_a_null_distributions.npz"
    write_csv(global_csv, global_rows)
    write_csv(species_csv, species_rows)
    np.savez_compressed(null_npz, **null_arrays)

    primary = next(row for row in global_rows if row["k"] == PRIMARY_K)
    sensitivities = [row for row in global_rows if row["k"] in SENSITIVITY_K]
    same_direction = all(row["direction"] == primary["direction"] for row in sensitivities)
    sensitivity_significance = all(
        (float(row["p_upper"]) <= 0.05) == (float(primary["p_upper"]) <= 0.05)
        for row in sensitivities
    )
    manifest = {
        "protocol": "jbi_ch1_stage_a_spatial_colour_test_v1",
        "evaluation_jsonl": str(args.evaluation_jsonl),
        "evaluation_jsonl_sha256": sha256_file(args.evaluation_jsonl),
        "n_evaluation": len(rows),
        "species_counts": dict(sorted(counts.items())),
        "continuous_colour_representation": representation.label,
        "representation_dimension": representation.dimension,
        "input_was_prestandardized": representation.standardized,
        "scaler_provenance": scaler_provenance,
        "calibration_files": calibration_files,
        "species_scalers": scalers,
        "primary_k": PRIMARY_K,
        "sensitivity_k": list(SENSITIVITY_K),
        "n_permutations": args.permutations,
        "seed": args.seed,
        "graph_construction": "colour_blind_within_species_spherical_knn",
        "edge_statistic": "continuous_vector_colour_discontinuity",
        "permutation_unit": "complete_colour_vector_within_species",
        "global_weighting": "equal_species_mean",
        "threshold_tuning": False,
        "final_colour_classes_created": False,
        "complete_case_rows": len(rows),
        "dropped_rows": 0,
        "primary_result": primary,
        "sensitivity_same_direction": same_direction,
        "sensitivity_same_0_05_decision": sensitivity_significance,
        "output_sha256": {
            global_csv.name: sha256_file(global_csv),
            species_csv.name: sha256_file(species_csv),
            null_npz.name: sha256_file(null_npz),
        },
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    manifest_path = args.output_dir / "stage_a_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(primary, indent=2, sort_keys=True))
    print(f"Stage A outputs written to {args.output_dir}")


if __name__ == "__main__":
    main()
