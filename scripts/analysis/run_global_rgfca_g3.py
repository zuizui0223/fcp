#!/usr/bin/env python3
"""Run the frozen observed-only G3 prevalence and heterogeneity analysis."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from fcp_pipeline.global_g3 import run_g3_prevalence
from fcp_pipeline.global_repeated_atlas import schedule_audit
from fcp_pipeline.global_rgfca_engine import COLOUR_COLUMNS

ROOT = Path(__file__).resolve().parents[2]
EXECUTION = ROOT / "docs/supporting/global_monte_carlo_inference_execution_contract_v1.json"
MEASUREMENT = ROOT / "docs/supporting/global_monte_carlo_measurement_result_v1.json"
MEASURED = ROOT / "data/derived/global_monte_carlo_measured_photos_v1.csv"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    execution = json.loads(EXECUTION.read_text(encoding="utf-8"))
    measurement = json.loads(MEASUREMENT.read_text(encoding="utf-8"))
    gate = execution["input_gate"]
    if measurement.get("status") != gate["measurement_status_required"]:
        raise RuntimeError("global measurement status does not authorize G3")
    post = measurement.get("postmeasurement_gate", {})
    if post.get("pass") is not True or post.get("decision") != gate["measurement_postgate_decision_required"]:
        raise RuntimeError("global measurement postgate did not authorize G3")
    if sha256_file(MEASURED) != str(measurement.get("lineage", {}).get("measured_table_sha256") or ""):
        raise RuntimeError("measured table lineage differs from frozen measurement result")

    frame = pd.read_csv(MEASURED)
    required = {"photo_id", "species", "latitude", "longitude", "global_classifiable", *COLOUR_COLUMNS}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise RuntimeError(f"measured table lacks G3 inputs: {missing}")
    classifiable = frame["global_classifiable"].astype(str).str.casefold().isin({"true", "1"})
    pool = frame.loc[classifiable, ["photo_id", "species", "latitude", "longitude", *COLOUR_COLUMNS]].copy()
    minimum_pool = int(gate["minimum_classifiable_photos_per_species"])
    counts = pool.groupby("species", observed=True).size()
    eligible = set(counts[counts >= minimum_pool].index.astype(str))
    pool = pool.loc[pool["species"].astype(str).isin(eligible)].reset_index(drop=True)
    if int(pool["species"].nunique()) != int(post["evaluable_species"]):
        raise RuntimeError("G3 pool species count differs from measurement postgate")

    outer = execution["outer_schedule"]
    g3 = execution["g3_prevalence_heterogeneity"]
    variance_floor = float(g3["hierarchical_heterogeneity"]["sampling_variance_floor"])
    result = run_g3_prevalence(
        pool,
        n_outer=int(outer["observed_resamples"]),
        species_per_outer=int(outer["species_per_resample"]),
        photos_per_species=int(outer["photos_per_species"]),
        minimum_pool_photos_per_species=minimum_pool,
        species_seed=int(outer["species_seed"]),
        photo_master_seed=int(outer["photo_master_seed"]),
        variance_floor=variance_floor,
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    outer_path = args.output_dir / "global_rgfca_g3_outer_v1.csv"
    species_path = args.output_dir / "global_rgfca_g3_species_v1.csv"
    manifest_path = args.output_dir / "global_rgfca_g3_result_v1.json"
    result.outer.to_csv(outer_path, index=False, lineterminator="\n")
    result.species.to_csv(species_path, index=False, lineterminator="\n")
    manifest = {
        "protocol": execution["protocol"],
        "status": "complete_global_rgfca_g3_prevalence_heterogeneity",
        "confirmatory_p_value": None,
        "cannot_rescue_null_g1": True,
        "classifiable_pool_rows": int(len(pool)),
        "eligible_species": int(pool["species"].nunique()),
        "median_outer_mean_rho": result.median_outer_mean_rho,
        "median_outer_median_rho": result.median_outer_median_rho,
        "median_outer_positive_fraction": result.median_outer_positive_fraction,
        "hierarchical_heterogeneity": {
            "tau2_fisher_z": result.tau2_fisher_z,
            "species_used": result.tau2_species_used,
            "sampling_variance_floor": variance_floor,
            "method": "DerSimonian-Laird on species mean Fisher-z effects using repeated-photo resampling variance",
        },
        "schedule_audit": schedule_audit(result.schedule),
        "lineage": {
            "execution_contract_sha256": sha256_file(EXECUTION),
            "measurement_manifest_sha256": sha256_file(MEASUREMENT),
            "measured_table_sha256": sha256_file(MEASURED),
            "outer_csv_sha256": sha256_file(outer_path),
            "species_csv_sha256": sha256_file(species_path),
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
