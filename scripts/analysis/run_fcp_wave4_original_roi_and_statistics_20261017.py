#!/usr/bin/env python3
"""Wave4 glue ONLY: original ROI-v4/EfficientSAM & frozen FCP statistics unchanged.

This module transports already-frozen functions from the original 149900-photo
study to additional outcome-blind, photo-ID-disjoint specimens. It does not
implement a new colour classifier or a modified ecological null hypothesis.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "analysis"))

from fcp_pipeline.photo_first_measurement_execution import (  # noqa: E402
    WORKER_FIELDS, ACQUISITION_FIELDS, compute_partition, semantic_shard,
    reassemble_complete_measurement,
)
from fcp_pipeline.third_cohort_prospective_execution_gate import (  # noqa: E402
    measurement_id, measurement_batch,
)

N_BATCHES = 2
N_SHARDS = 32
N_PARTS = 4
TOTAL_PARTS = N_BATCHES * N_SHARDS * N_PARTS
TARGET = 100
MORPHS = {"white", "yellow_orange", "red_pink", "blue_purple"}
BIO9 = ["white", "yellow", "orange", "red", "pink", "magenta", "purple", "blue", "bronze"]
FRACTIONS = ["flower_fraction_" + x for x in BIO9]
SOURCE_RUN = 38989666486


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def authenticated_source(source: Path) -> tuple[dict, dict, pd.DataFrame]:
    selection_path = source / "selection/selection_result.json"
    photo_path = source / "metadata/wave4_authorized_metadata.csv.gz"
    metadata_path = source / "metadata/metadata_result.json"
    for path in [selection_path, photo_path, metadata_path]:
        if not path.is_file():
            raise RuntimeError("wave4 original metadata source artifact incomplete: " + str(path))
    sel = json.loads(selection_path.read_text())
    meta = json.loads(metadata_path.read_text())
    if sel.get("wave4_selected_species") != 500 or sel.get("remaining_candidate_species") != 2730:
        raise RuntimeError("the original wave4 outcome-blind selection changed")
    if meta.get("selection_sha256") != sel.get("new_selection_manifest_sha256"):
        raise RuntimeError("wave4 selection manifest SHA mismatch")
    if meta.get("metadata_gate_pass") is not True or meta.get("image_pixels_opened") is not False:
        raise RuntimeError("unqualified original photo metadata must not open pixels")
    if meta.get("new_selected_species") != 500 or meta.get("original_photos") != 149900:
        raise RuntimeError("original study population changed")
    n = int(meta["new_exact_100_species"])
    if not 300 <= n <= 500 or int(meta["new_exact_100_photo_rows"]) != n * TARGET:
        raise RuntimeError("original 100-photo/spp support gate failed")
    if sha(photo_path) != meta["outputs"]["wave4_authorized_metadata.csv.gz"]:
        raise RuntimeError("new original photo metadata SHA drift")
    photo = pd.read_csv(photo_path, low_memory=False)
    required = {"photo_id", "observation_id", "photo_url_large", "photo_license", "inat_taxon_id",
                "latitude", "longitude", "query_species", "observer_id"}
    if required - set(photo.columns):
        raise RuntimeError("original metadata missing required columns: " + str(sorted(required - set(photo.columns))))
    if len(photo) != n * TARGET or photo.photo_id.nunique() != len(photo) or photo.observation_id.nunique() != len(photo):
        raise RuntimeError("wave4 original photograph photo-ID/observation-ID count mismatch")
    counts = photo.groupby("inat_taxon_id").size()
    if len(counts) != n or not counts.eq(TARGET).all():
        raise RuntimeError("wave4 original species are not exactly 100 original photos per species")
    if photo["photo_url_large"].isna().any() or photo["photo_license"].isna().any():
        raise RuntimeError("photo acquisition URL/license missing")
    if photo["query_species"].nunique() != n:
        raise RuntimeError("source photo taxonomy mismatch")
    return sel, meta, photo


def build_firewall(source: Path, dest: Path) -> dict:
    selected, source_receipt, photo = authenticated_source(source)
    dest.mkdir(parents=True, exist_ok=True)
    x = photo.reset_index(drop=True).copy()
    x.insert(0, "measurement_id", [measurement_id(p) for p in x["photo_id"]])
    x.insert(1, "measurement_batch", [measurement_batch(p) for p in x["measurement_id"]])
    if x.measurement_id.nunique() != len(x):
        raise RuntimeError("nonunique original frozen measurement ID")
    x["species"] = x["query_species"].astype(str)
    worker = pd.DataFrame({
        "measurement_id": x["measurement_id"].astype(str),
        "image_filename": x["measurement_id"].astype(str) + ".jpg",
        "photo_license": x["photo_license"].astype(str),
    })
    acq = pd.DataFrame({
        "measurement_id": worker["measurement_id"],
        "image_filename": worker["image_filename"],
        "photo_url_large": x["photo_url_large"].astype(str),
        "photo_license": worker["photo_license"],
    })
    if tuple(worker.columns) != WORKER_FIELDS or tuple(acq.columns) != ACQUISITION_FIELDS:
        raise RuntimeError("original three-field blind worker changed")
    (dest / "sealed_keys").mkdir(exist_ok=True)
    join = dest / "sealed_keys/metadata_join_key.csv"
    x.drop(columns=["photo_url_large"]).to_csv(join, index=False, lineterminator="\n")
    assignments = []
    for b in range(N_BATCHES):
        dirs = dest / f"batch_{b}"
        (dirs / "worker_packet").mkdir(parents=True, exist_ok=True)
        (dirs / "sealed_keys").mkdir(exist_ok=True)
        mask = x["measurement_batch"].eq(b)
        w = worker.loc[mask].reset_index(drop=True)
        a = acq.loc[mask].reset_index(drop=True)
        w.to_csv(dirs / "worker_packet/measurement_manifest.csv", index=False, lineterminator="\n")
        a.to_csv(dirs / "sealed_keys/acquisition_key.csv", index=False, lineterminator="\n")
        for m in w.measurement_id:
            assignments.append((m, b, semantic_shard(m, N_SHARDS), compute_partition(m, N_PARTS)))
    ass = pd.DataFrame(assignments, columns=["measurement_id", "batch", "semantic_shard", "compute_partition"])
    if len(ass) != len(x) or ass.measurement_id.nunique() != len(x):
        raise RuntimeError("blind partition assignment coverage incomplete")
    ass.to_csv(dest / "partition_assignments.csv", index=False, lineterminator="\n")
    record = {
        "schema": "fcp_same_149900_method_wave4_firewall_v1",
        "source_metadata_run_id": SOURCE_RUN,
        "original_frozen_results_changed": False,
        "original_photo_rows": 149900,
        "new_species": int(source_receipt["new_exact_100_species"]),
        "new_photo_rows": len(photo),
        "authorized_metadata_sha256": source_receipt["outputs"]["wave4_authorized_metadata.csv.gz"],
        "metadata_join_sha256": sha(join),
        "source_selection_sha256": selected["new_selection_manifest_sha256"],
        "partitions": TOTAL_PARTS,
        "blind_worker_columns": list(WORKER_FIELDS),
        "uses_original_third_cohort_blind_id_and_partition_functions": True,
        "source_photo_urls_in_worker": False,
        "source_coordinates_in_worker": False,
        "photo_pixels_opened": False,
        "original_inference_unmodified": True,
    }
    (dest / "measurement_firewall_manifest.json").write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    return record


def reassemble(firewall: Path, source_parts: Path, dest: Path) -> dict:
    record = json.loads((firewall / "measurement_firewall_manifest.json").read_text())
    if record.get("schema") != "fcp_same_149900_method_wave4_firewall_v1" or record.get("partitions") != TOTAL_PARTS:
        raise RuntimeError("new wave original ROI firewall not complete")
    workers = []
    for batch in range(N_BATCHES):
        workers.append(pd.read_csv(firewall / f"batch_{batch}/worker_packet/measurement_manifest.csv", dtype=str).fillna(""))
    worker = pd.concat(workers, ignore_index=True)
    join = pd.read_csv(firewall / "sealed_keys/metadata_join_key.csv", dtype={"measurement_id": str}, low_memory=False)
    if len(worker) != record["new_photo_rows"] or len(join) != len(worker) or sha(firewall / "sealed_keys/metadata_join_key.csv") != record["metadata_join_sha256"]:
        raise RuntimeError("new wave sealed metadata denominator changed")
    filenames = {}
    for p in source_parts.glob("b*_partition_s*_p*"):
        m = re.fullmatch(r"b(\d+)_partition_s(\d{2})_p(\d{2})\.(csv|json)", p.name)
        if m is None:
            raise RuntimeError("unknown terminal partition filename " + p.name)
        b, s, part, fmt = int(m[1]), int(m[2]), int(m[3]), m[4]
        if not (0 <= b < N_BATCHES and 0 <= s < N_SHARDS and 0 <= part < N_PARTS):
            raise RuntimeError("terminal partition index outside original matrix")
        key = (b, s, part, fmt)
        if key in filenames:
            raise RuntimeError("duplicate terminal source")
        filenames[key] = p
    expected = {(b, s, part) for b in range(N_BATCHES) for s in range(N_SHARDS) for part in range(N_PARTS)}
    for fmt in ("csv", "json"):
        if {k[:3] for k in filenames if k[3] == fmt} != expected:
            raise RuntimeError("exact 256 original worker terminal partitions incomplete")
    partition_data = []
    for b, s, part in sorted(expected):
        receipt = json.loads(filenames[(b, s, part, "json")].read_text())
        if (receipt.get("status") != "complete_random_photo_first_terminal_partition"
                or receipt.get("semantic_shard") != s or receipt.get("compute_partition") != part
                or receipt.get("source_urls_present") is not False
                or receipt.get("species_present") is not False
                or receipt.get("coordinates_present") is not False):
            raise RuntimeError("original sealed ROI worker terminal failed for " + str((b,s,part)))
        partition_data.append(pd.read_csv(filenames[(b,s,part,"csv")], dtype={"measurement_id": str}).fillna(""))
    assembled = reassemble_complete_measurement(
        partition_data, worker, join, expected_partition_receipts=TOTAL_PARTS
    )
    measured = assembled.joined_photos
    if len(measured) != record["new_photo_rows"] or measured.measurement_id.nunique() != len(measured):
        raise RuntimeError("photo pixel terminal left unmatched source photo IDs")
    measured["global_classifiable"] = measured["measurement_status"].eq("classified_four_state_morph") & measured["morph"].isin(MORPHS)
    missing = sorted(set(FRACTIONS) - set(measured.columns))
    if missing:
        raise RuntimeError("original nine-colour ROI fractions missing: " + str(missing))
    for col in FRACTIONS:
        measured[col] = pd.to_numeric(measured[col], errors="raise")
    palette = measured.loc[measured.global_classifiable, FRACTIONS].sum(axis=1)
    if len(palette) and not np.allclose(palette.to_numpy(float), 1.0, atol=1e-8):
        raise RuntimeError("original fixed palette normalization changed")
    by_species = measured.groupby(["species", "inat_taxon_id"]).agg(
        raw_photos=("measurement_id", "size"),
        classifiable_photos=("global_classifiable", "sum"),
    ).reset_index()
    if len(by_species) != record["new_species"] or not by_species.raw_photos.eq(TARGET).all():
        raise RuntimeError("new image terminal species 100-photo denominator lost")
    by_species["measurement_evaluable"] = by_species.classifiable_photos.ge(40)
    dest.mkdir(parents=True, exist_ok=True)
    csv_path = dest / "wave4_original_method_measured_photos.csv"
    support_path = dest / "wave4_species_measurement_support.csv"
    measured.to_csv(csv_path, index=False, lineterminator="\n")
    by_species.to_csv(support_path, index=False, lineterminator="\n")
    out = {
        "schema": "fcp_same_149900_method_wave4_real_photo_terminals_v1",
        "original_results_changed": False,
        "wave4_new_photo_ids": len(measured),
        "wave4_classifiable_photo_rows": int(measured.global_classifiable.sum()),
        "wave4_nonclassifiable_photo_rows": int((~measured.global_classifiable).sum()),
        "wave4_100photo_species": record["new_species"],
        "wave4_D_evaluable_species": int(by_species.measurement_evaluable.sum()),
        "terminal_partitions": TOTAL_PARTS,
        "measurement_original_fixed_colour_worker": True,
        "source_metadata_sha256": record["authorized_metadata_sha256"],
        "measured_sha256": sha(csv_path),
        "support_sha256": sha(support_path),
        "no_floral_colour_imputed": True,
        "same_iNaturalist_and_classifier_family": True,
    }
    (dest / "terminal_result.json").write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    return out


def original_spatial_primary_50km_n(summary: dict) -> int:
    """Read the unchanged original cohort_summary's string-keyed 50km result.

    cohort_summary() has always emitted radii as str(int(radius)). It never
    emits a numeric 50.0 key. This is reporting only, not an estimator change.
    """
    row = summary.get("50", {})
    if not row or row.get("evaluable") is False:
        return 0
    n = int(row["n_species"])
    if n < 0:
        raise RuntimeError("original spatial eligible species count invalid")
    return n


def analyze_original_methods(measured: Path, dest: Path) -> dict:
    """Call historical H1, 50km, and W functions without changing their mathematics."""
    from run_polymorphism_h1_observer_disjoint_reliability_20260913 import prepare_cohort, analyze_cohort
    from run_polymorphism_distributed_polymorphism_posthoc_20261007 import (
        analyse_species, cohort_summary, RADII_KM, as_bool
    )
    from run_polymorphism_delta_geometry_validation_20260912 import (
        load_classifiable, species_delta_vectors, nclass40_count
    )
    from run_polymorphism_h2_third_cohort_prospective_white_axis_20260917 import (
        white_contrast, run_threshold, SEEDS, THRESHOLDS
    )
    df = pd.read_csv(measured, low_memory=False)
    if df.photo_id.nunique() != len(df) or not df.groupby("inat_taxon_id").size().eq(TARGET).all():
        raise RuntimeError("wave4 original photo IDs/spp denominator changed")
    full_d, eligible_h1, profiles = prepare_cohort(measured, "wave4")
    partition, primary_h1, sensitivity_h1 = analyze_cohort(profiles, "wave4")
    work = df[as_bool(df.global_classifiable) & df.morph.isin(MORPHS)].copy()
    counts = work.groupby("inat_taxon_id").size()
    work = work.loc[work.inat_taxon_id.isin(counts[counts >= 40].index)].copy()
    by_radius = {r: [] for r in RADII_KM}
    metrics = []
    for _, g in work.groupby("inat_taxon_id", sort=True):
        g = g.dropna(subset=["latitude","longitude"])
        if len(g) < 40: continue
        row, null = analyse_species(g, "third")  # same historic seed/vertex label code
        metrics.append(row)
        for rad in RADII_KM:
            by_radius[rad].append(null[rad])
    spatial = cohort_summary(pd.DataFrame(metrics), by_radius) if metrics else {}
    original_h2_work = load_classifiable(measured)
    if nclass40_count(original_h2_work) != int(len(full_d)):
        raise RuntimeError("original H1 vs original W classifiable support mismatch")
    W = {}
    dest.mkdir(parents=True, exist_ok=True)
    for label, threshold in THRESHOLDS.items():
        selected, _, gates = species_delta_vectors(original_h2_work, threshold)
        estimate, null = run_threshold(original_h2_work, selected, white_contrast(), SEEDS[label])
        W[label] = {**estimate, "frozen_threshold": threshold, "old_gate_counts": gates}
        selected.to_csv(dest / (label + "_wave4_original_delta_vectors.csv"), index=False, lineterminator="\n")
        if null is not None:
            pd.DataFrame({"W": null}).to_csv(dest / (label + "_wave4_original_structured_null.csv"), index=False, lineterminator="\n")
    full_d.to_csv(dest / "wave4_species_D_original_method.csv", index=False, lineterminator="\n")
    pd.DataFrame(metrics).to_csv(dest / "wave4_species_50km_original_method.csv", index=False, lineterminator="\n")
    partition.to_csv(dest / "wave4_original_observer_disjoint_H1.csv", index=False, lineterminator="\n")
    out = {
        "schema": "fcp_same_149900_method_wave4_statistical_replication_v1",
        "original_frozen_results_unchanged": True,
        "wave4_measured_table_sha256": sha(measured),
        "wave4_D_evaluable_species": len(full_d),
        "wave4_50km_evaluable_species": original_spatial_primary_50km_n(spatial),
        "h1_original_method": {"primary": primary_h1, "sensitivity": sensitivity_h1},
        "spatial_original_method": {str(int(k)): v for k,v in spatial.items()},
        "h2_original_W_method": W,
        "methods": "historical original implementation imported verbatim; no alternative spatial weights, modified classifier or new null",
        "same_source_not_independent_provider": True,
        "post_hoc_50km_remains_post_hoc": True,
    }
    (dest / "wave4_original_149900_method_analysis_result.json").write_text(
        json.dumps(out, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="action", required=True)
    f = sub.add_parser("firewall")
    f.add_argument("--source", type=Path, required=True)
    f.add_argument("--out", type=Path, required=True)
    re = sub.add_parser("reassemble")
    re.add_argument("--firewall", type=Path, required=True)
    re.add_argument("--parts", type=Path, required=True)
    re.add_argument("--out", type=Path, required=True)
    inr = sub.add_parser("infer")
    inr.add_argument("--measured", type=Path, required=True)
    inr.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    if a.action == "firewall":
        result = build_firewall(a.source, a.out)
    elif a.action == "reassemble":
        result = reassemble(a.firewall, a.parts, a.out)
    else:
        result = analyze_original_methods(a.measured, a.out)
    print("FCP_WAVE4_ORIGINAL_METHOD_" + a.action.upper() + "_COMPLETE")
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
