#!/usr/bin/env python3
"""Freeze a source-only RGFCA balanced PHOTO robustness experiment.

The 42,111 taxon/93,209,351-photo dataset is only the opportunity universe.
Draw 1,000 taxa without looking at flower pixels or literature FCP labels.
Use precisely 200 x 250 species realizations: each of 1,000 taxa appears
50 times. Drawing 80 photos repeatedly from the SAME future ROI-screen pool
does not manufacture independent biological replicates.

This script never downloads, opens or classifies an original photograph.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

SOURCE_SHA = "ef26c16252dd785f2d08bde3b18938c7624d73a54386fa3d3974dceb09ed97a9"
FRAME_SHA = "5d871bd1f1190c68fd513bb527e6c93de1606dc35ff8cfff5f6187850dc382cc"
SOURCE_PHOTOS = 93_209_351
SOURCE_TAXA = 42_111
SALT = "FCP_RGFCA_BALANCED_PHOTO_POLYMORPHISM_20261016_V1"
STRATA = (("shallow_128_255", 128, 255, 334),
          ("medium_256_1023", 256, 1023, 333),
          ("deep_ge1024", 1024, None, 333))
N_SPECIES = 1000
REPS = 200
PER_REP = 250


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def order_hash(*values) -> str:
    return hashlib.sha256("|".join((SALT, *(str(x) for x in values))).encode()).hexdigest()


def source_frame(source_counts: Path, audit: Path) -> list[dict]:
    report = json.loads(audit.read_text(encoding="utf-8"))
    if (sha256(source_counts) != SOURCE_SHA
        or report.get("schema") != "fcp_42111_original_photo_opportunity_audit_v1"
        or report.get("status") != "VERIFIED_COMPLETE_COUNTS_ONLY_NOT_COLOUR_MEASUREMENT"
        or report.get("source_frame_species") != SOURCE_TAXA
        or report.get("source_frame_sha256") != FRAME_SHA
        or report.get("count_ledger_sha256") != SOURCE_SHA
        or report.get("all_source_photo_opportunity_ids") != SOURCE_PHOTOS
        or report.get("original_source_taxon_IDs_and_names_independently_matched") is not True
        or report.get("all_source_photo_id_manifests_frozen") is not False
        or report.get("all_original_image_pixels_measured") is not False):
        raise ValueError("Original 42111 count universe not independently audited")
    with source_counts.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames != ["inat_taxon_id", "species",
                                 "n_unique_source_photo_ids", "snapshot_date"]:
            raise ValueError("Original photo opportunity CSV schema changed")
        rows = list(reader)
    if len(rows) != SOURCE_TAXA:
        raise ValueError("Original taxon universe incomplete")
    seen = set()
    total = 0
    for r in rows:
        tid = r["inat_taxon_id"]
        n = r["n_unique_source_photo_ids"]
        if (not tid.isdecimal() or not n.isdecimal()
            or int(tid) <= 0 or int(tid) in seen
            or not r["species"].strip()
            or r["snapshot_date"] != "2026-09-27"):
            raise ValueError("Invalid or repeated source taxon/photo count")
        seen.add(int(tid))
        r["inat_taxon_id"] = int(tid)
        r["n_photo"] = int(n)
        total += int(n)
    if total != SOURCE_PHOTOS:
        raise ValueError("Source photo-count total changed")
    return rows


def choose_and_schedule(rows: list[dict]) -> tuple[list[dict], list[dict], dict]:
    cohort = []
    capacities = {}
    for label, low, high, quota in STRATA:
        pool = [r for r in rows if r["n_photo"] >= low
                and (high is None or r["n_photo"] <= high)]
        capacities[label] = len(pool)
        if len(pool) < quota:
            raise ValueError("Original opportunity stratum underfilled")
        pool.sort(key=lambda r: (order_hash("enrolment", r["inat_taxon_id"], r["species"]),
                                 r["inat_taxon_id"]))
        for p in pool[:quota]:
            cohort.append({**p, "tier": label})
    if len(cohort) != N_SPECIES or len({r["inat_taxon_id"] for r in cohort}) != N_SPECIES:
        raise ValueError("Original outcome-blind study cohort not 1,000 unique species")
    cohort.sort(key=lambda r: r["inat_taxon_id"])
    assignments = []
    totals = Counter()
    for cycle in range(50):
        randomized = sorted(cohort, key=lambda r: (
            order_hash("cycle", cycle, r["inat_taxon_id"]), r["inat_taxon_id"]))
        for j in range(4):
            realization = cycle * 4 + j + 1
            batch = randomized[j * PER_REP:(j + 1) * PER_REP]
            if len(batch) != PER_REP:
                raise ValueError("RGFCA realization does not contain exactly 250 taxa")
            for slot, r in enumerate(batch, 1):
                tid = r["inat_taxon_id"]
                assignments.append({"realization": realization, "cycle": cycle + 1,
                                    "slot": slot, "inat_taxon_id": tid,
                                    "species": r["species"],
                                    "source_opportunity_tier": r["tier"]})
                totals[tid] += 1
    if len(assignments) != REPS * PER_REP or set(totals.values()) != {50}:
        raise ValueError("Balanced RGFCA inclusion ledger invalid")
    checkpoints = {}
    for cutoff in (25, 50, 100, 150, 200):
        c = Counter(s["inat_taxon_id"] for s in assignments if s["realization"] <= cutoff)
        checkpoints[str(cutoff)] = {
            "n_species_seen": len(c),
            "min_inclusions": min(c.get(r["inat_taxon_id"], 0) for r in cohort),
            "max_inclusions": max(c.get(r["inat_taxon_id"], 0) for r in cohort)}
    return cohort, assignments, {"capacity_by_tier": capacities,
                                  "coverage_checkpoints": checkpoints}


def output(rows: list[dict], outdir: Path) -> dict:
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError("Do not overwrite frozen RGFCA realization schedule")
    cohort, assignment, diagnostics = choose_and_schedule(rows)
    outdir.mkdir(parents=True, exist_ok=True)
    cpath = outdir / "rgfca_species_1000_original_source_only.csv"
    apath = outdir / "rgfca_200x250_balanced_species_schedule.csv"
    with cpath.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["inat_taxon_id", "species", "source_photo_opportunities",
                         "opportunity_stratum", "enrolment_hash"])
        for r in cohort:
            writer.writerow([r["inat_taxon_id"], r["species"], r["n_photo"],
                             r["tier"], order_hash("enrolment", r["inat_taxon_id"], r["species"])])
    with apath.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(assignment[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(assignment)
    report = {
        "schema": "fcp_rgfca_photo_polymorphism_balanced_sample_v1",
        "status": "ORIGINAL_SOURCE_BALANCED_REPETITION_PLAN_NO_COLOUR_MEASURED",
        "original_source_snapshot_date": "2026-09-27",
        "source_frame_species": SOURCE_TAXA, "source_photo_opportunities": SOURCE_PHOTOS,
        "original_source_count_sha256": SOURCE_SHA, "original_frame_sha256": FRAME_SHA,
        "eligible_tier_counts": diagnostics["capacity_by_tier"],
        "eligible_source_taxa_ge128": sum(diagnostics["capacity_by_tier"].values()),
        "metadata_only_sampling_quotas": {label: quota for label, low, high, quota in STRATA},
        "selected_species": N_SPECIES,
        "realizations": REPS, "species_per_realization": PER_REP,
        "n_balanced_species_assignments": len(assignment),
        "species_inclusions_each": 50, "species_sampling_imbalance": 0,
        "primary_repeated_photo_draw_per_species": 120,
        "photo_draw_per_independent_observer_fold": 60,
        "photo_depth_sensitivities": [40, 80, 120, 128],
        "unique_photo_acquisition_pool_max_per_species": 128,
        "maximum_unique_original_image_pixels_if_pool_filled": N_SPECIES * 128,
        "repetition_resamples_original_photo_pool_not_independent_new_biological_data": True,
        "species_not_enrolled_are_NOT_monomorphic_controls": True,
        "sampled_1000_not_unbiased_all_flowering_plant_prevalence": True,
        "old_149900_flower_colours_used": False,
        "new_flower_images_opened": False,
        "photographic_poly_probabilities_measured": False,
        "convergence": diagnostics["coverage_checkpoints"],
        "selected_species_csv_sha256": sha256(cpath),
        "realization_schedule_csv_sha256": sha256(apath),
    }
    (outdir / "rgfca_sample_receipt.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-counts", type=Path, required=True)
    p.add_argument("--independent-audit", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    args = p.parse_args(argv)
    result = output(source_frame(args.source_counts, args.independent_audit), args.outdir)
    print("FCP_RGFCA_REAL_42111_PHOTO_SOURCE_BALANCED_SPECIES_SAMPLE",
          json.dumps(result, sort_keys=True))
    return result


if __name__ == "__main__":
    main()
