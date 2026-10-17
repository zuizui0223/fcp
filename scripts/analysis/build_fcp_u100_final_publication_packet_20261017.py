#!/usr/bin/env python3
"""Make a manuscript-ready FCP U100 result packet, only from authenticated REAL outputs.

Descriptive reporting only: never reclassifies a photograph, refits any model,
alters the 199-vertex null, or upgrades a post hoc 50-km radius. Fails closed
if a five-cohort full-pixel measurement and all final inferential files are
not present, mutually SHA-verified, and numerically consistent.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

COHORTS = ("discovery", "validation", "prospective", "wave4", "new_U100")
HISTORICAL = {
    "discovery": (166, 0.020529254583812922, 0.005),
    "validation": (181, 0.018672971642749295, 0.005),
    "prospective": (204, 0.01468491968437185, 0.005),
    "wave4": (190, 0.017293596294403593, 0.005),
}
ORIGINAL_SOURCE_SHAS = {
    "discovery": "ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4",
    "validation": "0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6",
    "prospective": "57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186",
}
EXCLUDED_TAXON_IDS = [62367, 75633, 1045161]


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def check(test: bool, reason: str) -> None:
    if not test:
        raise ValueError("FCP_U100_UNVERIFIED: " + reason)


def finite(x: object, name: str) -> float:
    check(type(x) in (int, float), f"{name}: nonnumeric or boolean")
    f = float(x)
    check(math.isfinite(f), f"{name}: nonfinite")
    return f


def source_hashes(paths: dict[str, Path], release: dict, combined: dict) -> dict:
    expected = {
        "photos": (release.get("source_archive_sha256"),
                   combined.get("merged_photo_table_sha256")),
        "combine": (release.get("source_join_receipt_sha256"),),
        "spatial": (release.get("original_spatial_result_sha256"),),
        "h1w": (release.get("original_H1_D_W_result_sha256"),),
    }
    found = {}
    for name, path in paths.items():
        check(path.is_file(), f"{name}: missing input")
        actual = digest(path)
        checks = expected[name]
        check(all(isinstance(v, str) and len(v) == 64 and v == actual for v in checks),
              f"{name}: SHA256 mismatch against full-run receipt")
        found[name] = actual
    return found


def spatial_row(row: dict, label: str) -> dict:
    check(isinstance(row, dict), f"{label}: not an object")
    n = row.get("n_species")
    check(type(n) is int and n >= 0, f"{label}: invalid eligible n")
    if n == 0:
        check(row.get("evaluable") is False, f"{label}: empty test not marked non-evaluable")
        return {"n_species": 0, "mean_depletion": None, "p_upper": None}
    check(row.get("evaluable") is True and finite(row.get("radius_km"), label + ".radius") == 50.0,
          f"{label}: invalid 50-km test")
    dep = finite(row.get("mean_depletion"), label + ".depletion")
    p = finite(row.get("mean_depletion_p_upper"), label + ".p")
    check(-1.0 <= dep <= 1.0, f"{label}: invalid diversity contrast")
    check(0.005 <= p <= 1.0 and abs(p * 200 - round(p * 200)) < 1e-7,
          f"{label}: not a plus-one p from exactly 199 nulls")
    for key in ("mean_D_pair", "mean_D_local"):
        value = finite(row.get(key), f"{label}.{key}")
        check(0.0 <= value <= 1.0, f"{label}: invalid pair diversity")
    check(abs(row["mean_D_pair"] - row["mean_D_local"] - dep) < 1e-9,
          f"{label}: mean pairwise depletion does not equal global minus local")
    return {"n_species": n, "mean_depletion": dep, "p_upper": p}


def verify_packet(receipt: dict, joined: dict, spatial: dict, h1w: dict,
                  paths: dict[str, Path]) -> dict:
    check(receipt.get("schema") == "fcp_all4730_permanent_measured_original_photo_results_v1",
          "not a genuine complete measured-results receipt")
    check(type(receipt.get("final_run_id")) is int and receipt["final_run_id"] > 0,
          "no authenticated full inferential workflow identity")
    check(receipt.get("full_original_eligible_species") == 4730
          and receipt.get("genuinely_measured_100photo_species") == 4475
          and receipt.get("genuinely_measured_photo_rows") == 447500
          and receipt.get("actual_100photo_species_analysis_omissions_total") == 255,
          "4,730 frame or 4,475 measured-species denominators invalid")
    check(receipt.get("old_149900_original_photos_unchanged") is True
          and receipt.get("original_149900_50km_three_effects_exactly_reproduced") is True
          and receipt.get("no_new_ecological_statistic_or_classifier") is True,
          "original method/measurement change")
    check(joined.get("schema") == "fcp_all4730_original_149900_plus_all_new_pixel_measurements_v1"
          and joined.get("total_original_100photo_measured_species") == 4475
          and joined.get("total_original_pixel_terminal_photo_rows") == 447500
          and joined.get("original_149900_photo_rows_unchanged") is True
          and joined.get("old_50km_and_W_estimators_unchanged") is True
          and joined.get("same_original_iNaturalist_source_not_independent_provider") is True,
          "five-cohort original source or method inconsistent")
    check(joined.get("new_U100_source_reuse_excluded_taxa") == 3
          and joined.get("new_U100_source_reuse_excluded_taxon_IDs") == EXCLUDED_TAXON_IDS
          and joined.get("new_U100_source_reuse_excluded_original_photo_rows") == 300
          and joined.get("remaining_U100_under_100_original_H9_photos") == 252,
          "three whole-taxon reuse exclusions or 252 under-100 taxa not applied")
    check(all(joined.get("source_sha256", {}).get(k) == v
              for k, v in ORIGINAL_SOURCE_SHAS.items()),
          "one of the frozen old 149,900 source inputs changed")
    cohorts = joined.get("cohorts")
    check(isinstance(cohorts, list) and len(cohorts) == 5
          and [c.get("cohort") for c in cohorts] == list(COHORTS),
          "missing or reordered five source cohorts")
    photo_by_cohort = {}
    for c in cohorts:
        name = c["cohort"]
        n = c.get("named_species")
        p = c.get("photo_rows")
        check(type(n) is int and type(p) is int and n >= 1 and p == n * 100,
              f"{name}: original per-species 100-photo denominator missing")
        photo_by_cohort[name] = {"species": n, "photo_rows": p,
                                  "classifiable_rows": c.get("four_state_classifiable_rows")}
    check(sum(c["species"] for c in photo_by_cohort.values()) == 4475
          and sum(c["photo_rows"] for c in photo_by_cohort.values()) == 447500
          and [photo_by_cohort[k]["species"] for k in COHORTS] == [500,500,499,498,2478],
          "cohort total or original independent new-U100 photo count incorrect")
    check(spatial.get("schema") == "fcp_all4730_same_original_50km_vertex_null_v1"
          and spatial.get("old_three_original_149900_effects_exactly_reproduced") is True
          and spatial.get("original_50km_radius_changed") is False
          and spatial.get("original_species_equal_weighting_changed") is False
          and spatial.get("original_199_vertex_null_changed") is False
          and spatial.get("still_posthoc_50km") is True,
          "historical spatial test/null not authenticated")
    all_sources = spatial.get("all_sources")
    check(isinstance(all_sources, dict) and set(all_sources) == set(COHORTS),
          "spatial test missing one of five cohorts")
    rows = {k: spatial_row(all_sources[k].get("50"), k) for k in COHORTS}
    for name, (n, dep, p) in HISTORICAL.items():
        measured = rows[name]
        check(measured["n_species"] == n
              and abs(measured["mean_depletion"] - dep) < 1e-10
              and measured["p_upper"] == p,
              f"{name}: previously frozen spatial result changed")
    pooled = spatial_row(spatial.get("pooled_original_50km"), "pooled")
    total_n = sum(x["n_species"] for x in rows.values())
    check(total_n == pooled["n_species"]
          and total_n == spatial.get("primary_50km_evaluable_species")
          and spatial.get("primary_50km_evaluable_species_by_cohort") ==
              {k: rows[k]["n_species"] for k in COHORTS},
          "pooled true local-pair eligible n differs from cohort denominators")
    check(pooled["n_species"] > 0, "no final spatially evaluable species")
    expected_effect = sum(x["n_species"] * x["mean_depletion"] for x in rows.values()
                          if x["n_species"] > 0) / pooled["n_species"]
    check(abs(pooled["mean_depletion"] - expected_effect) < 1e-9,
          "pooled species-equal effect inconsistent with five cohorts")
    check(h1w.get("schema") == "fcp_all4730_unchanged_historical_H1_D_H2_W_v1"
          and h1w.get("all_real_photo_rows") == 447500
          and h1w.get("existing_method_called_directly") is True
          and h1w.get("old_frozen_149900_inference_not_overwritten") is True
          and h1w.get("same_iNaturalist_source") is True
          and isinstance(h1w.get("H1_original_observer_disjoint"), dict)
          and bool(h1w["H1_original_observer_disjoint"])
          and isinstance(h1w.get("H2_W_original_fixed_white_axis"), dict)
          and bool(h1w["H2_W_original_fixed_white_axis"]),
          "H1/D/W did not complete with the genuine old method")
    source = source_hashes(paths, receipt, joined)
    return {
        "schema": "fcp_u100_completed_only_manuscript_packet_v1",
        "source_release_run_id": receipt["final_run_id"],
        "actual_measured_photos": 447500,
        "actual_measured_species": 4475,
        "source_frame_species": 4730,
        "missing_under_100_species": 252,
        "excluded_reused_photo_species": 3,
        "same_source_not_independent_provider": True,
        "radius_km": 50,
        "post_hoc_radius": True,
        "original_null_permutations": 199,
        "primary_spatial_species": pooled["n_species"],
        "cohorts": {k: {**photo_by_cohort[k], **rows[k]} for k in COHORTS},
        "pooled_50km": pooled,
        "pooled_50km_positive_under_null": (
            pooled["mean_depletion"] > 0 and pooled["p_upper"] < 0.05
        ),
        "observer_D_H1_and_fixed_white_axis_W_authentically_executed": True,
        "h1w_requires_separate_reporting_of_actual_full_numeric_results": True,
        "source_sha256": source,
        "source_inference_complete": True,
    }


def manuscript_text(packet: dict) -> str:
    def value(x: float | None, fmt: str) -> str:
        return format(x, fmt) if x is not None else "not estimable"
    lines = [
        "# FCP five-cohort manuscript-ready Results — completed authentic measurement only",
        "",
        f"Authenticated final GitHub Actions inferential run: **{packet['source_release_run_id']}**. "
        "This draft may be copied into a manuscript only alongside the linked permanent release receipt.",
        "",
        f"**Actual measurement:** **{packet['actual_measured_species']:,} distinct species / "
        f"{packet['actual_measured_photos']:,} original photograph terminals**, "
        f"out of an outcome-blind **{packet['source_frame_species']:,}-taxon** U100 source frame. "
        "**252** taxa lacked 100 source photographs and **three entire new taxa** were "
        "excluded for genuine cross-cohort reused image/observation identities. "
        "No source-eligible but unmeasured photographs have been counted.",
        "",
        "## Unchanged 50-km post hoc within-species colour-structure comparison",
        "",
        "| Cohort | Actually measured species | Photo terminals | 50-km eligible species | "
        "Global-minus-local pair dissimilarity | Plus-one upper p (199 nulls) |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for cohort, label in (
        ("discovery", "Discovery (original)"),
        ("validation", "Validation (original)"),
        ("prospective", "Prospective (original)"),
        ("wave4", "Fourth cohort (post-freeze)"),
        ("new_U100", "Fifth cohort (new U100)"),
    ):
        row = packet["cohorts"][cohort]
        lines.append(
            f"| {label} | {row['species']:,} | {row['photo_rows']:,} | "
            f"{row['n_species']:,} | {value(row['mean_depletion'], '.7f')} | "
            f"{value(row['p_upper'], '.3f')} |"
        )
    pooled = packet["pooled_50km"]
    lines.extend([
        f"| **All five cohorts (species equal weight)** | **4,475** | **447,500** | "
        f"**{pooled['n_species']:,}** | **{pooled['mean_depletion']:.7f}** | "
        f"**{pooled['p_upper']:.3f}** |",
        "",
        ("The pooled local-depletion contrast is positive under the unchanged vertex-label "
         "permutation null (one-sided plus-one p < 0.05)."
         if packet["pooled_50km_positive_under_null"]
         else "The pooled local-depletion contrast is **not supported** as positive under "
              "the unchanged vertex-label permutation null (one-sided plus-one p >= 0.05 "
              "or the observed contrast is not positive)."),
        "This inference conditions on each species' photographed colour composition and "
        "fixed coordinates, requires at least 40 classifiable photos and 30 local pairs, "
        "uses 199 within-species vertex permutations and retains the **post hoc 50-km** radius. "
        "The minimum possible upper p is 0.005.",
        "",
        "## Additional unchanged method outputs and limits",
        "",
        "The original observer-disjoint H1/D and structured-null fixed white-axis W analyses "
        "completed on the same 447,500 measured photographs. **Their detailed coefficients, "
        "null baselines and eligibility denominators must be reported from the archived "
        "H1/D/W result JSON; no H1/W support verdict is inferred by this packet.** "
        "Squared W cannot establish direction of evolutionary change toward white.",
        "",
        "These are species-conditional **photographic patterns** from the same iNaturalist "
        "source, not natural-population morph-frequency estimates or independent-provider "
        "replication. They do not identify pollinator-mediated selection, local adaptation, "
        "gene flow, drift or the evolutionary direction of flower colour.",
        "",
        "All inputs were individually SHA256 matched to the genuinely completed permanent "
        "five-cohort release receipt; the original 149,900 source photos and three frozen "
        "50-km cohort results were independently checked unchanged.",
        "",
    ])
    return "\n".join(lines)


def main() -> None:
    p = argparse.ArgumentParser()
    for name in ("receipt", "combine", "spatial", "h1w", "photos", "out"):
        p.add_argument("--" + name, required=True, type=Path)
    a = p.parse_args()
    paths = {"photos": a.photos, "combine": a.combine,
             "spatial": a.spatial, "h1w": a.h1w}
    for name, path in {**paths, "receipt": a.receipt}.items():
        check(path.is_file(), name + ": missing final source file")
    packet = verify_packet(
        json.loads(a.receipt.read_text(encoding="utf-8")),
        json.loads(a.combine.read_text(encoding="utf-8")),
        json.loads(a.spatial.read_text(encoding="utf-8")),
        json.loads(a.h1w.read_text(encoding="utf-8")),
        paths,
    )
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "fcp_all4730_verified_manuscript_results_packet.json").write_text(
        json.dumps(packet, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    (a.out / "fcp_all4730_verified_manuscript_results_packet.md").write_text(
        manuscript_text(packet), encoding="utf-8",
    )
    print("FCP_ALL4730_VERIFIED_PUBLICATION_PACKET_COMPLETE",
          packet["primary_spatial_species"],
          packet["pooled_50km"]["mean_depletion"],
          packet["pooled_50km"]["p_upper"])


if __name__ == "__main__":
    main()
