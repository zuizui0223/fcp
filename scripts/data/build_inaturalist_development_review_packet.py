#!/usr/bin/env python3
"""Build the observation-level, morph-blind iNaturalist development review packet.

The upstream development sample stores one representative image per selected
observation.  Human flower-visibility and morph review, however, is performed
at the observation grain and must expose every reusable photo belonging to that
observation.  This script expands only the preassigned ``development_40``
partition, keeps spatial and observer metadata outside the reviewer sheets, and
profiles technical image properties without turning them into biological
labels or post-hoc exclusion thresholds.

This remains a moving-API workflow-development artifact.  It is not the frozen
Open Dataset pilot and cannot support morph-frequency or spatial inference.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import urllib.parse
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

from PIL import Image, ImageFilter, ImageStat

if __package__:
    from .build_inaturalist_photo_development_sample import (
        REUSABLE_LICENSES,
        medium_url,
        request_bytes,
        sha256,
        stable_rank,
    )
else:
    from build_inaturalist_photo_development_sample import (  # type: ignore[no-redef]
        REUSABLE_LICENSES,
        medium_url,
        request_bytes,
        sha256,
        stable_rank,
    )


DEFAULT_SEED = "fcp-inaturalist-development-review-v1"
REVIEW_FIELDS = [
    "review_order",
    "canonical_name",
    "encounter_blind_id",
    "image_files",
    "n_photos",
    "flower_visible",
    "single_encounter_usable",
    "anonymous_morph_code",
    "multiple_morphs_same_frame",
    "classification_confidence",
    "exclusion_reason",
    "reviewer_id",
    "reviewed_at_utc",
]
PROHIBITED_REVIEW_FIELDS = {
    "observation_id",
    "photo_id",
    "observer_id",
    "observed_on",
    "flowering_week",
    "latitude",
    "longitude",
    "grid_x_50km_epsg6933",
    "grid_y_50km_epsg6933",
    "photo_medium_url",
    "photo_attribution",
}
REVIEW_OUTCOME_FIELDS = {
    "flower_visible",
    "single_encounter_usable",
    "anonymous_morph_code",
    "multiple_morphs_same_frame",
    "classification_confidence",
    "exclusion_reason",
    "reviewer_id",
    "reviewed_at_utc",
}
CODEBOOK_OUTCOME_FIELDS = {
    "diagnostic_flower_region",
    "retained_morph_codes",
    "disallowed_lighting_or_age_cues",
    "multiple_morph_rule",
    "not_classifiable_rule",
    "literature_basis",
    "codebook_version",
    "frozen_at_utc",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def assert_no_human_outcomes_before_rebuild(output: Path) -> None:
    """Protect reviewer work and frozen codebooks from an accidental rebuild."""

    for name in ["reviewer_A_annotation_sheet.csv", "reviewer_B_annotation_sheet.csv"]:
        path = output / name
        if not path.exists():
            continue
        rows = read_csv(path)
        if any(str(row.get(field, "")).strip() for row in rows for field in REVIEW_OUTCOME_FIELDS):
            raise RuntimeError(f"refusing to overwrite completed human outcomes in {path}")
    codebook_path = output / "species_codebook_template.csv"
    if codebook_path.exists():
        rows = read_csv(codebook_path)
        if any(str(row.get(field, "")).strip() for row in rows for field in CODEBOOK_OUTCOME_FIELDS):
            raise RuntimeError(f"refusing to overwrite a populated species codebook in {codebook_path}")


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str] | None = None) -> None:
    if not rows:
        raise ValueError(f"cannot write empty table: {path}")
    fieldnames = fields or list(rows[0])
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def load_observation_details(raw_dir: Path) -> dict[int, dict[str, Any]]:
    details: dict[int, dict[str, Any]] = {}
    for path in sorted(raw_dir.glob("selected_*.json")):
        decoded = json.loads(path.read_text(encoding="utf-8"))
        for observation in decoded.get("results", []):
            observation_id = int(observation["id"])
            previous = details.get(observation_id)
            if previous is not None and previous != observation:
                raise ValueError(f"conflicting cached detail for observation {observation_id}")
            details[observation_id] = observation
    return details


def reusable_photos(observation: dict[str, Any]) -> list[dict[str, Any]]:
    photos = [
        photo
        for photo in observation.get("photos", [])
        if str(photo.get("license_code", "")).casefold() in REUSABLE_LICENSES
    ]
    return sorted(photos, key=lambda photo: int(photo["id"]))


def histogram_quantile(histogram: list[int], probability: float) -> int:
    if not 0 <= probability <= 1:
        raise ValueError("probability must be between zero and one")
    total = sum(histogram)
    if total <= 0:
        raise ValueError("empty histogram")
    threshold = probability * max(total - 1, 0)
    cumulative = 0
    for value, count in enumerate(histogram):
        cumulative += count
        if cumulative > threshold:
            return value
    return 255


def technical_profile(path: Path) -> dict[str, Any]:
    """Return descriptive technical metrics; no metric is an exclusion label."""

    with Image.open(path) as source:
        source.verify()
    with Image.open(path) as source:
        rgb = source.convert("RGB")
        width, height = rgb.size
        gray = rgb.convert("L")
        histogram = gray.histogram()
        pixels = width * height
        edge = gray.filter(ImageFilter.FIND_EDGES)
        if width > 2 and height > 2:
            edge = edge.crop((1, 1, width - 1, height - 1))
        edge_variance = float(ImageStat.Stat(edge).var[0])
        mean_rgb = [float(value) for value in ImageStat.Stat(rgb).mean]
        return {
            "image_bytes": path.stat().st_size,
            "image_width": width,
            "image_height": height,
            "aspect_ratio": width / height,
            "luminance_mean": sum(i * count for i, count in enumerate(histogram)) / pixels,
            "luminance_p01": histogram_quantile(histogram, 0.01),
            "luminance_p99": histogram_quantile(histogram, 0.99),
            "fraction_luminance_le_5": sum(histogram[:6]) / pixels,
            "fraction_luminance_ge_250": sum(histogram[250:]) / pixels,
            "mean_red": mean_rgb[0],
            "mean_green": mean_rgb[1],
            "mean_blue": mean_rgb[2],
            "edge_variance_descriptive": edge_variance,
        }


def file_extension(photo: dict[str, Any]) -> str:
    suffix = Path(urllib.parse.urlparse(str(photo["url"])).path).suffix.lower()
    return suffix if suffix in {".jpg", ".jpeg", ".png"} else ".jpg"


def materialize_photo(url: str, target: Path, reusable_source: Path | None) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        technical_profile(target)
        return
    partial = target.with_suffix(target.suffix + ".partial")
    if reusable_source is not None and reusable_source.exists():
        shutil.copy2(reusable_source, partial)
    else:
        payload, _headers = request_bytes(url)
        partial.write_bytes(payload)
    technical_profile(partial)
    partial.replace(target)


def reviewer_rows(
    encounters: Iterable[dict[str, Any]], reviewer: str, seed: str
) -> list[dict[str, Any]]:
    ranked = sorted(
        encounters,
        key=lambda row: stable_rank(
            seed, reviewer, row["canonical_name"], row["encounter_blind_id"]
        ),
    )
    rows: list[dict[str, Any]] = []
    for order, encounter in enumerate(ranked, start=1):
        row = {
            "review_order": order,
            "canonical_name": encounter["canonical_name"],
            "encounter_blind_id": encounter["encounter_blind_id"],
            "image_files": "|".join(encounter["image_files"]),
            "n_photos": len(encounter["image_files"]),
            "flower_visible": "",
            "single_encounter_usable": "",
            "anonymous_morph_code": "",
            "multiple_morphs_same_frame": "",
            "classification_confidence": "",
            "exclusion_reason": "",
            "reviewer_id": "",
            "reviewed_at_utc": "",
        }
        if PROHIBITED_REVIEW_FIELDS.intersection(row):
            raise AssertionError("reviewer row contains prohibited provenance fields")
        rows.append(row)
    return rows


def quantile(values: list[float], probability: float) -> float:
    if not values:
        raise ValueError("cannot summarize an empty vector")
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-artifact", type=Path, required=True)
    parser.add_argument("--output-artifact", type=Path, required=True)
    parser.add_argument(
        "--public-manifest",
        type=Path,
        default=Path(
            "docs/supporting/jbi_inaturalist_development_review_packet_manifest_v1.json"
        ),
    )
    parser.add_argument("--seed", default=DEFAULT_SEED)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    source = args.source_artifact.resolve()
    output = args.output_artifact.resolve()
    output.mkdir(parents=True, exist_ok=True)
    assert_no_human_outcomes_before_rebuild(output)
    images_out = output / "images"

    provenance = read_csv(source / "private_provenance.csv")
    development = [row for row in provenance if row["annotation_partition"] == "development_40"]
    details = load_observation_details(source / "raw_api")
    if len(development) != 480:
        raise RuntimeError(f"expected 480 development encounters, found {len(development)}")
    if len({row["blind_id"] for row in development}) != len(development):
        raise RuntimeError("development encounter blind IDs are not unique")

    encounters: list[dict[str, Any]] = []
    photo_ledger: list[dict[str, Any]] = []
    technical_rows: list[dict[str, Any]] = []
    species_stats: dict[str, Counter[str]] = defaultdict(Counter)
    acquired_at = datetime.now(timezone.utc).isoformat()

    for row in sorted(development, key=lambda value: value["blind_id"]):
        observation_id = int(row["observation_id"])
        observation = details.get(observation_id)
        if observation is None:
            raise RuntimeError(f"missing cached detail for development observation {observation_id}")
        photos = reusable_photos(observation)
        if not photos:
            raise RuntimeError(f"development observation {observation_id} has no reusable photo")
        representative_id = int(row["photo_id"])
        if representative_id not in {int(photo["id"]) for photo in photos}:
            raise RuntimeError(
                f"representative photo is absent from reusable photo set for {observation_id}"
            )

        image_files: list[str] = []
        for index, photo in enumerate(photos, start=1):
            photo_blind_id = f"{row['blind_id']}-P{index:02d}"
            relative = Path("images") / row["blind_id"] / (
                photo_blind_id + file_extension(photo)
            )
            target = output / relative
            existing = None
            if int(photo["id"]) == representative_id:
                existing = source / "images" / row["image_file"]
            url = medium_url(photo)
            materialize_photo(url, target, existing)
            profile = technical_profile(target)
            image_files.append(relative.as_posix())
            photo_ledger.append(
                {
                    "encounter_blind_id": row["blind_id"],
                    "photo_blind_id": photo_blind_id,
                    "canonical_name": row["canonical_name"],
                    "observation_id": observation_id,
                    "photo_id": int(photo["id"]),
                    "photo_license": str(photo.get("license_code", "")).casefold(),
                    "photo_attribution": photo.get("attribution", ""),
                    "photo_medium_url": url,
                    "image_file": relative.as_posix(),
                    "image_sha256": sha256(target),
                    "acquired_at_utc": acquired_at,
                }
            )
            technical_rows.append(
                {
                    "encounter_blind_id": row["blind_id"],
                    "photo_blind_id": photo_blind_id,
                    "canonical_name": row["canonical_name"],
                    "image_file": relative.as_posix(),
                    **profile,
                    "technical_status": "profile_complete_no_biological_admission",
                }
            )

        encounters.append(
            {
                "canonical_name": row["canonical_name"],
                "encounter_blind_id": row["blind_id"],
                "image_files": image_files,
            }
        )
        counter = species_stats[row["canonical_name"]]
        counter["development_encounters"] += 1
        counter["reusable_photos"] += len(photos)
        counter["multi_photo_encounters"] += int(len(photos) > 1)

    reviewer_a = reviewer_rows(encounters, "reviewer_A", args.seed)
    reviewer_b = reviewer_rows(encounters, "reviewer_B", args.seed)
    if [row["encounter_blind_id"] for row in reviewer_a] == [
        row["encounter_blind_id"] for row in reviewer_b
    ]:
        raise AssertionError("independent reviewer orderings unexpectedly match")

    photo_ledger_path = output / "private_photo_provenance.csv"
    technical_path = output / "technical_image_profile.csv"
    reviewer_a_path = output / "reviewer_A_annotation_sheet.csv"
    reviewer_b_path = output / "reviewer_B_annotation_sheet.csv"
    codebook_path = output / "species_codebook_template.csv"
    write_csv(photo_ledger_path, photo_ledger)
    write_csv(technical_path, technical_rows)
    write_csv(reviewer_a_path, reviewer_a, REVIEW_FIELDS)
    write_csv(reviewer_b_path, reviewer_b, REVIEW_FIELDS)
    write_csv(
        codebook_path,
        [
            {
                "canonical_name": species,
                "diagnostic_flower_region": "",
                "retained_morph_codes": "",
                "disallowed_lighting_or_age_cues": "",
                "multiple_morph_rule": "",
                "not_classifiable_rule": "",
                "literature_basis": "",
                "codebook_version": "",
                "frozen_at_utc": "",
            }
            for species in sorted(species_stats)
        ],
    )

    public_species: list[dict[str, Any]] = []
    for species in sorted(species_stats):
        counter = species_stats[species]
        profiles = [row for row in technical_rows if row["canonical_name"] == species]
        public_species.append(
            {
                "canonical_name": species,
                "development_encounters": counter["development_encounters"],
                "reusable_photos": counter["reusable_photos"],
                "multi_photo_encounters": counter["multi_photo_encounters"],
                "photos_per_encounter_mean": counter["reusable_photos"]
                / counter["development_encounters"],
                "image_width_min": min(row["image_width"] for row in profiles),
                "image_height_min": min(row["image_height"] for row in profiles),
                "luminance_p01_median": quantile(
                    [float(row["luminance_p01"]) for row in profiles], 0.5
                ),
                "luminance_p99_median": quantile(
                    [float(row["luminance_p99"]) for row in profiles], 0.5
                ),
            }
        )

    artifact_manifest = {
        "status": "complete_development_review_packet_pending_two_human_reviews",
        "protocol": "fcp-inaturalist-development-review-v1",
        "created_at_utc": acquired_at,
        "source_artifact": str(source),
        "source_artifact_manifest_sha256": sha256(source / "artifact_manifest.json"),
        "annotation_partition": "development_40_only",
        "review_grain": "one observation with all reusable CC0/CC-BY photos shown together",
        "seed": args.seed,
        "development_encounters": len(encounters),
        "development_photos": len(photo_ledger),
        "multi_photo_encounters": sum(
            counter["multi_photo_encounters"] for counter in species_stats.values()
        ),
        "admission_counts": {
            "review_ready_pending_two_human_reviews": len(encounters),
            "human_reviewed_encounters": 0,
            "biologically_admitted_encounters": 0,
            "human_excluded_encounters": 0,
            "anonymous_colour_marks_admitted": 0,
        },
        "species_summary": public_species,
        "technical_metrics_use": (
            "descriptive data-quality diagnostics only; no flower visibility, morph, "
            "or biological exclusion is assigned automatically"
        ),
        "private_files_sha256": {
            path.name: sha256(path)
            for path in [
                photo_ledger_path,
                technical_path,
                reviewer_a_path,
                reviewer_b_path,
                codebook_path,
            ]
        },
        "claim_ceiling": (
            "Observation-grain review readiness and technical image availability only. "
            "Zero flower-visibility or morph labels have been admitted; no morph frequency, "
            "spatial randomness, boundary concordance, or organizer effect may be estimated."
        ),
    }
    artifact_manifest_path = output / "artifact_manifest.json"
    artifact_manifest_path.write_text(
        json.dumps(artifact_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    public_manifest = {
        key: value
        for key, value in artifact_manifest.items()
        if key not in {"private_files_sha256", "source_artifact"}
    }
    public_manifest["source_artifact_reference"] = (
        "private moving-API development artifact documented by "
        "source_artifact_manifest_sha256"
    )
    public_manifest["artifact_manifest_sha256"] = sha256(artifact_manifest_path)
    public_manifest["public_data_exclusions"] = [
        "observation and photo identifiers",
        "exact coordinates and dates",
        "observer identifiers",
        "image URLs, files, and attribution ledger",
        "reviewer sheets before adjudication",
    ]
    args.public_manifest.parent.mkdir(parents=True, exist_ok=True)
    args.public_manifest.write_text(
        json.dumps(public_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(public_manifest, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
