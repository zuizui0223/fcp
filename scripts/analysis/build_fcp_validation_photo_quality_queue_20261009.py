#!/usr/bin/env python3
"""Frozen photo-quality review queue from source-verified 2026-10-08 metadata.

No pixels, photo URLs, colour classes or new iNaturalist requests are used.
The historical 10-km source-only feasibility HOLD remains unchanged.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

INPUT_SHA = {
    "gap1": "a64959733c77e5c6775f2e6a6568f24a543d7a74bf413d12f924d149e4ed653e",
    "gap2": "7ecb805d5fa7066d12f9f2964aa9acf72ede90d4a58fa1fdfdb3e0902f65390d",
    "status": "b808db25310a74da5968e19d19a2ef44db9005426dd2302536cfca27c75b6371",
    "queue": "a69df5ad49cc8cabde2217fdd7e792b57e9ce972abebefd37b044b3fd5e60bb0",
}
LICENSES = {"cc0", "cc-by", "cc-by-sa", "cc-by-nc", "cc-by-nc-sa"}
QUALITY_FIELDS = ("flower_visible", "species_identity_confirmed", "target_flower_organ",
                  "exposure_acceptable", "colour_measurable", "individual_identity_verified")
COLUMNS = ("gap_class", "inat_taxon_id", "species", "target_year", "calendar_month",
           "source_anchor_photo_id", "required_distinct_observers_in_year",
           "observation_id", "photo_id", "observer_id", "public_distance_km",
           "position_accuracy_m", "photo_license", *QUALITY_FIELDS, "quality_review_status")


def checked(path: Path, key: str):
    if hashlib.sha256(path.read_bytes()).hexdigest() != INPUT_SHA[key]:
        raise ValueError(f"Source SHA mismatch: {key}")
    if key == "queue":
        with path.open(newline="", encoding="utf-8") as f:
            return list(csv.DictReader(f))
    return json.loads(path.read_text(encoding="utf-8"))


def assemble(gap1: list, gap2: list, statuses: list, queue: list) -> dict:
    original = {}
    for r in queue:
        if r["cohort"] != "validation":
            continue
        key = (str(r["inat_taxon_id"]), int(r["target_year"]))
        if key in original:
            raise ValueError("Duplicate frozen taxon-year target")
        original[key] = r
    good1 = {str(r["inat_taxon_id"]): r for r in gap1
             if r["query_status"] == "ELIGIBLE_PUBLIC_METADATA_CANDIDATE_EXISTS"}
    good2 = {str(r["inat_taxon_id"]) for r in statuses
             if r["status"] == "POSSIBLE_COMPLETE_METADATA_ONLY"}
    if len(gap1) != 22 or len(statuses) != 60 or len(good1) != 15 or len(good2) != 10:
        raise ValueError("Unexpected historic 22/60/15/10 taxon counts")
    if set(good1) & good2:
        raise ValueError("Gap-one/two species overlap")
    by2 = defaultdict(list)
    for row in gap2:
        by2[str(row["inat_taxon_id"])].append(row)

    photos = []
    species = []
    seen_obs = set()
    seen_photo = set()
    for gap, taxa in ((1, sorted(good1, key=int)), (2, sorted(good2, key=int))):
        for tid in taxa:
            targets = [good1[tid]] if gap == 1 else by2[tid]
            if not targets:
                raise ValueError("Successful species without year metadata")
            requirements = {}
            months, anchors = set(), set()
            for cell in sorted(targets, key=lambda z: int(z["target_year"])):
                year = int(cell["target_year"])
                meta = original.get((tid, year))
                if meta is None or int(meta["priority_tier"]) != gap:
                    raise ValueError("Photograph year not in frozen opportunity")
                month = int(meta["calendar_month"])
                anchor = str(meta["source_anchor_photo_id"])
                if int(cell["month"]) != month:
                    raise ValueError("Calendar month mismatch")
                supplied_anchor = cell["original_anchor_photo_id"] if gap == 1 else cell["source_anchor_photo_id"]
                if str(supplied_anchor) != anchor:
                    raise ValueError("Photo-centred site anchor mismatch")
                months.add(month); anchors.add(anchor)
                required = int(meta["required_additional_distinct_observer_photos"])
                candidates = cell["eligible_photo_observation_ids"] if gap == 1 else cell["metadata_candidates"]
                observers = set()
                for candidate in candidates:
                    photo_id = str(candidate["photo_id"])
                    obs_id = str(candidate["observation_id"])
                    observer = str(candidate["observer_id"])
                    if not all(z.isdecimal() for z in (photo_id, obs_id, observer)):
                        raise ValueError("Malformed metadata ID")
                    if photo_id in seen_photo or obs_id in seen_obs or observer in observers:
                        raise ValueError("Duplicate candidate observer, image or observation")
                    d = float(candidate["public_distance_km"])
                    accuracy = float(candidate["position_accuracy_m"])
                    if not (0 <= d <= 5.000001 and 0 <= accuracy <= 5000):
                        raise ValueError("Unqualified site distance or accuracy")
                    if candidate["photo_license"] not in LICENSES:
                        raise ValueError("Unqualified photo licence")
                    seen_photo.add(photo_id); seen_obs.add(obs_id); observers.add(observer)
                    p = {
                        "gap_class": gap, "inat_taxon_id": tid,
                        "species": meta["species"], "target_year": year,
                        "calendar_month": month, "source_anchor_photo_id": anchor,
                        "required_distinct_observers_in_year": required,
                        "observation_id": obs_id, "photo_id": photo_id,
                        "observer_id": observer, "public_distance_km": d,
                        "position_accuracy_m": accuracy,
                        "photo_license": candidate["photo_license"],
                        "quality_review_status": "UNREVIEWED",
                    }
                    p.update({name: None for name in QUALITY_FIELDS})
                    photos.append(p)
                requirements[str(year)] = {
                    "required": required, "metadata_candidates": len(candidates)
                }
                if len(observers) < required:
                    raise ValueError("Insufficient year-specific observer candidates")
            if len(months) != 1 or len(anchors) != 1:
                raise ValueError("Same species changes month or site")
            if sum(v["required"] for v in requirements.values()) != gap:
                raise ValueError("One or two gap slots not retained")
            if gap == 2 and len(requirements) not in (1, 2):
                raise ValueError("Gap-two has wrong number of year cells")
            species.append({
                "inat_taxon_id": tid, "species": original[(tid, int(targets[0]["target_year"]))]["species"],
                "gap_class": gap, "calendar_month": next(iter(months)),
                "source_anchor_photo_id": next(iter(anchors)),
                "year_requirements": requirements,
                "species_quality_pass": None, "photo_review_completed": False
            })
    photos.sort(key=lambda x: (x["gap_class"], int(x["inat_taxon_id"]), x["target_year"], int(x["photo_id"])))
    species.sort(key=lambda x: (x["gap_class"], int(x["inat_taxon_id"])))
    if (len(species), len(photos)) != (25, 115):
        raise ValueError("Source-ledger taxon/photo count drift")
    if Counter(p["gap_class"] for p in photos) != {1: 39, 2: 76}:
        raise ValueError("Unbalanced historical 39+76 photo source ledger")
    return {
        "schema": "fcp_validation_frozen_image_quality_queue_v1",
        "status": "UNREVIEWED_METADATA_NO_PIXELS",
        "date_jst": "2026-10-09",
        "source_runs": [37798611129, 37801452055],
        "source_sha256": INPUT_SHA,
        "n_candidate_species": 25, "n_photo_ids": 115,
        "gap_one_species": 15, "gap_two_species": 10,
        "gap_one_photo_ids": 39, "gap_two_photo_ids": 76,
        "n_images_inspected": 0, "n_flower_images_confirmed": 0,
        "n_colour_labels_measured": 0,
        "original_10km_gate": "HOLD_UNCHANGED",
        "max_metadata_ceiling_validation_species": 35,
        "n_candidate_species_required_to_pass_for_gate": 20,
        "n_candidate_species_allowed_to_fail": 5,
        "notes": [
            "Photo metadata only: review entire 115 ID ledger without morphology-based replacement",
            "A gap-two species requires every missing calendar-year cell to pass image quality checks",
            "No confirmed polymorphism, pigment, genotype, fitness or year-specific climate claim",
            "Prospective 2000+730 taxa and photo identities remain unopened",
        ],
        "species": species, "photo_candidates": photos,
    }


def main():
    ap = argparse.ArgumentParser()
    for x in ("gap1", "gap2", "status", "queue"):
        ap.add_argument("--" + x, required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    a = ap.parse_args()
    d = assemble(*(checked(getattr(a, name), name) for name in ("gap1", "gap2", "status", "queue")))
    a.outdir.mkdir(parents=True, exist_ok=True)
    (a.outdir / "result.json").write_text(json.dumps(d, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    with (a.outdir / "candidate_photo_quality_review.csv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(d["photo_candidates"])
    print(json.dumps({key: d[key] for key in ("schema", "n_candidate_species", "n_photo_ids",
                                              "n_images_inspected", "original_10km_gate")}, sort_keys=True))


if __name__ == "__main__":
    main()
