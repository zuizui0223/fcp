"""Metadata-only mapping of author-indexed filenames; never request photographs."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import subprocess

AUTHOR_SHA = "94c045bc9b03049830761c42c4265b04e790aafea6e10a13223b5541ed2d18a5"
ARCHIVE_SHA = "7c9d213514e319b11cc15623747277f88eb37c453722f665b8570b40b9f0dd0d"
OCCURRENCE_SHA = "39d2d1a4b75f044ba1435637c909458a06b24fd778fca1226f6ff9dd244765e8"
INTAKE_SHA = "1f5acb5af939b83c76b75d6513a0f404019ea88d4a4924d6d1322c8b8cde70f8"
FRAME_HASHES = {
    "29584f3ad7ae0cd99a1d8f43459252af38f615da:data/derived/global_monte_carlo_measured_photos_v1.csv":
        "ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4",
    "9f5abe7b45fcdc20ba83adf75d1a8d4a640f622c:data/derived/rgfca_reserve_replication_measured_photos_v1.csv":
        "0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6",
}


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for part in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(part)
    return h.hexdigest()


def identifier(value):
    value = value.strip()
    if not re.fullmatch(r"[0-9]+(?:\.0+)?", value):
        raise ValueError("Non-integral source identifier")
    return value.split(".")[0]


def frame_metadata(git_ref):
    photo_ids, observation_ids = set(), set()
    rows = 0
    data = subprocess.check_output(["git", "show", git_ref])
    if hashlib.sha256(data).hexdigest() != FRAME_HASHES[git_ref]:
        raise ValueError("Immutable comparison frame hash mismatch")
    with io.StringIO(data.decode("utf-8-sig"), newline="") as stream:
        for row in csv.DictReader(stream):
            photo_ids.add(identifier(row["photo_id"]))
            observation_ids.add(identifier(row["observation_id"]))
            rows += 1
    if rows != 50_000:
        raise ValueError("Comparison must include the complete 50,000-photo frame")
    return {"sha256": hashlib.sha256(data).hexdigest(), "git_ref": git_ref,
            "rows": rows, "photo_ids": photo_ids, "observation_ids": observation_ids}


def audit(intake_path, author_path, occurrence_path):
    if sha(intake_path) != INTAKE_SHA:
        raise ValueError("Pinned complete intake receipt mismatch")
    intake = json.loads(intake_path.read_text(encoding="utf-8"))
    if intake["archive_sha256"] != ARCHIVE_SHA or sha(author_path) != AUTHOR_SHA:
        raise ValueError("Pinned intake/author source mismatch")
    if len(intake["image_inventory"]) != 110:
        raise ValueError("The complete 110-image census is required")
    with author_path.open(encoding="utf-8", newline="") as stream:
        author_rows = list(csv.DictReader(stream, delimiter="\t"))
    if len(author_rows) != 41069:
        raise ValueError("Author multimedia census changed")
    if sha(occurrence_path) != OCCURRENCE_SHA:
        raise ValueError("Pinned author occurrence source mismatch")
    observation_index = {}
    with occurrence_path.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            gbif_id = identifier(row["gbifID"])
            if gbif_id in observation_index:
                raise ValueError("Duplicate GBIF join key")
            observation_index[gbif_id] = row["occurrenceID"]
    # Read the immutable Git objects, not a newline-translated working copy.
    frames = {
        "discovery": frame_metadata("29584f3ad7ae0cd99a1d8f43459252af38f615da:data/derived/global_monte_carlo_measured_photos_v1.csv"),
        "reserve": frame_metadata("9f5abe7b45fcdc20ba83adf75d1a8d4a640f622c:data/derived/rgfca_reserve_replication_measured_photos_v1.csv")
    }
    mapped, licenses = [], Counter()
    photos, observations = defaultdict(list), defaultdict(list)
    for image in intake["image_inventory"]:
        name = image["original_name"]
        if not re.fullmatch(r"[0-9]+\.jpg", name):
            raise ValueError("Original filename does not match author indexing convention")
        source_index = int(name[:-4])
        if source_index >= len(author_rows):
            raise ValueError("Source index outside the author multimedia census")
        source = author_rows[source_index]
        photo = re.search(r"/photos/([0-9]+)/", source["identifier"])
        gbif_id = identifier(source["gbifID"])
        observation_url = observation_index[gbif_id]
        observation = re.fullmatch(r"https?://(?:www\.)?inaturalist\.org/observations/([0-9]+)/?", observation_url)
        if photo is None or observation is None:
            raise ValueError("Original photo/observation identifier could not be parsed")
        photo_id, observation_id = photo.group(1), observation.group(1)
        row = {"path": image["path"], "split": image["split"], "source_row_index_zero_based": source_index,
               "photo_id": photo_id, "observation_id": observation_id, "gbif_id": gbif_id,
               "original_photo_license_from_2025_archive": source["license"],
               "source_image_url": source["identifier"], "source_photo_reference": source["references"], "observation_url": observation_url,
               "overlap": {key: {"photo_id": photo_id in frame["photo_ids"], "observation_id": observation_id in frame["observation_ids"]} for key, frame in frames.items()}}
        mapped.append(row)
        licenses[source["license"]] += 1
        photos[photo_id].append(image["path"])
        observations[observation_id].append(image["path"])
    return {
        "schema_version": "fcp-monarda-author-index-mapping-v1", "archive_sha256": ARCHIVE_SHA,
        "intake_sha256": sha(intake_path), "author_multimedia_sha256": AUTHOR_SHA,
        "author_occurrence_sha256": OCCURRENCE_SHA,
        "author_repository_commit": "49812a431abce78c871f23aff6bc7978573dd502",
        "author_multimedia_git_blob": "84b64ccb89019174a738bdf380528aa79e804fe7",
        "author_occurrence_git_blob": "192b2874235ea4f2457b519cff0f0144fd190081",
        "filename_rule_source_git_blob": "24bd0653f5a3d8cf74cbf89364c91100a94e22b4",
        "filename_mapping_rule": "Notebook 1 reads multimedia.txt with default pandas row index, downloads table.identifier, and writes each file as table.index[idx].jpg. COCO extra.name preserves this filename.",
        "source_selection_caution": "Notebook 3 at blob e579b1c4817c4bd47f472dc3575e14b233a88773 samples 200 GPT flower-present indices via np.random.choice without an explicit seed or replace=False. This does not reproduce or explain selection of the released 110 annotations.",
        "reference_images_mapped": len(mapped), "unique_photo_ids": len(photos), "unique_observation_ids": len(observations),
        "source_license_counts": dict(sorted(licenses.items())),
        "duplicate_photo_id_groups": {k: v for k, v in photos.items() if len(v) > 1},
        "duplicate_observation_id_groups": {k: v for k, v in observations.items() if len(v) > 1},
        "fcp_comparison_frames": {key: {"rows": frame["rows"], "immutable_input_sha256": frame["sha256"], "git_ref": frame["git_ref"],
                                        "photo_id_matches": sum(r["overlap"][key]["photo_id"] for r in mapped),
                                        "observation_id_matches": sum(r["overlap"][key]["observation_id"] for r in mapped)} for key, frame in frames.items()},
        "image_source_mapping": mapped,
        "limits": ["Mapping is documentary linkage through indexed names, not byte identity to reacquired originals.",
                   "The 2025 source licence field is historical, not a current rights verification or permission to relicense photos.",
                   "The raw occurrence source includes locations; only gbifID and occurrenceID were retained for this mapping.",
                   "No coordinate, colour, phenology or observer fields were joined.",
                   "ID disjointness from two FCP frames does not prove absence from JRC, foundation-model training, legacy frames, nearby observations or other participants' datasets.",
                   "No image was decoded or measured; no reference was admitted for model evaluation."],
        "decision": "author_filename_to_source_id_mapping_completed_with_bounded_fcp_id_overlap_check"
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("intake", "author", "occurrence", "output"):
        parser.add_argument(arg, type=Path)
    args = parser.parse_args()
    result = audit(args.intake, args.author, args.occurrence)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({k: result[k] for k in ("reference_images_mapped", "unique_photo_ids", "unique_observation_ids", "source_license_counts", "duplicate_observation_id_groups", "fcp_comparison_frames", "decision")}, indent=2))
