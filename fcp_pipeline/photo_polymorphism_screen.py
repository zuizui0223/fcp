"""Photo-first FCP screening from NEW, frozen ROI-v4 image measurements.

42,111 is an outcome-blind opportunity frame, never the number of FCP species.
Classify natural photographic within-species colour-state variation ONLY from
location-blind original flower ROI-v4/EfficientSAM terminal outcomes, not from
literature lists, historic 149,900 colour rows, or metadata-only counts.

The positive tier requires independent observer-disjoint photo evidence.
Single-state detections and undersampled species are NEVER called monomorphic.
Photographic variation alone is not a genetic polymorphism or same-population
coexistence claim.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

MORPHS = ("white", "yellow_orange", "red_pink", "blue_purple")
PALETTE = ("white", "yellow", "orange", "red", "pink", "magenta", "purple", "blue", "bronze")
FRAME_SHA = "5d871bd1f1190c68fd513bb527e6c93de1606dc35ff8cfff5f6187850dc382cc"
SOURCE_COUNT_SHA = "ef26c16252dd785f2d08bde3b18938c7624d73a54386fa3d3974dceb09ed97a9"
SOURCE_PHOTOS = 93_209_351
SNAPSHOT = "2026-09-27"
PROTOCOL = "PHOTO_FCP_OBSERVER_DISJOINT_SCREEN_V1"
METHOD_STATUS = "classified_four_state_morph"
MIN_CLASSIFIABLE = 40
MIN_SPLIT = 20
MIN_MINOR_FRAC = 0.10
MIN_EVIDENCE_OBSERVATIONS_PER_MORPH = 3
MIN_OBSERVERS_PER_MORPH = 2
MIN_OBSERVATIONS_PER_MORPH_PER_SPLIT = 2
POSITIVE = "photo_supported_two_colour_states_observer_disjoint"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _records(path: Path) -> tuple[list[str], list[dict]]:
    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fields = list(reader.fieldnames or ())
        return fields, list(reader)


def load_frame(frame: Path, counts_file: Path, *, synthetic: bool = False) -> dict[int, dict]:
    if not synthetic and (sha256(frame) != FRAME_SHA or sha256(counts_file) != SOURCE_COUNT_SHA):
        raise ValueError("42,111 taxon source frame or 93M source-photo census changed")
    heads, rows = _records(frame)
    if not {"inat_taxon_id", "species"}.issubset(heads):
        raise ValueError("Original 42,111 taxon ID/name columns missing")
    result = {}
    for row in rows:
        idtext = row["inat_taxon_id"]
        if not idtext.isdecimal() or int(idtext) <= 0 or not row["species"].strip():
            raise ValueError("Invalid source taxon")
        ident = int(idtext)
        if ident in result:
            raise ValueError("Duplicate source taxon ID")
        result[ident] = {"species": row["species"].strip(), "n_source_photo_ids": None}
    if len(result) != (len(rows) if synthetic else 42111):
        raise ValueError("42,111 source taxa are not complete")
    cols, source_rows = _records(counts_file)
    if cols != ["inat_taxon_id", "species", "n_unique_source_photo_ids", "snapshot_date"]:
        raise ValueError("Original source photo-count ledger schema drift")
    seen = set()
    for row in source_rows:
        raw, n = row["inat_taxon_id"], row["n_unique_source_photo_ids"]
        if not raw.isdecimal() or not n.isdecimal() or row["snapshot_date"] != SNAPSHOT:
            raise ValueError("Invalid original provider photo count")
        ident = int(raw)
        if ident in seen or ident not in result or row["species"] != result[ident]["species"]:
            raise ValueError("Original photo count taxon ID/name mismatch")
        seen.add(ident)
        result[ident]["n_source_photo_ids"] = int(n)
    if seen != set(result):
        raise ValueError("Original photo opportunity count does not cover all taxa")
    if not synthetic and sum(v["n_source_photo_ids"] for v in result.values()) != SOURCE_PHOTOS:
        raise ValueError("Source photo opportunity total changed")
    return result


def validate_measurement_receipt(receipt: Path, measured: Path, *,
                                 source_taxa: int, synthetic: bool = False) -> dict:
    doc = json.loads(receipt.read_text(encoding="utf-8"))
    schema = doc.get("schema")
    # Screening is separated from a complete all-photo census: it does NOT
    # change the full source universe or imply that unscreened taxa are negative.
    if schema == "fcp_fresh_complete_remeasurement_v1":
        ok = (doc.get("status") == "COMPLETE_FRESH_MEASURED_PHOTO_CENSUS"
              and doc.get("new_photo_table_sha256") == sha256(measured)
              and doc.get("n_source_taxa") == source_taxa
              and doc.get("historical_measurements_reused") is False
              and doc.get("n_terminal_measurements") == doc.get("n_frozen_licensed_photo_ids"))
        kind = "complete_original_photo_census"
    elif schema == "fcp_sealed_photo_screen_terminal_v1":
        ok = (doc.get("status") == "SEALED_PHOTO_SCREEN_UNBIASED_NOT_FINAL"
              and doc.get("source_date") == SNAPSHOT
              and doc.get("source_frame_sha256") == FRAME_SHA
              and doc.get("source_species") == source_taxa
              and doc.get("protocol") == PROTOCOL
              and doc.get("new_photo_table_sha256") == sha256(measured)
              and doc.get("original_ROI_v4_method_attested") is True
              and doc.get("photo_ID_selection_pre_colour") is True
              and doc.get("source_photo_IDs_sealed") is True
              and doc.get("historical_measured_flower_colours_reused") is False
              and doc.get("n_terminal_measurements", -1) >= 0)
        kind = "predeclared_bounded_photo_screen"
    else:
        raise ValueError("No genuine sealed original-photo FCP measurement receipt")
    if not ok:
        raise ValueError("Incomplete or changed ORIGINAL image measurements; refuse photo-polymorphism claim")
    if synthetic:
        return {"kind": kind, "total": doc["n_terminal_measurements"],
                "receipt_sha256": sha256(receipt)}
    if (source_taxa != 42111 or
        (kind == "complete_original_photo_census" and
         doc.get("n_frozen_licensed_photo_ids") != SOURCE_PHOTOS)):
        raise ValueError("Original 42111-species source photo census does not agree with measured table")
    return {"kind": kind, "total": doc["n_terminal_measurements"],
            "receipt_sha256": sha256(receipt)}


def observer_split(taxon: int, observer_id: str) -> int:
    key = f"FCP_PHOTO_SCREEN_20261016_V1|{taxon}|{observer_id}"
    return int(hashlib.sha256(key.encode("utf-8")).hexdigest()[:16], 16) & 1


def _as_classifiable(row: dict) -> bool:
    flag = str(row["global_classifiable"]).lower()
    if flag not in ("true", "false"):
        raise ValueError("Invalid frozen photo classifiability flag")
    state = str(row["morph"])
    status = str(row["measurement_status"])
    if (flag == "true") != (status == METHOD_STATUS and state in MORPHS):
        raise ValueError("Photo morph/classifiable/ROI status mismatch")
    if flag == "false" and state != "mixed_uncertain":
        raise ValueError("Failed flower ROI must not be a biological colour")
    return flag == "true"


def _finalize(taxon: int, name: str, n_source: int,
              photos: int, good: int, counts: Counter,
              obs_to_labels: dict[str, set[str]],
              obs_to_observer: dict[str, str],
              split_photo: list[Counter]) -> dict:
    obs_count = Counter()
    observers = defaultdict(set)
    split_obs = [Counter(), Counter()]
    contradictory = 0
    # Two different coarse labels within one iNaturalist observation can mean
    # multiple pictured flowers, false masks or different development stages.
    # They remain in photo counts but cannot alone support an independent morph.
    for obs, states in obs_to_labels.items():
        if len(states) != 1:
            contradictory += 1
            continue
        state = next(iter(states))
        owner = obs_to_observer[obs]
        obs_count[state] += 1
        observers[state].add(owner)
        split_obs[observer_split(taxon, owner)][state] += 1
    states = sorted((m for m in MORPHS if counts[m] > 0),
                    key=lambda m: (-counts[m], MORPHS.index(m)))
    second = (counts[states[1]] / good) if len(states) > 1 else 0.0
    p = [counts[m] / good for m in MORPHS] if good else [0.0] * 4
    d = 1.0 - sum(x * x for x in p) if good >= MIN_CLASSIFIABLE else None
    fold_photos = [sum(x.values()) for x in split_photo]
    nominated = [
        pair for i, a in enumerate(states) for b in states[i + 1:]
        if (min(counts[a], counts[b]) / good >= MIN_MINOR_FRAC and
            all(obs_count[m] >= MIN_EVIDENCE_OBSERVATIONS_PER_MORPH and
                len(observers[m]) >= MIN_OBSERVERS_PER_MORPH for m in (a, b)))
        for pair in [(a, b)]
    ] if good else []
    replicated = [
        pair for pair in nominated
        if all(fold_photos[i] >= MIN_SPLIT and
               all(split_obs[i][m] >= MIN_OBSERVATIONS_PER_MORPH_PER_SPLIT
                   for m in pair) for i in (0, 1))
    ]
    if photos == 0:
        status = "not_screened_no_terminal_photo"
    elif good < MIN_CLASSIFIABLE:
        status = "insufficient_classifiable_not_monomorphic"
    elif replicated:
        status = POSITIVE
    elif nominated:
        status = "photo_variation_candidate_holdout_unresolved"
    elif len(states) > 1:
        status = "photo_colour_signal_weak_or_unreplicated"
    else:
        status = "single_photo_state_detected_not_monomorphic"
    return {
        "inat_taxon_id": taxon, "species": name,
        "source_photo_opportunities": n_source,
        "n_terminal_photo_records": photos, "n_classifiable": good,
        "n_classifiable_observations_consistent": sum(obs_count.values()),
        "n_observations_conflicting_photo_labels": contradictory,
        "n_observer_split_0_classifiable_photos": fold_photos[0],
        "n_observer_split_1_classifiable_photos": fold_photos[1],
        "n_detected_coarse_states": len(states),
        "coarse_state_counts": json.dumps({m: counts[m] for m in MORPHS}, sort_keys=True),
        "coarse_state_observation_counts": json.dumps(dict(obs_count), sort_keys=True),
        "coarse_state_distinct_observers": json.dumps(
            {m: len(observers[m]) for m in MORPHS}, sort_keys=True),
        "minor_coarse_photo_fraction": round(second, 12),
        "D_photo_equal": round(d, 12) if d is not None else "",
        "photographic_colour_status": status,
        "heldout_replicated_pair": json.dumps(replicated[0]) if replicated else "",
        "genetic_polymorphism_proven": False,
        "local_population_coexistence_proven": False,
    }


def screen(measured: Path, frame: Path, counts: Path, receipt: Path,
           outdir: Path, *, _synthetic: bool = False) -> dict:
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError("Cannot overwrite existing photo-based polymorphism cohort")
    source = load_frame(frame, counts, synthetic=_synthetic)
    evidence = validate_measurement_receipt(receipt, measured, source_taxa=len(source),
                                            synthetic=_synthetic)
    result, positive, unresolved = [], [], []
    current = None
    photo_count = good_count = 0
    morph_counts = Counter()
    obs_to_labels = {}
    obs_to_observer = {}
    split_photo = [Counter(), Counter()]
    photo_ids = set()
    processed = 0

    def flush() -> None:
        if current is None:
            return
        row = _finalize(current, source[current]["species"],
                        source[current]["n_source_photo_ids"],
                        photo_count, good_count, morph_counts, obs_to_labels,
                        obs_to_observer, split_photo)
        result.append(row)
        (positive if row["photographic_colour_status"] == POSITIVE
         else unresolved).append(row)

    with measured.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            "inat_taxon_id", "species", "photo_id", "observation_id", "observer_id",
            "measurement_status", "morph", "global_classifiable",
            *(f"flower_fraction_{p}" for p in PALETTE),
        }
        if not required.issubset(reader.fieldnames or ()):
            raise ValueError("Measurement must contain ORIGINAL flower ROI-v4 nine-palette columns")
        for row in reader:
            raw = str(row["inat_taxon_id"])
            if not raw.isdecimal() or int(raw) not in source:
                raise ValueError("Unknown source photo species")
            taxon = int(raw)
            if row["species"] != source[taxon]["species"]:
                raise ValueError("Fresh photo taxon name drift")
            if current is not None and taxon < current:
                raise ValueError("Original screen must be taxon sorted; no opportunistic reselection")
            if current != taxon:
                flush()
                current = taxon
                photo_count = good_count = 0
                morph_counts = Counter()
                obs_to_labels = {}
                obs_to_observer = {}
                split_photo = [Counter(), Counter()]
                photo_ids = set()
            pid = str(row["photo_id"]).strip()
            if not pid or pid in photo_ids:
                raise ValueError("Duplicate/blank original photo ID in species")
            photo_ids.add(pid)
            photo_count += 1
            processed += 1
            if not _as_classifiable(row):
                continue
            obs = str(row["observation_id"]).strip()
            observer = str(row["observer_id"]).strip()
            if not obs or not observer:
                raise ValueError("Missing independent observation or photographer ID")
            previous = obs_to_observer.setdefault(obs, observer)
            if previous != observer:
                raise ValueError("One original observation attributed to multiple observers")
            morph = row["morph"]
            vals = [float(row[f"flower_fraction_{p}"]) for p in PALETTE]
            if any(not math.isfinite(v) or v < 0 or v > 1 for v in vals):
                raise ValueError("Original colour fractions invalid")
            good_count += 1
            morph_counts[morph] += 1
            obs_to_labels.setdefault(obs, set()).add(morph)
            split_photo[observer_split(taxon, observer)][morph] += 1
    flush()
    if processed != evidence["total"]:
        raise ValueError("New original photo terminal count does not match sealed receipt")
    visited = {x["inat_taxon_id"] for x in result}
    for taxon, data in source.items():
        if taxon not in visited:
            row = _finalize(taxon, data["species"], data["n_source_photo_ids"],
                            0, 0, Counter(), {}, {}, [Counter(), Counter()])
            result.append(row)
            unresolved.append(row)
    result.sort(key=lambda x: x["inat_taxon_id"])
    positive.sort(key=lambda x: x["inat_taxon_id"])
    unresolved.sort(key=lambda x: x["inat_taxon_id"])
    from csv import DictWriter
    outdir.mkdir(parents=True, exist_ok=True)
    fields = list(result[0])
    for name, rows in (
        ("all_42111_photo_polymorphism_screen.csv", result),
        ("photo_supported_polymorphic_species.csv", positive),
        ("not_proven_polymorphic_do_not_call_monomorphic.csv", unresolved),
    ):
        with (outdir / name).open("w", newline="", encoding="utf-8") as dest:
            writer = DictWriter(dest, fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    classes = dict(Counter(x["photographic_colour_status"] for x in result))
    summary = {
        "schema": "fcp_photo_derived_polymorphism_screen_v1",
        "status": "IMAGE_DERIVED_PHENOTYPIC_VARIATION_ONLY_NOT_GENETIC_FCP",
        "method": "original_ROI_v4_EfficientSAM_four_state_location_blind",
        "source_measurement_kind": evidence["kind"],
        "sealed_original_measurement_receipt_sha256": evidence["receipt_sha256"],
        "original_colour_table_sha256": sha256(measured),
        "original_source_frame_sha256": sha256(frame),
        "source_photo_counts_sha256": sha256(counts),
        "source_species": len(source),
        "original_photo_terminals_screened": processed,
        "photo_supported_polymorphic_species": len(positive),
        "photo_states_not_literature_sources_decide_inclusion": True,
        "heldout_observer_disjoint_validation_required": True,
        "classifiable_min": MIN_CLASSIFIABLE,
        "min_photo_minor_fraction": MIN_MINOR_FRAC,
        "min_independent_observations_per_morph": MIN_EVIDENCE_OBSERVATIONS_PER_MORPH,
        "min_distinct_observers_per_morph": MIN_OBSERVERS_PER_MORPH,
        "per_split_classifiable_min": MIN_SPLIT,
        "genetic_FCP_established": False,
        "photographic_single_state_equals_biological_monomorphism": False,
        "previous_34_literature_C_S_labels_used": False,
        "historical_149900_measured_colour_rows_reused": False,
        "complete_42111_photo_census": evidence["kind"] == "complete_original_photo_census",
        "source_status_counts": classes,
        "warning": ("A photo-positive requires two coarse morphs with independent "
                    "observer-held-out replication. Nonpositives are never evidence "
                    "of biological monomorphism. Local coexistence requires separate "
                    "geographically and population-aware verification."),
    }
    (outdir / "photo_screen_receipt.json").write_text(
        json.dumps(summary, ensure_ascii=False, sort_keys=True, indent=2) + "\n",
        encoding="utf-8")
    return summary


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--measured", required=True, type=Path)
    p.add_argument("--source-frame", required=True, type=Path)
    p.add_argument("--source-counts", required=True, type=Path)
    p.add_argument("--measurement-receipt", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    args = p.parse_args(argv)
    print(json.dumps(screen(args.measured, args.source_frame, args.source_counts,
                            args.measurement_receipt, args.outdir),
                     ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
