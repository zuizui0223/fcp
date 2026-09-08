"""Independent stdlib-only verification of a completed Monarda result bundle.

Reads retained counts/statuses only; never imports a model or opens an image.
Mathematical consistency is not annotation truth or model accuracy validation.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from decimal import Decimal, InvalidOperation
import hashlib
import json
import math
from pathlib import Path
import statistics


METRICS = ("prediction_precision", "reference_recall", "iou", "dice")
COUNTS = ("reference_pixels", "predicted_pixels", "intersection_pixels", "union_pixels")


def check(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def integer(value):
    # Pandas writes integral count columns containing the unknown row as N.0.
    # Parse exactly with Decimal; never truncate fractional or nonfinite counts.
    try:
        number = Decimal(str(value))
    except InvalidOperation as exc:
        raise ValueError("invalid retained count") from exc
    check(number.is_finite() and number >= 0 and number == number.to_integral_value(),
          "non-integer or negative retained count")
    return int(number)


def near(left, right, name):
    value = float(left)
    check(math.isfinite(value) and math.isclose(value, right, rel_tol=0, abs_tol=1e-12),
          f"metric mismatch: {name}")


def positive_summary(rows):
    check(len(rows) == 109, "positive-reference denominator is not 109")
    vectors = {key: [] for key in METRICS}
    ref_total = pred_total = inter_total = nonempty = 0
    for row in rows:
        ref, pred, inter, union = (integer(row[key]) for key in COUNTS)
        check(ref > 0 and inter <= min(ref, pred) and union == ref + pred - inter,
              "inconsistent retained pixel counts")
        expected = (inter / pred if pred else 0.0, inter / ref,
                    inter / union, 2 * inter / (ref + pred))
        for key, number in zip(METRICS, expected):
            near(row[key], number, key)
            vectors[key].append(number)
        if row["model_status"].startswith("model_runtime_failure:"):
            check(pred == 0 and inter == 0, "runtime failure did not retain frozen zero-mask penalty")
        else:
            check(row["model_status"] == "model_success", "unknown runtime terminal status")
        ref_total += ref
        pred_total += pred
        inter_total += inter
        nonempty += pred > 0
    result = {
        "positive_reference_images": 109,
        "reference_pixels_total": ref_total, "predicted_pixels_total": pred_total,
        "intersection_pixels_total": inter_total,
        "pooled_prediction_precision": inter_total / pred_total if pred_total else 0.0,
        "pooled_reference_recall": inter_total / ref_total,
        "median_positive_image_prediction_precision": statistics.median(vectors["prediction_precision"]),
        "median_positive_image_reference_recall": statistics.median(vectors["reference_recall"]),
        "median_positive_image_iou": statistics.median(vectors["iou"]),
        "median_positive_image_dice": statistics.median(vectors["dice"]),
        "fraction_positive_images_iou_ge_0_25": sum(x >= .25 for x in vectors["iou"]) / 109,
        "fraction_positive_images_with_nonempty_prediction": nonempty / 109,
    }
    gates = {
        "pooled_prediction_precision_ge_0_70": result["pooled_prediction_precision"] >= .70,
        "pooled_reference_recall_ge_0_35": result["pooled_reference_recall"] >= .35,
        "median_positive_image_prediction_precision_ge_0_70": result["median_positive_image_prediction_precision"] >= .70,
    }
    return {**result, "gate_components": gates, "limited_gate_pass": all(gates.values())}


def csv_rows(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def row_key(row):
    return row["split"], int(row["coco_image_id"])


def verify_bundle(directory: Path, intake_path: Path) -> dict:
    summary_path = directory / "monarda_region_agreement_result_v1.json"
    rows_path = directory / "monarda_region_agreement_rows_v1.csv"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    check(summary["status"] == "complete_monarda_limited_region_agreement_diagnostic",
          "incomplete execution is not eligible for complete-result verification")
    rows = csv_rows(rows_path)
    alignment = csv_rows(directory / "alignment_rows_v1.csv")
    checkpoint = json.loads((directory / "all_110_checkpoint.json").read_text(encoding="utf-8"))
    terminal = json.loads((directory / "execution_terminal.json").read_text(encoding="utf-8"))
    start = json.loads((directory / "execution_start.json").read_text(encoding="utf-8"))
    intake = json.loads(intake_path.read_text(encoding="utf-8"))
    expected_keys = {row_key(r) for r in intake["image_inventory"]}
    expected_rows = {row_key(r): r for r in intake["image_inventory"]}
    for frame in (rows, alignment, checkpoint):
        check(len(frame) == 110 and {row_key(r) for r in frame} == expected_keys,
              "missing, duplicate or replaced image in fixed 110-row census")
    check(Counter(r["split"] for r in rows) == {"train": 77, "valid": 22, "test": 11}, "split census mismatch")
    for row in alignment:
        check(row["alignment_status"] == "alignment_pass", "incomplete reference alignment")
        check(integer(row["decoded_width"]) == integer(row["width_declared"]) and
              integer(row["decoded_height"]) == integer(row["height_declared"]), "dimension mismatch")
    check(all(r["scoring_status"] == "terminal" for r in checkpoint), "partial ledger")
    allowed = {"split", "coco_image_id", "member_path", "annotation_count", "reference_role",
               "model_status", "retained_instances", "automated_colour_state_status_descriptive_only", *COUNTS, *METRICS}
    check(all(set(r) == allowed for r in rows), "unexpected output fields / potential colour leakage")
    known = [r for r in rows if r["reference_role"] == "positive_generic_flower_region"]
    unknown = [r for r in rows if r["reference_role"] == "reference_unknown_empty"]
    check(len(known) == 109 and len(unknown) == 1, "reference-role denominator mismatch")
    check(row_key(unknown[0]) == ("train", 69) and integer(unknown[0]["annotation_count"]) == 0,
          "unknown reference identity changed")
    check(all(unknown[0][k] == "" for k in (*METRICS, "intersection_pixels", "union_pixels")),
          "unknown was scored as a verified negative")
    aligned = {row_key(r): r for r in alignment}
    retained = {row_key(r): r for r in checkpoint}
    for row in rows:
        key = row_key(row)
        source = expected_rows[key]
        check(row["member_path"] == source["path"] and integer(row["annotation_count"]) == source["annotation_count"],
              "image/annotation census identity changed")
        area = source["width_declared"] * source["height_declared"]
        check(integer(row["reference_pixels"]) <= area and integer(row["predicted_pixels"]) <= area,
              "mask pixel count exceeds declared image area")
        check(integer(row["reference_pixels"]) == integer(aligned[key]["reference_pixels"]), "reference count changed after alignment")
        for field in (*COUNTS, *METRICS, "model_status"):
            saved = retained[key][field]
            if saved is None:
                check(row[field] == "", "unknown checkpoint/CSV mismatch")
            elif field == "model_status":
                check(row[field] == saved, "runtime status mismatch")
            else:
                near(row[field], saved, "checkpoint " + field)
    recomputed = positive_summary(known)
    check(set(summary["aggregate"]) == set(recomputed), "aggregate field drift")
    for field, value in recomputed.items():
        if isinstance(value, (dict, bool)):
            check(summary["aggregate"][field] == value, f"aggregate mismatch: {field}")
        else:
            near(summary["aggregate"][field], value, field)
    check(summary["limited_gate_pass"] == recomputed["limited_gate_pass"], "overall gate mismatch")
    for split in ("train", "valid", "test"):
        group = [r for r in known if r["split"] == split]
        reported = summary["provider_split_descriptive_only"][split]
        check(reported["positive_reference_images"] == len(group), "split positive denominator mismatch")
        for short, metric in (("iou", "iou"), ("precision", "prediction_precision"), ("recall", "reference_recall")):
            near(reported["median_" + short], statistics.median(float(r[metric]) for r in group), "split " + metric)
    reported_unknown = summary["unknown_zero_annotation_image"]
    check(reported_unknown["counted_as_verified_negative"] is False and row_key(reported_unknown) == ("train", 69),
          "unknown summary role mismatch")
    check(reported_unknown["predicted_pixels_descriptive_only"] == integer(unknown[0]["predicted_pixels"]) and
          reported_unknown["model_status"] == unknown[0]["model_status"], "unknown descriptive summary mismatch")
    failures = sum(r["model_status"].startswith("model_runtime_failure:") for r in rows)
    check(summary["model_runtime_failure_images"] == failures, "failure count mismatch")
    check(summary["continuous_colour_values_saved_or_analyzed"] is False and summary["ecological_results_changed"] is False,
          "scope/claim mismatch")
    check(summary["lineage"]["archive_sha256"] == intake["archive_sha256"] == start["archive_sha256"], "archive lineage mismatch")
    check(summary["lineage"]["rows_sha256"] == sha(rows_path), "row file hash mismatch")
    check(summary["lineage"]["intake_receipt_sha256"] == sha(intake_path), "intake hash mismatch")
    check(summary["all_110_images_accounted_for"] is True and summary["reference_alignment_complete"] is True,
          "complete-result flags mismatch")
    check(terminal["status"] == "terminal_result_retained" and terminal["automatic_restart_permitted"] is False,
          "terminal execution status mismatch")
    for field, path in (("result_sha256", summary_path), ("checkpoint_sha256", directory / "all_110_checkpoint.json"),
                        ("events_sha256", directory / "events.jsonl")):
        check(terminal[field] == sha(path), f"terminal hash mismatch: {field}")
    events = [json.loads(s) for s in (directory / "events.jsonl").read_text().splitlines()]
    for phase in ("alignment_started", "alignment_terminal", "score_started", "score_terminal"):
        event_keys = [tuple(e["image_key"]) for e in events if e["phase"] == phase]
        check(len(event_keys) == 110 and set(event_keys) == expected_keys, f"incomplete/duplicate event phase: {phase}")
    positions = {(e["phase"], tuple(e["image_key"])): i for i, e in enumerate(events) if e["image_key"] is not None}
    for key in expected_keys:
        check(positions[("alignment_started", key)] < positions[("alignment_terminal", key)] <
              positions[("score_started", key)] < positions[("score_terminal", key)], "per-image event ordering mismatch")
    check(max(positions[("alignment_terminal", k)] for k in expected_keys) <
          min(positions[("score_started", k)] for k in expected_keys), "scoring preceded complete alignment")
    return {"status": "verified_complete_110_row_region_result", "rows": 110,
            "positive_reference_images": 109, "unknown_reference_images": 1,
            "model_runtime_failure_images": failures, "recomputed": recomputed,
            "result_sha256": sha(summary_path), "rows_sha256": sha(rows_path),
            "independent_method": "standard-library integer counts and statistics.median; no runner aggregation import",
            "claim_limit": "Mathematical/provenance consistency, not independent annotation truth or ecological evidence."}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--intake", type=Path, default=Path("docs/supporting/rgfca_monarda_archive_intake_v1.json"))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = verify_bundle(args.directory, args.intake)
    if args.output:
        with args.output.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
