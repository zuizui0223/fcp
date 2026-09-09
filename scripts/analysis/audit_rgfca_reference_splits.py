"""Audit saved public filename lists; never open images, models or the network.

This is a retrospective source-integrity audit, not a benchmark authorization.
Exact filename disjointness is necessary but insufficient for image, event,
observer or foundation-training independence. Original list rows are retained.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re


MAX_LIST_BYTES = 1_000_000


def parse_names(raw: bytes) -> list[str]:
    """Parse one basename per line, without repairing or deduplicating input."""
    if not raw or len(raw) > MAX_LIST_BYTES:
        raise ValueError("empty or oversized filename list")
    rows = raw.decode("ascii").replace("\r\n", "\n").split("\n")
    if rows[-1] == "":
        rows.pop()
    if any(not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*\.[A-Za-z0-9]+", row)
           for row in rows):
        raise ValueError("expected unpadded ASCII basenames, one per line")
    return rows


def summarize_splits(train: list[str], validation: list[str]) -> dict:
    """Keep exact case-sensitive keys; report casefold collisions separately."""
    groups = {"train": train, "validation": validation}
    profile = {}
    for name, rows in groups.items():
        if not rows:
            raise ValueError("both splits must contain rows")
        counts = Counter(rows)
        profile[name] = {
            "row_count": len(rows),
            "unique_filename_count": len(counts),
            "duplicate_extra_row_count": len(rows) - len(counts),
            "duplicate_filenames": {key: count for key, count in sorted(counts.items())
                                    if count > 1},
        }
    train_keys, validation_keys = set(train), set(validation)
    overlap = sorted(train_keys & validation_keys)
    casefold_overlap = sorted({key.casefold() for key in train_keys}
                             & {key.casefold() for key in validation_keys})
    variants: dict[str, set[str]] = {}
    for key in train_keys | validation_keys:
        variants.setdefault(key.casefold(), set()).add(key)
    collisions = {key: sorted(values) for key, values in sorted(variants.items())
                  if len(values) > 1}
    return {
        "profiles": profile,
        "exact_shared_filenames": overlap,
        "exact_shared_filename_count": len(overlap),
        "exact_union_filename_count": len(train_keys | validation_keys),
        "casefold_shared_keys": casefold_overlap,
        "casefold_spelling_collisions": collisions,
        "filename_split_disjoint": not overlap,
        "filename_integrity_gate_passed": (
            not casefold_overlap and not collisions
            and all(p["duplicate_extra_row_count"] == 0 for p in profile.values())
        ),
    }


def audit_bundle(bundle: Path) -> dict:
    """Hash and audit only the two exact text files named by the local manifest."""
    bundle = bundle.resolve()
    manifest_raw = (bundle / "source_manifest.json").read_bytes()
    manifest = json.loads(manifest_raw)
    if manifest.get("schema_version") != "rgfca-reference-split-source-v1":
        raise ValueError("unknown source manifest schema")
    sources = manifest["files"]
    if not isinstance(sources, list) or len(sources) != 2:
        raise ValueError("exactly two source lists required")
    if {source["role"] for source in sources} != {"train", "validation"}:
        raise ValueError("unique train and validation roles required")
    rows, receipts = {}, []
    for source in sources:
        name = source["local_filename"]
        if not isinstance(name, str) or not re.fullmatch(r"[A-Za-z0-9_.-]+\.txt", name):
            raise ValueError("source must be a local text basename")
        path = (bundle / name).resolve()
        if path.parent != bundle:
            raise ValueError("source escaped bundle directory")
        if path.stat().st_size > MAX_LIST_BYTES:
            raise ValueError("oversized source")
        raw = path.read_bytes()
        sha256 = hashlib.sha256(raw).hexdigest()
        md5 = hashlib.md5(raw, usedforsecurity=False).hexdigest()
        if len(raw) != source["bytes"] or sha256 != source["sha256"]:
            raise ValueError("source size or SHA-256 mismatch")
        if md5 != source["provider_md5"]:
            raise ValueError("provider MD5 mismatch")
        rows[source["role"]] = parse_names(raw)
        receipts.append({"role": source["role"], "local_filename": name,
                         "bytes": len(raw), "sha256": sha256,
                         "provider_md5_verified": True})
    return {
        "schema_version": "rgfca-reference-split-audit-v1",
        "audit_kind": "retrospective_public_filename_split_integrity",
        "source_record": manifest["source_record"],
        "source_manifest_sha256": hashlib.sha256(manifest_raw).hexdigest(),
        "source_receipts": receipts,
        **summarize_splits(rows["train"], rows["validation"]),
        "image_content_identity_verified": False,
        "event_observer_training_independence_verified": False,
        "benchmark_execution_authorized": False,
        "image_pixels_decoded": False,
        "model_run": False,
        "coordinates_joined": False,
        "ecological_inference_performed": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path)
    parser.add_argument("--verify", type=Path,
                        help="compare all fields with a saved audit; no writes")
    args = parser.parse_args()
    result = audit_bundle(args.bundle)
    if args.verify is not None:
        if json.loads(args.verify.read_bytes()) != result:
            raise ValueError("saved audit differs from independently recomputed result")
    # Success means the audit ran faithfully, NOT that the split gate passed.
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
