#!/usr/bin/env python3
"""Validate the immutable iNaturalist snapshot identity receipt."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    path = ROOT / "docs/supporting/jbi_atlas_inaturalist_snapshot_receipt_v1.json"
    receipt = json.loads(path.read_text(encoding="utf-8"))
    amendment_path = ROOT / receipt.get("parent_amendment", "")
    amendment = json.loads(amendment_path.read_text(encoding="utf-8"))
    source = receipt.get("source", {})
    integrity = receipt.get("archive_integrity", {})
    if (
        receipt.get("protocol") != "jbi-atlas-inaturalist-snapshot-receipt-v1"
        or receipt.get("status") != "pass_exact_archive_identity_and_tar_integrity"
        or receipt.get("candidate_image_pixels_opened") is not False
        or receipt.get("continuous_colour_used") is not False
        or receipt.get("moving_latest_used") is not False
        or source.get("content_length_bytes") != 35093052336
        or source.get("sha256")
        != "c98202c07796b275fe41fc1518fc394ac09caf2dede370a4ee64ce6d68b0c50d"
        or not str(source.get("official_url", "")).endswith("20260827.tar.gz")
        or integrity.get("full_tar_gzip_listing_exit_code") != 0
        or set(integrity.get("required_reconciliation_tables", ()))
        != {"observations.csv", "observers.csv", "photos.csv", "taxa.csv"}
    ):
        raise RuntimeError("iNaturalist dated snapshot receipt changed")
    members = {Path(value).name for value in integrity.get("members", ())}
    if not set(integrity["required_reconciliation_tables"]).issubset(members):
        raise RuntimeError("dated snapshot lacks a required reconciliation table")
    for field in (
        "snapshot_date",
        "official_url",
        "object_key",
        "content_length_bytes",
        "etag",
        "last_modified_utc",
    ):
        if source.get(field) != amendment.get("snapshot", {}).get(field):
            raise RuntimeError(f"snapshot receipt and amendment disagree: {field}")
    print(json.dumps({"status": "pass", "sha256": source["sha256"]}))


if __name__ == "__main__":
    main()
