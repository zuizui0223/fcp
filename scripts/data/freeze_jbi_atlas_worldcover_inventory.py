#!/usr/bin/env python3
"""Freeze the complete official ESA WorldCover 2021 v200 map inventory."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET


BUCKET_URL = "https://esa-worldcover.s3.eu-central-1.amazonaws.com/"
PREFIX = "v200/2021/map/"
USER_AGENT = "FCP-image-first-atlas/3.0 (environment inventory freeze)"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = []
    continuation = ""
    while True:
        query = {"list-type": "2", "prefix": PREFIX, "max-keys": "1000"}
        if continuation:
            query["continuation-token"] = continuation
        request = Request(BUCKET_URL + "?" + urlencode(query), headers={"User-Agent": USER_AGENT})
        with urlopen(request, timeout=60) as response:
            root = ET.fromstring(response.read())
        namespace = {"s3": "http://s3.amazonaws.com/doc/2006-03-01/"}
        for item in root.findall("s3:Contents", namespace):
            key = str(item.findtext("s3:Key", namespaces=namespace))
            if not key.endswith("_Map.tif"):
                raise RuntimeError(f"unexpected object in WorldCover map prefix: {key}")
            rows.append(
                {
                    "key": key,
                    "size_bytes": int(item.findtext("s3:Size", namespaces=namespace)),
                    "etag": str(item.findtext("s3:ETag", namespaces=namespace)).strip('"'),
                    "last_modified": str(
                        item.findtext("s3:LastModified", namespaces=namespace)
                    ),
                    "https_url": BUCKET_URL + key,
                }
            )
        truncated = root.findtext("s3:IsTruncated", namespaces=namespace) == "true"
        if not truncated:
            break
        continuation = str(
            root.findtext("s3:NextContinuationToken", namespaces=namespace) or ""
        )
        if not continuation:
            raise RuntimeError("truncated S3 inventory omitted its continuation token")
    rows.sort(key=lambda row: row["key"])
    if len({row["key"] for row in rows}) != len(rows):
        raise RuntimeError("WorldCover inventory contains duplicate keys")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    inventory = args.output_dir / "worldcover_2021_v200_map_inventory.csv"
    with inventory.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        "status": "pass_worldcover_inventory_freeze",
        "source": "ESA WorldCover 2021 v200 public AWS COG prefix",
        "bucket_url": BUCKET_URL,
        "prefix": PREFIX,
        "objects": len(rows),
        "total_size_bytes": sum(row["size_bytes"] for row in rows),
        "inventory_sha256": sha256(inventory),
        "scaleout_colour_opened": False,
        "claim_ceiling": "Source inventory only; no land-cover pixels or flower-colour outcomes were read.",
    }
    manifest_path = args.output_dir / "worldcover_2021_v200_inventory_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
