#!/usr/bin/env python3
"""Submit the prepared broad GBIF occurrence-download request.

Credentials are never written to disk. Supply the username and notification email
as arguments and the password through the GBIF_PASSWORD environment variable.
The script replaces the placeholder notification email, submits the request using
HTTP Basic authentication and saves the returned download key in a user-selected
output file.
"""
from __future__ import annotations

import argparse
import base64
import json
import os
import urllib.error
import urllib.request
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--request",
        type=Path,
        default=Path("docs/supporting/jbi_gbif_doi_bundle/jbi_gbif_broad_download_request.json"),
    )
    parser.add_argument("--username", required=True)
    parser.add_argument("--email", required=True)
    parser.add_argument("--out", type=Path, default=Path("jbi_gbif_download_key.txt"))
    args = parser.parse_args()

    password = os.environ.get("GBIF_PASSWORD")
    if not password:
        raise SystemExit("GBIF_PASSWORD is required and must not be committed")
    request = json.loads(args.request.read_text(encoding="utf-8"))
    request["notificationAddresses"] = [args.email]

    credentials = base64.b64encode(f"{args.username}:{password}".encode()).decode()
    http_request = urllib.request.Request(
        "https://api.gbif.org/v1/occurrence/download/request",
        data=json.dumps(request).encode("utf-8"),
        method="POST",
        headers={
            "Authorization": f"Basic {credentials}",
            "Content-Type": "application/json",
            "Accept": "text/plain, application/json",
            "User-Agent": "fcp-jbi-gbif-download/1.0",
        },
    )
    try:
        with urllib.request.urlopen(http_request, timeout=120) as response:
            key = response.read().decode("utf-8").strip()
    except urllib.error.HTTPError as error:
        detail = error.read().decode("utf-8", errors="replace")
        raise SystemExit(f"GBIF request failed ({error.code}): {detail}") from error
    if not key:
        raise SystemExit("GBIF returned an empty download key")
    args.out.write_text(key + "\n", encoding="utf-8")
    print(f"GBIF download key: {key}")
    print(f"Saved to: {args.out}")
    print("Wait for SUCCEEDED status, then record the GBIF download DOI in the manuscript and derived-dataset metadata.")


if __name__ == "__main__":
    main()
