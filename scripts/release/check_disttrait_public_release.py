#!/usr/bin/env python3
"""Hard gate for a public standalone disttrait release."""
from __future__ import annotations

import argparse
import json
import re
import tomllib
from pathlib import Path


PLACEHOLDERS = ("OWNER", "CHOOSE-SPDX-ID", "FIRST", "LAST", "TODO", "TBD")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package-dir", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    args = parser.parse_args()

    package = args.package_dir.resolve()
    metadata_path = args.metadata.resolve()
    data = tomllib.loads(metadata_path.read_text(encoding="utf-8"))
    release = data.get("release", {})

    errors: list[str] = []
    if release.get("release_ready") is not True:
        errors.append("release.release_ready is not true")

    required_files = (
        "LICENSE",
        "CITATION.cff",
        "RELEASE_METADATA.json",
        "README.md",
        "CHANGELOG.md",
        "RELEASING.md",
        "pyproject.toml",
    )
    for name in required_files:
        p = package / name
        if not p.exists() or not p.read_text(encoding="utf-8").strip():
            errors.append(f"missing or empty required file: {name}")

    raw_metadata = metadata_path.read_text(encoding="utf-8")
    for token in PLACEHOLDERS:
        if token in raw_metadata:
            errors.append(f"placeholder remains in metadata: {token}")

    pyproject = (package / "pyproject.toml").read_text(encoding="utf-8")
    if "https://github.com/zuizui0223/fcp" in pyproject:
        errors.append("pyproject still points to the FCP monorepo")
    for field in ("authors =", "maintainers =", "license ="):
        if field not in pyproject:
            errors.append(f"pyproject missing {field.strip()}")

    citation = (package / "CITATION.cff").read_text(encoding="utf-8") if (package / "CITATION.cff").exists() else ""
    for token in PLACEHOLDERS:
        if token in citation:
            errors.append(f"placeholder remains in CITATION.cff: {token}")

    resolved_path = package / "RELEASE_METADATA.json"
    if resolved_path.exists():
        resolved = json.loads(resolved_path.read_text(encoding="utf-8"))
        if resolved.get("release_ready") is not True:
            errors.append("resolved release metadata is not release_ready")

    version_match = re.search(r'^version\s*=\s*"([^"]+)"', pyproject, flags=re.MULTILINE)
    expected_version = str(release.get("version") or "")
    if not version_match or version_match.group(1) != expected_version:
        errors.append("metadata version does not match pyproject version")

    if errors:
        raise RuntimeError("public release gate failed:\n- " + "\n- ".join(errors))

    print(json.dumps({
        "status": "PASS",
        "version": expected_version,
        "repository_url": release.get("repository_url"),
        "license_spdx": release.get("license_spdx"),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
