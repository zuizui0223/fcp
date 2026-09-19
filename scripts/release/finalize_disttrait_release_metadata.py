#!/usr/bin/env python3
"""Finalize public release metadata for a standalone disttrait package.

This script does not choose a licence or authorship. It only validates explicit
metadata and renders pyproject/CITATION metadata from those decisions.
"""
from __future__ import annotations

import argparse
import json
import re
import tomllib
from pathlib import Path


PLACEHOLDERS = (
    "OWNER",
    "CHOOSE-SPDX-ID",
    "FIRST",
    "LAST",
    "TODO",
    "TBD",
)


def _q(value: str) -> str:
    return json.dumps(str(value), ensure_ascii=False)


def _require_text(value: object, field: str) -> str:
    text = str(value or "").strip()
    if not text:
        raise RuntimeError(f"missing required metadata: {field}")
    if any(token in text for token in PLACEHOLDERS):
        raise RuntimeError(f"placeholder remains in {field}: {text!r}")
    return text


def _person_name(person: dict, label: str) -> str:
    given = _require_text(person.get("given_names"), f"{label}.given_names")
    family = _require_text(person.get("family_names"), f"{label}.family_names")
    return f"{given} {family}"


def _load_metadata(path: Path) -> dict:
    data = tomllib.loads(path.read_text(encoding="utf-8"))
    release = data.get("release", {})
    if release.get("release_ready") is not True:
        raise RuntimeError("release.release_ready must be true before finalization")
    authors = data.get("authors", [])
    maintainers = data.get("maintainers", [])
    if not authors:
        raise RuntimeError("at least one [[authors]] entry is required")
    if not maintainers:
        raise RuntimeError("at least one [[maintainers]] entry is required")
    return data


def _render_people_pyproject(people: list[dict], label: str) -> str:
    rendered = []
    for i, person in enumerate(people, start=1):
        name = _person_name(person, f"{label}[{i}]")
        email = str(person.get("email") or "").strip()
        entry = f"{{name = {_q(name)}"
        if email:
            entry += f", email = {_q(email)}"
        entry += "}"
        rendered.append(entry)
    return "[" + ", ".join(rendered) + "]"


def _update_pyproject(package_dir: Path, data: dict) -> None:
    path = package_dir / "pyproject.toml"
    text = path.read_text(encoding="utf-8")
    release = data["release"]
    package_name = _require_text(release.get("package_name"), "release.package_name")
    version = _require_text(release.get("version"), "release.version")
    repository_url = _require_text(release.get("repository_url"), "release.repository_url")

    name_match = re.search(r'^name\s*=\s*"([^"]+)"', text, flags=re.MULTILINE)
    version_match = re.search(r'^version\s*=\s*"([^"]+)"', text, flags=re.MULTILINE)
    if not name_match or name_match.group(1) != package_name:
        raise RuntimeError("package name does not match pyproject")
    if not version_match or version_match.group(1) != version:
        raise RuntimeError("package version does not match pyproject")
    if repository_url.rstrip("/") == "https://github.com/zuizui0223/fcp":
        raise RuntimeError("standalone repository URL must not remain the FCP monorepo URL")

    authors_line = "authors = " + _render_people_pyproject(data["authors"], "authors")
    maintainers_line = "maintainers = " + _render_people_pyproject(
        data["maintainers"], "maintainers"
    )
    license_line = 'license = { file = "LICENSE" }'

    for key in ("authors", "maintainers", "license"):
        text = re.sub(
            rf"^{key}\s*=.*\n",
            "",
            text,
            flags=re.MULTILINE,
        )

    readme_match = re.search(r'^readme\s*=.*$', text, flags=re.MULTILINE)
    if readme_match is None:
        raise RuntimeError("pyproject readme field not found")
    insert = "\n".join((authors_line, maintainers_line, license_line))
    pos = readme_match.end()
    text = text[:pos] + "\n" + insert + text[pos:]

    issue_url = repository_url.rstrip("/") + "/issues"
    if "[project.urls]" not in text:
        text += "\n[project.urls]\n"
    text = re.sub(
        r'^Repository\s*=.*$',
        f"Repository = {_q(repository_url)}",
        text,
        flags=re.MULTILINE,
    )
    text = re.sub(
        r'^Issues\s*=.*$',
        f"Issues = {_q(issue_url)}",
        text,
        flags=re.MULTILINE,
    )
    if not re.search(r'^Repository\s*=', text, flags=re.MULTILINE):
        text += f"Repository = {_q(repository_url)}\n"
    if not re.search(r'^Issues\s*=', text, flags=re.MULTILINE):
        text += f"Issues = {_q(issue_url)}\n"

    path.write_text(text, encoding="utf-8")


def _yaml_scalar(value: str) -> str:
    return json.dumps(str(value), ensure_ascii=False)


def _render_citation(package_dir: Path, data: dict) -> None:
    release = data["release"]
    title = _require_text(release.get("citation_title"), "release.citation_title")
    version = _require_text(release.get("version"), "release.version")
    date = _require_text(release.get("date_released"), "release.date_released")
    repository = _require_text(release.get("repository_url"), "release.repository_url")
    license_spdx = _require_text(release.get("license_spdx"), "release.license_spdx")
    doi = str(release.get("archival_doi") or "").strip()

    lines = [
        'cff-version: "1.2.0"',
        'message: "If you use this software, please cite it using this metadata."',
        f"title: {_yaml_scalar(title)}",
        "type: software",
        f"version: {_yaml_scalar(version)}",
        f"date-released: {_yaml_scalar(date)}",
        f"repository-code: {_yaml_scalar(repository)}",
        f"license: {_yaml_scalar(license_spdx)}",
        "authors:",
    ]
    for i, person in enumerate(data["authors"], start=1):
        given = _require_text(person.get("given_names"), f"authors[{i}].given_names")
        family = _require_text(person.get("family_names"), f"authors[{i}].family_names")
        lines.append(f"  - given-names: {_yaml_scalar(given)}")
        lines.append(f"    family-names: {_yaml_scalar(family)}")
        affiliation = str(person.get("affiliation") or "").strip()
        orcid = str(person.get("orcid") or "").strip()
        if affiliation:
            lines.append(f"    affiliation: {_yaml_scalar(affiliation)}")
        if orcid:
            lines.append(f"    orcid: {_yaml_scalar(orcid)}")
    if doi:
        lines.append(f"doi: {_yaml_scalar(doi)}")
    (package_dir / "CITATION.cff").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--package-dir", type=Path, required=True)
    args = parser.parse_args()

    metadata_path = args.metadata.resolve()
    package_dir = args.package_dir.resolve()
    data = _load_metadata(metadata_path)
    release = data["release"]

    licence = package_dir / _require_text(release.get("license_file"), "release.license_file")
    if not licence.exists() or not licence.read_text(encoding="utf-8").strip():
        raise RuntimeError(f"selected licence file is missing or empty: {licence}")

    _update_pyproject(package_dir, data)
    _render_citation(package_dir, data)

    resolved = {
        "schema": "disttrait_public_release_metadata_v1",
        "package_name": release["package_name"],
        "version": release["version"],
        "repository_url": release["repository_url"],
        "repository_visibility": release["repository_visibility"],
        "license_spdx": release["license_spdx"],
        "release_tag": release["release_tag"],
        "archival_doi": release.get("archival_doi", ""),
        "authors": data["authors"],
        "maintainers": data["maintainers"],
        "release_ready": True,
    }
    (package_dir / "RELEASE_METADATA.json").write_text(
        json.dumps(resolved, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"status": "PASS", "release_ready": True}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
