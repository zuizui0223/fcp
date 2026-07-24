#!/usr/bin/env python3
"""Build a deterministic Journal of Biogeography submission package and manifest.

Preview mode permits explicit author-controlled placeholders and reports them.
Strict mode fails if journal-facing files retain unresolved placeholders or if the
working tree has substantive uncommitted changes. Known generated validation and
release outputs are excluded from that dirty-tree decision. The repository release
itself should still be archived through a permanent service such as Zenodo.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
FIXED_ZIP_TIME = (1980, 1, 1, 0, 0, 0)
GENERATED_REPORTS = {"docs/jbi_submission_validation_report.json"}

EXPLICIT_FILES = [
    "docs/jbi_manuscript_editorial_revision.md",
    "docs/jbi_title_page_template.md",
    "docs/jbi_cover_letter_template.md",
    "docs/jbi_supporting_information_index.md",
    "docs/jbi_literature_search_provenance.md",
    "docs/jbi_taxon_image_candidate.md",
    "docs/jbi_gbif_doi_protocol.md",
    "docs/jbi_editorial_submission_check.md",
    "docs/jbi_submission_completion_checklist.md",
    "docs/jbi_author_confirmation_form.md",
    "docs/jbi_archive_release_protocol.md",
    "docs/jbi_zenodo_metadata_template.json",
    "analysis_phylogenetic_spatial_scale.R",
    "analysis_vphylomaker2_dated_spatial_scale.R",
    "build_gbif_occurrence_sample_paginated.py",
    "prepare_jbi_gbif_doi_bundle.py",
    "submit_jbi_gbif_download.py",
    "validate_jbi_submission_package.py",
    "validate_jbi_gbif_doi_bundle.py",
]

GLOB_PATTERNS = [
    "docs/figures/*",
    "docs/supporting/**/*",
]

STRICT_TEXT_FILES = {
    "docs/jbi_manuscript_editorial_revision.md": [r"\bNot verified\b"],
    "docs/jbi_title_page_template.md": [r"\bNot verified\b", r"REPLACE_", r"\[[^\]]+\]"],
    "docs/jbi_cover_letter_template.md": [r"REPLACE_", r"\[[^\]]+\]"],
    "docs/jbi_zenodo_metadata_template.json": [r"REPLACE_"],
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def git_text(*args: str) -> str:
    result = subprocess.run(
        ["git", *args],
        cwd=ROOT,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return result.stdout.strip()


def collect_files() -> list[Path]:
    files: set[Path] = set()
    for relative in EXPLICIT_FILES:
        files.add(ROOT / relative)
    for pattern in GLOB_PATTERNS:
        for path in ROOT.glob(pattern):
            if path.is_file():
                files.add(path)
    return sorted(files, key=lambda p: p.relative_to(ROOT).as_posix())


def placeholder_findings(files: list[Path]) -> list[dict[str, object]]:
    findings: list[dict[str, object]] = []
    for path in files:
        relative = path.relative_to(ROOT).as_posix()
        patterns = STRICT_TEXT_FILES.get(relative, [])
        if not patterns:
            continue
        text = path.read_text(encoding="utf-8")
        for pattern in patterns:
            matches = sorted(set(re.findall(pattern, text, flags=re.IGNORECASE)))
            if matches:
                findings.append({"path": relative, "pattern": pattern, "matches": matches[:20]})
    return findings


def validate_required_files(files: list[Path]) -> None:
    missing = [p.relative_to(ROOT).as_posix() for p in files if not p.is_file() or p.stat().st_size == 0]
    if missing:
        raise SystemExit(f"Missing or empty release files: {missing}")


def porcelain_path(line: str) -> str:
    value = line[3:].strip() if len(line) >= 4 else line.strip()
    if " -> " in value:
        value = value.rsplit(" -> ", 1)[1]
    return value.strip('"')


def substantive_status(lines: list[str], output_dir_relative: str) -> list[str]:
    output_prefix = output_dir_relative.rstrip("/") + "/"
    retained: list[str] = []
    for line in lines:
        path = porcelain_path(line)
        if path in GENERATED_REPORTS:
            continue
        if path == output_dir_relative.rstrip("/") or path.startswith(output_prefix):
            continue
        retained.append(line)
    return retained


def write_deterministic_zip(files: list[Path], destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(destination, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for path in files:
            relative = path.relative_to(ROOT).as_posix()
            info = zipfile.ZipInfo(relative, FIXED_ZIP_TIME)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            archive.writestr(info, path.read_bytes(), compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="analysis_outputs/jbi_release")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--strict", action="store_true")
    mode.add_argument("--allow-placeholders", action="store_true")
    args = parser.parse_args()

    files = collect_files()
    validate_required_files(files)

    commit_sha = git_text("rev-parse", "HEAD")
    raw_status = [line for line in git_text("status", "--porcelain").splitlines() if line]
    effective_status = substantive_status(raw_status, Path(args.output_dir).as_posix())
    findings = placeholder_findings(files)

    if args.strict and effective_status:
        raise SystemExit(f"Strict release requires no substantive uncommitted changes: {effective_status}")
    if args.strict and findings:
        raise SystemExit(f"Strict release blocked by unresolved placeholders: {findings}")

    output_dir = ROOT / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    archive_path = output_dir / "jbi_submission_files.zip"
    manifest_path = output_dir / "jbi_release_manifest.json"

    write_deterministic_zip(files, archive_path)

    file_records = [
        {
            "path": path.relative_to(ROOT).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in files
    ]
    manifest = {
        "status": "strict-pass" if args.strict else "preview",
        "source_commit": commit_sha,
        "working_tree_clean_for_release": not effective_status,
        "raw_working_tree_status": raw_status,
        "substantive_working_tree_status": effective_status,
        "ignored_generated_paths": sorted(GENERATED_REPORTS | {Path(args.output_dir).as_posix()}),
        "file_count": len(file_records),
        "files": file_records,
        "placeholder_findings": findings,
        "submission_archive": archive_path.name,
        "submission_archive_bytes": archive_path.stat().st_size,
        "submission_archive_sha256": sha256_file(archive_path),
        "archive_scope": "Journal-facing submission and provenance bundle; not a substitute for the permanent repository release.",
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
