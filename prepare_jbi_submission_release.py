#!/usr/bin/env python3
"""Build a deterministic Journal of Biogeography submission package and manifest.

Preview mode permits explicit author-controlled placeholders and reports them.
Strict mode fails if journal-facing files retain unresolved placeholders or if the
working tree contains substantive changes. Known validator-generated reports are
ignored for cleanliness because they are deterministic build products, not source.
The repository release itself should still be archived through a permanent service
such as Zenodo; this script builds the journal-facing bundle and an auditable
manifest from the frozen source commit.
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
GENERATED_STATUS_PATHS = {
    "docs/jbi_submission_validation_report.json",
}

EXPLICIT_FILES = [
    "docs/jbi_manuscript_editorial_revision.md",
    "docs/jbi_title_page_template.md",
    "docs/jbi_cover_letter_template.md",
    "docs/jbi_supporting_information_index.md",
    "docs/jbi_literature_search_provenance.md",
    "docs/jbi_author_metadata_provenance.md",
    "docs/jbi_classification_review_protocol.md",
    "docs/jbi_classification_rule_audit.md",
    "docs/jbi_taxon_image_candidate.md",
    "docs/jbi_gbif_doi_protocol.md",
    "docs/jbi_editorial_submission_check.md",
    "docs/jbi_submission_completion_checklist.md",
    "docs/jbi_author_confirmation_form.md",
    "docs/jbi_archive_release_protocol.md",
    "docs/jbi_zenodo_metadata_template.json",
    "analysis_evidence_spatial_scale.py",
    "analysis_enriched_spatial_scale.py",
    "data/resolved_inputs/global_flower_colour_review_queue_resolved.csv",
    "build_jbi_blinded_classification_review.py",
    "update_jbi_review_supporting_index.py",
    "analysis_phylogenetic_spatial_scale.R",
    "analysis_vphylomaker2_dated_spatial_scale.R",
    "build_gbif_occurrence_sample_paginated.py",
    "prepare_jbi_gbif_doi_bundle.py",
    "submit_jbi_gbif_download.py",
    "validate_jbi_submission_package.py",
    "validate_jbi_gbif_doi_bundle.py",
    "validate_jbi_author_prefill.py",
    "validate_jbi_classification_transparency.py",
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


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


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


def status_path(line: str) -> str:
    payload = line[3:].strip()
    if " -> " in payload:
        payload = payload.split(" -> ", 1)[1]
    return payload.strip('"')


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
    raw_status_lines = [line for line in git_text("status", "--porcelain").splitlines() if line]
    ignored_status_lines = [line for line in raw_status_lines if status_path(line) in GENERATED_STATUS_PATHS]
    substantive_status_lines = [line for line in raw_status_lines if status_path(line) not in GENERATED_STATUS_PATHS]
    findings = placeholder_findings(files)

    if args.strict and substantive_status_lines:
        raise SystemExit(f"Strict release requires no substantive working-tree changes: {substantive_status_lines}")
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
        "working_tree_clean_for_release": not substantive_status_lines,
        "substantive_working_tree_status": substantive_status_lines,
        "ignored_generated_working_tree_status": ignored_status_lines,
        "ignored_generated_paths": sorted(GENERATED_STATUS_PATHS),
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
