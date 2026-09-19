#!/usr/bin/env python3
"""Export a standalone disttrait release-candidate tree from the FCP monorepo.

This is an internal dry-run exporter. It deliberately does not invent licence,
authorship, citation metadata, a standalone repository URL or package-index
credentials.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / "packages" / "disttrait"

ROOT_FILES = (
    "pyproject.toml",
    "README.md",
    "CHANGELOG.md",
    "RELEASING.md",
)

COPY_DIRS = (
    "src",
    "examples",
    "benchmarks",
    "scripts",
    "fixtures",
)

STANDALONE_TESTS = (
    "test_core.py",
    "test_reliability.py",
    "test_spatial.py",
    "test_meta.py",
    "test_multivariate.py",
)

UPSTREAM_ONLY_TESTS = (
    "test_fcp_equivalence.py",
    "test_full_fcp_replay_optional.py",
    "test_nonflower_example.py",
    "test_species_conditioning_benchmark.py",
    "test_performance_surface.py",
    "test_comparator_surface.py",
    "test_model_comparator_surface.py",
    "test_sf_empirical_fixture.py",
    "test_direction_heterogeneity_receipt.py",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def package_version(pyproject: Path) -> str:
    text = pyproject.read_text(encoding="utf-8")
    match = re.search(r'^version\s*=\s*"([^"]+)"', text, flags=re.MULTILINE)
    if match is None:
        raise RuntimeError("could not parse package version")
    return match.group(1)


def copy_tree(src: Path, dst: Path) -> None:
    shutil.copytree(
        src,
        dst,
        ignore=shutil.ignore_patterns(
            "__pycache__",
            "*.pyc",
            ".pytest_cache",
            "dist",
            "build",
            "*.egg-info",
        ),
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--source-sha", default="UNKNOWN")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    out = args.out.resolve()
    if out.exists():
        if not args.force:
            raise RuntimeError(f"output already exists: {out}")
        shutil.rmtree(out)
    out.mkdir(parents=True)

    for name in ROOT_FILES:
        src = PACKAGE / name
        if not src.exists():
            raise RuntimeError(f"missing required package file: {src}")
        shutil.copy2(src, out / name)

    for name in COPY_DIRS:
        src = PACKAGE / name
        if src.exists():
            copy_tree(src, out / name)

    test_out = out / "tests"
    test_out.mkdir()
    for name in STANDALONE_TESTS:
        src = PACKAGE / "tests" / name
        if not src.exists():
            raise RuntimeError(f"missing standalone test: {src}")
        shutil.copy2(src, test_out / name)

    version = package_version(out / "pyproject.toml")
    provenance = f"""# disttrait standalone candidate provenance

This tree was exported automatically from the FCP monorepo as an internal
release-candidate dry run.

- package version: {version}
- source repository: https://github.com/zuizui0223/fcp
- source commit: {args.source_sha}
- source package path: packages/disttrait/
- scientific validation status: see the FCP frozen receipts and
  docs/DISTTRAIT_METHODS_PAPER_ARCHITECTURE_20260919.md

The standalone candidate intentionally excludes monorepo regression tests that
read frozen FCP/result receipts by relative repository paths. Those upstream
tests remain part of the source-repository validation contract.

No claim is made that this candidate is ready for public release until software
licence, authorship/maintainer metadata, citation metadata and repository
location are frozen.
"""
    (out / "PROVENANCE.md").write_text(provenance, encoding="utf-8")

    pending = """# Public-release metadata still pending

The following items require an explicit ownership/release decision and are not
filled by automation:

- software licence;
- final package authors and maintainer ordering;
- final CITATION.cff;
- standalone repository name and visibility;
- repository URLs in pyproject.toml;
- release tag and archival DOI;
- optional TestPyPI/PyPI publication.

This file is expected to disappear or become a completed release record before
the first public standalone release.
"""
    (out / "RELEASE_METADATA_PENDING.md").write_text(pending, encoding="utf-8")

    files = sorted(
        p for p in out.rglob("*")
        if p.is_file()
    )
    manifest = {
        "schema": "disttrait_standalone_candidate_manifest_v1",
        "package_version": version,
        "source_repository": "zuizui0223/fcp",
        "source_commit": args.source_sha,
        "release_ready": False,
        "reason_not_release_ready": [
            "software licence not frozen",
            "author/maintainer metadata not frozen",
            "CITATION.cff not frozen",
            "standalone repository location not frozen",
            "archival DOI not available",
        ],
        "standalone_tests": list(STANDALONE_TESTS),
        "upstream_monorepo_tests_not_exported": list(UPSTREAM_ONLY_TESTS),
        "files": {
            str(p.relative_to(out)): sha256(p)
            for p in files
        },
    }
    (out / "standalone_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print(json.dumps({
        "status": "PASS",
        "version": version,
        "files": len(files) + 1,
        "standalone_tests": len(STANDALONE_TESTS),
        "release_ready": False,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
