#!/usr/bin/env python3
"""Recover the pre-measurement 6 x 200 Chapter 1 photo manifest from a worktree.

This scanner is designed for the interrupted local Work-mode state where the acquired
1,200-photo CSV exists somewhere under the worktree but was never pushed. It searches
CSV/TSV files, identifies files that satisfy the frozen 6 species x 200 photographs
contract, and writes a canonical source manifest only when the candidate ID set is
unambiguous.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys
from typing import Iterable

import pandas as pd

from fcp_pipeline.photo_split import find_outcome_columns

SPECIES_CANDIDATES = (
    "species",
    "taxon_name",
    "scientific_name",
    "accepted_species",
    "taxon",
)
PHOTO_ID_CANDIDATES = (
    "photo_id",
    "image_id",
    "media_id",
    "inat_photo_id",
    "inaturalist_photo_id",
)
SKIP_DIRS = {
    ".git",
    ".venv",
    "venv",
    "node_modules",
    "__pycache__",
}


def _norm(name: object) -> str:
    return str(name).strip().lower().replace("-", "_").replace(" ", "_")


def _resolve(columns: Iterable[object], candidates: tuple[str, ...]) -> str | None:
    mapping = {_norm(c): str(c) for c in columns}
    for candidate in candidates:
        if candidate in mapping:
            return mapping[candidate]
    return None


def _read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix in {".tsv", ".tab"}:
        return pd.read_csv(path, sep="\t")
    raise ValueError(f"unsupported table type: {path}")


def canonical_pair_hash(frame: pd.DataFrame, species_col: str, photo_id_col: str) -> str:
    pairs = sorted(
        (str(s).strip(), str(p).strip())
        for s, p in frame[[species_col, photo_id_col]].itertuples(index=False, name=None)
    )
    payload = "\n".join(f"{s}\t{p}" for s, p in pairs).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def inspect_candidate(path: Path) -> dict[str, object]:
    result: dict[str, object] = {
        "path": str(path),
        "status": "rejected",
        "reason": None,
    }
    try:
        frame = _read_table(path)
    except Exception as exc:
        result["reason"] = f"read_error:{type(exc).__name__}:{exc}"
        return result

    species_col = _resolve(frame.columns, SPECIES_CANDIDATES)
    photo_id_col = _resolve(frame.columns, PHOTO_ID_CANDIDATES)
    result["rows"] = int(len(frame))
    result["species_col"] = species_col
    result["photo_id_col"] = photo_id_col
    result["outcome_columns"] = find_outcome_columns(frame.columns)

    if species_col is None:
        result["reason"] = "no_species_column"
        return result
    if photo_id_col is None:
        result["reason"] = "no_photo_id_column"
        return result
    if len(frame) != 1200:
        result["reason"] = f"wrong_row_count:{len(frame)}"
        return result
    if result["outcome_columns"]:
        result["reason"] = "contains_downstream_outcome_columns"
        return result
    if frame[species_col].isna().any() or frame[photo_id_col].isna().any():
        result["reason"] = "missing_species_or_photo_id"
        return result

    species = frame[species_col].astype(str).str.strip()
    photo_ids = frame[photo_id_col].astype(str).str.strip()
    if (species == "").any() or (photo_ids == "").any():
        result["reason"] = "blank_species_or_photo_id"
        return result
    if photo_ids.duplicated().any():
        result["reason"] = "duplicate_photo_id"
        return result

    counts = species.value_counts().sort_index()
    result["species_count"] = int(len(counts))
    result["per_species_counts"] = {str(k): int(v) for k, v in counts.items()}
    if len(counts) != 6:
        result["reason"] = f"wrong_species_count:{len(counts)}"
        return result
    if not (counts == 200).all():
        result["reason"] = "not_200_per_species"
        return result

    result["pair_hash"] = canonical_pair_hash(frame, species_col, photo_id_col)
    result["status"] = "eligible"
    result["reason"] = None
    return result


def iter_tables(root: Path, extra_paths: list[Path]) -> list[Path]:
    paths: set[Path] = set()
    for extra in extra_paths:
        if extra.is_file() and extra.suffix.lower() in {".csv", ".tsv", ".tab"}:
            paths.add(extra.resolve())
    for path in root.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in {".csv", ".tsv", ".tab"}:
            continue
        try:
            rel = path.relative_to(root)
        except ValueError:
            rel = path
        if any(part in SKIP_DIRS for part in rel.parts):
            continue
        paths.add(path.resolve())
    return sorted(paths, key=lambda p: str(p).lower())


def choose_candidate(results: list[dict[str, object]]) -> tuple[dict[str, object], list[dict[str, object]]]:
    eligible = [r for r in results if r.get("status") == "eligible"]
    if not eligible:
        raise RuntimeError("no eligible 6x200 pre-measurement photo manifest found")

    by_hash: dict[str, list[dict[str, object]]] = {}
    for item in eligible:
        by_hash.setdefault(str(item["pair_hash"]), []).append(item)
    if len(by_hash) > 1:
        summary = "; ".join(
            f"{digest[:12]}: {[x['path'] for x in items]}" for digest, items in by_hash.items()
        )
        raise RuntimeError(
            "multiple conflicting eligible manifests found; refuse to guess. " + summary
        )

    duplicates = next(iter(by_hash.values()))
    # Identical ID sets are equivalent for the frozen split. Prefer the shallowest path,
    # then the lexicographically first path to keep recovery deterministic.
    selected = sorted(
        duplicates,
        key=lambda r: (len(Path(str(r["path"])).parts), str(r["path"]).lower()),
    )[0]
    return selected, duplicates


def canonicalize_source(path: Path, species_col: str, photo_id_col: str) -> pd.DataFrame:
    frame = _read_table(path).copy()
    if species_col != "species":
        frame = frame.rename(columns={species_col: "species"})
    if photo_id_col != "photo_id":
        frame = frame.rename(columns={photo_id_col: "photo_id"})
    frame["species"] = frame["species"].astype(str).str.strip()
    frame["photo_id"] = frame["photo_id"].astype(str).str.strip()
    # Keep all acquisition metadata, but place the two split-critical columns first.
    remaining = [c for c in frame.columns if c not in {"species", "photo_id"}]
    return frame[["species", "photo_id", *remaining]].sort_values(
        ["species", "photo_id"], kind="mergesort"
    ).reset_index(drop=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--extra", type=Path, action="append", default=[])
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/frozen/jbi_ch1_photo_source_manifest.csv"),
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=Path("docs/supporting/jbi_ch1_photo_source_recovery_report.json"),
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    root = args.root.resolve()
    candidates = iter_tables(root, args.extra)
    results = [inspect_candidate(path) for path in candidates]

    report: dict[str, object] = {
        "root": str(root),
        "tables_scanned": len(candidates),
        "eligible_count": sum(r.get("status") == "eligible" for r in results),
        "results": results,
    }

    try:
        selected, duplicates = choose_candidate(results)
    except RuntimeError as exc:
        report["status"] = "blocked"
        report["error"] = str(exc)
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(str(exc), file=sys.stderr)
        return 2

    report["status"] = "found"
    report["selected"] = selected
    report["equivalent_duplicate_paths"] = [item["path"] for item in duplicates]

    if not args.dry_run:
        selected_path = Path(str(selected["path"]))
        frame = canonicalize_source(
            selected_path,
            str(selected["species_col"]),
            str(selected["photo_id_col"]),
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(args.output, index=False, lineterminator="\n")
        report["canonical_output"] = str(args.output)
        report["canonical_pair_hash"] = canonical_pair_hash(frame, "species", "photo_id")

    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
