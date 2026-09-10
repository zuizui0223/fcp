#!/usr/bin/env python3
"""Validate the assembled FCP polymorphism paper v0.1 science bundle."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "paper" / "polymorphism_v0_1"
MANIFEST = PAPER / "manifest.json"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    m = json.loads(MANIFEST.read_text(encoding="utf-8"))
    if m["paper_version"] != "v0.1-post-step9":
        raise RuntimeError("paper version drift")

    manuscript_check = json.loads((PAPER / m["manuscript_check"]).read_text(encoding="utf-8"))
    if manuscript_check.get("status") != "PASS":
        raise RuntimeError("manuscript contract has not passed")

    required_text = [
        m["manuscript"],
        m["manuscript_check"],
        m["figure_qc"],
        m["figure_captions"],
        m["human_readable_number_ledger"],
        m["number_ledger"],
        m["figure_manifest"],
    ]
    for rel in required_text:
        path = PAPER / rel
        if not path.exists() or path.stat().st_size == 0:
            raise RuntimeError(f"missing/empty paper surface: {rel}")

    number_ledger = json.loads((PAPER / m["number_ledger"]).read_text(encoding="utf-8"))
    if number_ledger.get("paper_version") != m["paper_version"]:
        raise RuntimeError("number ledger version mismatch")

    figure_manifest = json.loads((PAPER / m["figure_manifest"]).read_text(encoding="utf-8"))
    if figure_manifest.get("paper_version") != m["paper_version"]:
        raise RuntimeError("figure manifest version mismatch")

    declared_outputs = {entry["path"]: entry for entry in figure_manifest["outputs"]}
    checked_figures = []
    for rel in m["main_figures"]:
        path = PAPER / rel
        root_rel = str(path.relative_to(ROOT))
        if not path.exists() or path.stat().st_size < 10_000:
            raise RuntimeError(f"main figure missing/too small: {rel}")
        if root_rel not in declared_outputs:
            raise RuntimeError(f"main figure absent from figure manifest: {rel}")
        expected = declared_outputs[root_rel]["sha256"]
        actual = sha256(path)
        if expected != actual:
            raise RuntimeError(f"main figure hash drift: {rel}")
        checked_figures.append({"path": rel, "bytes": path.stat().st_size, "sha256": actual})

    for name in m["figure_data"]:
        path = PAPER / "figure_data" / name
        if not path.exists() or path.stat().st_size == 0:
            raise RuntimeError(f"figure-data file missing/empty: {name}")

    for rel in m["source_receipts"]:
        path = ROOT / rel
        if not path.exists() or path.stat().st_size == 0:
            raise RuntimeError(f"source receipt missing/empty: {rel}")

    blockers = list(m.get("submission_blockers", []))
    receipt = {
        "paper_version": m["paper_version"],
        "science_bundle_ready": True,
        "submission_ready": len(blockers) == 0,
        "main_text_word_count": manuscript_check["main_text_word_count"],
        "main_figures_checked": checked_figures,
        "figure_data_files_checked": len(m["figure_data"]),
        "source_receipts_checked": len(m["source_receipts"]),
        "remaining_submission_blockers": blockers,
        "status": "PASS",
    }
    (PAPER / "BUNDLE_READINESS.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
