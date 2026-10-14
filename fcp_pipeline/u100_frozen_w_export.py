"""Export per-species W_i without refitting the prospectively frozen FCP H2 modes.

This bridge consumes the *original* prospective 0.10 delta-vector artifact.
It does not re-open photographs or recalculate the structured null distribution.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

PALETTE = ("white", "yellow", "orange", "red", "pink", "magenta", "purple", "blue", "bronze")
TARGET = np.array([1.0] + [-1.0 / 8.0] * 8, dtype=float)
TARGET /= np.linalg.norm(TARGET)
FROZEN_VECTOR_GIT_BLOB = "4e5e0edbeb38e5c6ee5d654a272fc62e7849fbff"
FROZEN_N = 158
FROZEN_W = 0.5172457461053418


def blob_sha1(path: Path) -> str:
    """Calculate Git's SHA-1 *blob* identity, not bare file SHA-1."""
    raw = path.read_bytes()
    return hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()


def export_rows(vectors: pd.DataFrame, frame: pd.DataFrame) -> tuple[pd.DataFrame, float]:
    required = {"species", "delta_norm", "n_classifiable", "H1_coarse_second_fraction",
                "continuous_minor_fraction"} | {f"delta_{k}" for k in PALETTE}
    missing = sorted(required - set(vectors))
    if missing:
        raise ValueError(f"H2 vector columns missing: {missing}")
    if not {"inat_taxon_id", "species"}.issubset(frame):
        raise ValueError("U100 must contain inat_taxon_id and species")
    if frame["inat_taxon_id"].duplicated().any() or frame["species"].duplicated().any():
        raise ValueError("U100 duplicate taxon ID or species name")
    if vectors["species"].duplicated().any() or vectors["species"].isna().any():
        raise ValueError("duplicate/missing H2 vector species")
    values = vectors[[f"delta_{k}" for k in PALETTE]].apply(pd.to_numeric, errors="raise").to_numpy(dtype=float)
    norms = np.linalg.norm(values, axis=1)
    expected_norms = pd.to_numeric(vectors["delta_norm"], errors="raise").to_numpy(float)
    ncls = pd.to_numeric(vectors["n_classifiable"], errors="raise").to_numpy(float)
    second = pd.to_numeric(vectors["H1_coarse_second_fraction"], errors="raise").to_numpy(float)
    minor = pd.to_numeric(vectors["continuous_minor_fraction"], errors="raise").to_numpy(float)
    if (not np.isfinite(values).all() or not np.isfinite(norms).all()
            or (norms <= 0).any() or not np.allclose(norms, expected_norms, rtol=1e-12, atol=1e-12)):
        raise ValueError("frozen H2 vector/norm inconsistency")
    if (not np.isfinite(ncls).all() or not np.isfinite(second).all() or not np.isfinite(minor).all()
            or (ncls < 40).any() or (second < 0.1 - 1e-12).any() or (minor < 0.1 - 1e-12).any()):
        raise ValueError("H2 primary 0.10 qualification not met")
    w = np.square((values / norms[:, None]) @ TARGET)
    if (not np.isfinite(w).all()) or ((w < 0) | (w > 1 + 1e-12)).any():
        raise ValueError("W projection out of range")
    names = frame[["inat_taxon_id", "species"]].copy()
    joined = vectors[["species"]].merge(names, on="species", how="left", validate="one_to_one", indicator=True)
    if not joined["_merge"].eq("both").all():
        missing_names = joined.loc[~joined["_merge"].eq("both"), "species"].tolist()
        raise ValueError(f"H2 species missing in U100: {missing_names[:8]}")
    out = joined[["inat_taxon_id"]].copy()
    out["cohort"] = "prospective"
    out["w_alignment_sq"] = w
    out["method"] = "frozen_study_specific_H2"
    out = out.sort_values("inat_taxon_id", kind="stable").reset_index(drop=True)
    return out, float(np.mean(w))


def verified_export(u100: Path, vectors: Path, out: Path) -> dict:
    if blob_sha1(vectors) != FROZEN_VECTOR_GIT_BLOB:
        raise ValueError("prospective delta-vector Git blob identity mismatch")
    frame = pd.read_csv(u100, low_memory=False)
    if len(frame) != 4730:
        raise ValueError("U100 opportunity frame must contain exactly 4,730 species")
    rows, mean_w = export_rows(pd.read_csv(vectors, low_memory=False), frame)
    if len(rows) != FROZEN_N or not np.isclose(mean_w, FROZEN_W, atol=1e-12, rtol=0):
        raise ValueError("per-species W does not reproduce the prospective H2 frozen result")
    out.parent.mkdir(parents=True, exist_ok=True)
    rows.to_csv(out, index=False, float_format="%.17g")
    result = {
        "status": "validated_reuse_of_frozen_prospective_vector_rows",
        "n_species": len(rows),
        "cohort": "prospective",
        "observed_equal_species_W": mean_w,
        "authoritative_observed_W": FROZEN_W,
        "vector_git_blob": FROZEN_VECTOR_GIT_BLOB,
        "vector_sha256": hashlib.sha256(vectors.read_bytes()).hexdigest(),
        "u100_sha256": hashlib.sha256(u100.read_bytes()).hexdigest(),
        "output_sha256": hashlib.sha256(out.read_bytes()).hexdigest(),
        "H2_refit": False,
        "structured_null_recomputed": False,
        "individual_species_significance_tested": False,
        "original_H2_decisions_changed": False,
    }
    out.with_suffix(".audit.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--u100", type=Path, required=True)
    p.add_argument("--vectors", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args(argv)
    report = verified_export(a.u100, a.vectors, a.out)
    print(json.dumps({"exported_W_species": report["n_species"],
                      "W": report["observed_equal_species_W"], "refitted_H2": False}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
