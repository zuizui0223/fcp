#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RESOLVED_ROOT = ROOT / "results/rgfca_42111_anchor_resolution_full_step8d_20260911"
RESOLVED = RESOLVED_ROOT / "anchor_metadata_42111.csv.gz"
RESOLVED_RESULT = RESOLVED_ROOT / "result.json"
V1 = ROOT / "data/frozen/global_monte_carlo_species_discovery_observation_index_v1.csv.gz"
V2 = ROOT / "data/frozen/global_monte_carlo_species_discovery_v2_observation_index_v1.csv.gz"
OUT = ROOT / "results/rgfca_42111_breadth_prepixel_step8e_20260911"
N_SPECIES = 42111
TRANSPORT_GATE = 0.95


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def canonical_frame_sha256(frame: pd.DataFrame) -> str:
    """Stable content identity independent of gzip container metadata."""
    if "breadth_rank" not in frame.columns:
        raise RuntimeError("canonical Step-8E frame lacks breadth_rank")
    canonical = frame.sort_values("breadth_rank", kind="mergesort").to_csv(
        index=False,
        lineterminator="\n",
    )
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    if not RESOLVED.exists() or not RESOLVED_RESULT.exists():
        raise RuntimeError("Step 8D full-resolution result is missing")
    r = json.loads(RESOLVED_RESULT.read_text(encoding="utf-8"))
    if r.get("status") != "complete_metadata_only_full_anchor_resolution":
        raise RuntimeError("Step 8D full-resolution status is not complete")
    if r.get("species_universe") != N_SPECIES or r.get("all_species_retained") is not True:
        raise RuntimeError("Step 8D denominator drifted")
    if r.get("image_pixels_opened") is not False or r.get("flower_colour_used") is not False:
        raise RuntimeError("Step 8D violated the metadata-only firewall")

    x = pd.read_csv(RESOLVED)
    if len(x) != N_SPECIES or x["inat_taxon_id"].nunique() != N_SPECIES:
        raise RuntimeError("resolved anchor table is not exact 42,111 species")
    if x["breadth_rank"].tolist() != list(range(1, N_SPECIES + 1)):
        raise RuntimeError("breadth rank is not exact 1..42,111")

    idx_parts = []
    for source, path in (("v1", V1), ("v2", V2)):
        q = pd.read_csv(path, usecols=["cell_id", "observation_id", "photo_id", "inat_taxon_id"])
        q["source"] = source
        idx_parts.append(q)
    idx = pd.concat(idx_parts, ignore_index=True)
    for c in ["cell_id", "observation_id", "photo_id", "inat_taxon_id"]:
        idx[c] = pd.to_numeric(idx[c], errors="raise").astype(int)
    idx = idx.drop_duplicates(["observation_id", "photo_id"], keep="first")

    joined = x.merge(
        idx[["cell_id", "observation_id", "photo_id", "inat_taxon_id"]],
        on=["observation_id", "photo_id", "inat_taxon_id"],
        how="left",
        validate="one_to_one",
    )
    if joined["cell_id"].isna().any():
        n = int(joined["cell_id"].isna().sum())
        raise RuntimeError(f"{n} frozen anchors lack their original discovery cell")
    joined["cell_id"] = joined["cell_id"].astype(int)
    joined = joined.sort_values("breadth_rank", kind="mergesort").reset_index(drop=True)

    resolved = joined["status"].astype(str).eq("resolved")
    resolved_n = int(resolved.sum())
    resolved_fraction = resolved_n / float(N_SPECIES)
    gate_pass = bool(resolved_fraction >= TRANSPORT_GATE)
    candidate = joined.loc[resolved].copy().sort_values("breadth_rank", kind="mergesort").reset_index(drop=True)
    if candidate["photo_url_large"].fillna("").astype(str).str.len().eq(0).any():
        raise RuntimeError("resolved candidate contains an empty current photo URL")

    candidate_path = OUT / "breadth_source_manifest_resolved.csv.gz"
    denominator_path = OUT / "breadth_denominator_42111.csv.gz"
    candidate.to_csv(candidate_path, index=False, compression="gzip", lineterminator="\n")
    joined.to_csv(denominator_path, index=False, compression="gzip", lineterminator="\n")

    result = {
        "analysis": "rgfca_42111_breadth_prepixel_step8e",
        "status": "transport_gate_pass" if gate_pass else "transport_gate_fail",
        "species_denominator": N_SPECIES,
        "resolved_exact_anchor_species": resolved_n,
        "unresolved_anchor_species": int(N_SPECIES - resolved_n),
        "resolved_fraction": resolved_fraction,
        "required_resolved_fraction": TRANSPORT_GATE,
        "transport_evaluable": gate_pass,
        "occupied_discovery_cells_in_denominator": int(joined["cell_id"].nunique()),
        "occupied_discovery_cells_in_resolved": int(candidate["cell_id"].nunique()),
        "image_pixels_opened": False,
        "flower_colour_used": False,
        "pixel_opening_authorized": False,
        "lineage": {
            "step8d_result_sha256": sha256_file(RESOLVED_RESULT),
            "step8d_table_sha256": sha256_file(RESOLVED),
            "source_manifest_canonical_sha256": canonical_frame_sha256(candidate),
            "denominator_canonical_sha256": canonical_frame_sha256(joined),
            "source_manifest_gzip_sha256_diagnostic_only": sha256_file(candidate_path),
            "denominator_gzip_sha256_diagnostic_only": sha256_file(denominator_path),
            "canonical_hash_rule": "sort breadth_rank mergesort; pandas to_csv(index=False, lineterminator='\\n'); SHA256 UTF-8 text"
        },
        "boundary": "This gate authorizes nothing by itself. A separate canonical-content-hash-bound authorization is required before Step-8F image pixels may open."
    }
    (OUT / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (OUT / "RESULT.md").write_text(
        "# RGFCA Step 8E — pre-pixel transport gate\n\n"
        f"- denominator: **{N_SPECIES:,} species**\n"
        f"- exact anchors resolved: **{resolved_n:,}**\n"
        f"- unresolved: **{N_SPECIES-resolved_n:,}**\n"
        f"- resolved fraction: **{resolved_fraction:.6f}**\n"
        f"- required fraction: **{TRANSPORT_GATE:.2f}**\n"
        f"- transport evaluable: **{gate_pass}**\n"
        f"- occupied discovery cells retained: **{int(candidate['cell_id'].nunique())} / {int(joined['cell_id'].nunique())}**\n"
        "- lineage authorization identity: **canonical uncompressed table content hash**\n"
        "- raw gzip hashes: **diagnostic only**\n"
        "- image pixels opened: **false**\n"
        "- pixel opening authorized: **false**\n\n"
        "All 42,111 species remain in the denominator; unresolved anchors are not replaced.\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
