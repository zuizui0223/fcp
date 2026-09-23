#!/usr/bin/env python3
"""Build FCP v2 biological Pass-B firewall after a frozen technical seal.

This script may read species/coordinate metadata only to create a sealed join key.
Biological workers receive no species, taxon, coordinates, observer or source photo ID.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

EXPECTED_ROWS = 40_000
EXPECTED_SPECIES = 400
EXPECTED_HEAVY = 8_000
BATCHES = 2
SHARDS = 32
PARTS = 4
ID_SALT = "FCP_V2_MEASUREMENT_ID_20260923"

WORKER_FIELDS = (
    "measurement_id",
    "image_filename",
    "photo_license",
    "heavy_counterfactual",
    "expected_source_sha256",
)
ACQUISITION_FIELDS = (
    "measurement_id",
    "image_filename",
    "photo_url_large",
    "photo_license",
    "heavy_counterfactual",
    "expected_source_sha256",
)


def measurement_id(photo_id: int) -> str:
    raw = f"{ID_SALT}|{int(photo_id)}".encode("utf-8")
    return "FCPV2-" + hashlib.sha256(raw).hexdigest()[:24].upper()


def slot(identifier: str) -> tuple[int, int, int]:
    digest = hashlib.sha256(identifier.encode("utf-8")).digest()
    return digest[0] % BATCHES, digest[1] % SHARDS, digest[2] % PARTS


def _bool_series(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.strip().str.casefold().isin(
        {"1", "true", "yes"}
    )


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--authorized-metadata", type=Path, required=True)
    p.add_argument("--technical-table", type=Path, required=True)
    p.add_argument("--technical-summary", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()

    summary = json.loads(args.technical_summary.read_text(encoding="utf-8"))
    if summary.get("status") != "RESPONSE_BLIND_TECHNICAL_SEAL_FROZEN":
        raise RuntimeError("technical seal status is not frozen")
    for key in (
        "biological_outcomes_opened",
        "species_opened",
        "morph_opened",
        "q_white_or_W_opened",
    ):
        if summary.get(key) is not False:
            raise RuntimeError(f"technical seal firewall opened: {key}")

    metadata = pd.read_csv(args.authorized_metadata, compression="gzip")
    technical = pd.read_csv(
        args.technical_table,
        compression="gzip",
        dtype={"measurement_id": str},
        low_memory=False,
    )
    if len(metadata) != EXPECTED_ROWS or len(technical) != EXPECTED_ROWS:
        raise RuntimeError(
            f"Pass-B denominator drift: metadata={len(metadata)} technical={len(technical)}"
        )

    required_meta = {
        "panel",
        "inat_taxon_id",
        "species",
        "observation_id",
        "photo_id",
        "photo_url_large",
        "photo_license",
        "latitude",
        "longitude",
        "observer_id",
    }
    missing = required_meta - set(metadata.columns)
    if missing:
        raise RuntimeError(f"authorized metadata missing fields: {sorted(missing)}")

    metadata = metadata.copy()
    metadata["measurement_id"] = [
        measurement_id(v) for v in metadata["photo_id"].astype(int)
    ]
    if metadata["measurement_id"].nunique() != EXPECTED_ROWS:
        raise RuntimeError("Pass-B measurement IDs are not unique")
    if technical["measurement_id"].astype(str).nunique() != EXPECTED_ROWS:
        raise RuntimeError("technical measurement IDs are not unique")
    if set(metadata["measurement_id"].astype(str)) != set(
        technical["measurement_id"].astype(str)
    ):
        raise RuntimeError("metadata and technical seal denominators differ")

    required_technical = {
        "measurement_id",
        "source_image_sha256",
        "heavy_counterfactual",
        "technical_metrics_available",
    }
    missing = required_technical - set(technical.columns)
    if missing:
        raise RuntimeError(f"technical table missing fields: {sorted(missing)}")
    if technical["source_image_sha256"].fillna("").astype(str).str.len().lt(64).any():
        raise RuntimeError("technical table contains missing source SHA256")

    technical = technical.copy()
    technical["heavy_counterfactual"] = _bool_series(
        technical["heavy_counterfactual"]
    )
    if int(technical["heavy_counterfactual"].sum()) != EXPECTED_HEAVY:
        raise RuntimeError("heavy-counterfactual denominator drift")

    joined = metadata.merge(
        technical[
            [
                "measurement_id",
                "source_image_sha256",
                "heavy_counterfactual",
            ]
        ],
        on="measurement_id",
        how="inner",
        validate="one_to_one",
    )
    if len(joined) != EXPECTED_ROWS:
        raise RuntimeError("Pass-B metadata/technical merge lost rows")

    joined["image_filename"] = joined["measurement_id"].astype(str) + ".img"
    slots = [slot(v) for v in joined["measurement_id"].astype(str)]
    joined["batch"] = [v[0] for v in slots]
    joined["shard"] = [v[1] for v in slots]
    joined["part"] = [v[2] for v in slots]

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)

    # Sealed join key is withheld from workers until every biological partition
    # has terminal coverage.
    join_columns = [
        "measurement_id",
        "panel",
        "inat_taxon_id",
        "species",
        "observation_id",
        "photo_id",
        "latitude",
        "longitude",
        "observer_id",
        "heavy_counterfactual",
    ]
    if "observer" in joined.columns:
        join_columns.append("observer")
    if "observed_on" in joined.columns:
        join_columns.append("observed_on")
    joined[join_columns].to_csv(
        out / "sealed_biological_join_key.csv",
        index=False,
        lineterminator="\n",
    )

    partition_counts = []
    for batch in range(BATCHES):
        for shard in range(SHARDS):
            for part in range(PARTS):
                sub = joined.loc[
                    (joined["batch"] == batch)
                    & (joined["shard"] == shard)
                    & (joined["part"] == part)
                ].copy()
                root = out / "partitions" / f"b{batch}_s{shard}_p{part}"
                root.mkdir(parents=True, exist_ok=True)
                worker = pd.DataFrame(
                    {
                        "measurement_id": sub["measurement_id"].astype(str),
                        "image_filename": sub["image_filename"].astype(str),
                        "photo_license": sub["photo_license"].astype(str),
                        "heavy_counterfactual": sub[
                            "heavy_counterfactual"
                        ].astype(bool),
                        "expected_source_sha256": sub[
                            "source_image_sha256"
                        ].astype(str),
                    }
                )
                acquisition = worker.copy()
                acquisition.insert(
                    2,
                    "photo_url_large",
                    sub["photo_url_large"].astype(str).to_numpy(),
                )
                worker.to_csv(
                    root / "worker_manifest.csv",
                    index=False,
                    lineterminator="\n",
                )
                acquisition.to_csv(
                    root / "acquisition_key.csv",
                    index=False,
                    lineterminator="\n",
                )
                partition_counts.append(
                    {
                        "batch": batch,
                        "shard": shard,
                        "part": part,
                        "rows": int(len(sub)),
                        "heavy_rows": int(
                            sub["heavy_counterfactual"].astype(bool).sum()
                        ),
                    }
                )

    manifest = {
        "schema": "fcp_v2_biological_firewall_v1",
        "status": "PASS_B_FIREWALL_FROZEN_AFTER_TECHNICAL_SEAL",
        "rows": EXPECTED_ROWS,
        "species_sealed": EXPECTED_SPECIES,
        "heavy_counterfactual_rows": EXPECTED_HEAVY,
        "partition_receipts_expected": BATCHES * SHARDS * PARTS,
        "worker_fields": list(WORKER_FIELDS),
        "acquisition_fields": list(ACQUISITION_FIELDS),
        "worker_contains_species": False,
        "worker_contains_taxon_id": False,
        "worker_contains_photo_id": False,
        "worker_contains_coordinates": False,
        "worker_contains_observer": False,
        "worker_contains_prior_D": False,
        "worker_contains_prior_W": False,
        "source_sha_match_required": True,
        "species_join_opened": False,
        "biological_outcomes_opened": False,
        "replacement_allowed": False,
        "partition_counts": partition_counts,
    }
    (out / "firewall_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
