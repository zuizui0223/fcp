#!/usr/bin/env python3
"""Freeze outcome-blind FCP v2 species queues from identity/capacity metadata only."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

SALT = "FCP_V2_MEASUREMENT_VALIDITY_20260923_V1"
QUEUE_PER_PANEL = 300
TARGET_PER_PANEL = 200
EXPECTED_LEGACY = 1000
EXPECTED_P500 = 500
EXPECTED_THIRD = 500
EXPECTED_U100 = 4730


def _read(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
    return pd.read_csv(path, sep=sep)


def _identity(frame: pd.DataFrame, label: str) -> pd.DataFrame:
    needed = {"inat_taxon_id", "species"}
    missing = needed - set(frame.columns)
    if missing:
        raise ValueError(f"{label} missing identity fields: {sorted(missing)}")
    out = frame[["inat_taxon_id", "species"]].copy()
    out["inat_taxon_id"] = pd.to_numeric(out["inat_taxon_id"], errors="raise").astype("int64")
    out["species"] = out["species"].astype(str).str.strip()
    if out["species"].eq("").any():
        raise ValueError(f"{label} contains blank species")
    out = out.drop_duplicates()
    if out["inat_taxon_id"].duplicated().any():
        raise ValueError(f"{label} has taxon id mapped to multiple rows")
    if out["species"].duplicated().any():
        raise ValueError(f"{label} has species name mapped to multiple rows")
    return out.sort_values(["inat_taxon_id", "species"], kind="mergesort").reset_index(drop=True)


def _hash(panel: str, taxon_id: int, species: str) -> str:
    raw = f"{SALT}|{panel}|{int(taxon_id)}|{species}".encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _queue(frame: pd.DataFrame, panel: str) -> pd.DataFrame:
    x = frame.copy()
    x["selection_hash"] = [
        _hash(panel, taxon_id, species)
        for taxon_id, species in zip(x["inat_taxon_id"], x["species"], strict=True)
    ]
    x = x.sort_values(
        ["selection_hash", "inat_taxon_id", "species"], kind="mergesort"
    ).head(QUEUE_PER_PANEL).reset_index(drop=True)
    if len(x) != QUEUE_PER_PANEL:
        raise ValueError(f"panel {panel} has only {len(x)} queue candidates")
    x.insert(0, "queue_rank", range(1, len(x) + 1))
    x.insert(1, "panel", panel)
    return x


def _csv_sha(frame: pd.DataFrame) -> str:
    payload = frame.to_csv(index=False, lineterminator="\n").encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--u100", type=Path, required=True)
    p.add_argument("--legacy", type=Path, required=True)
    p.add_argument("--p500", type=Path, required=True)
    p.add_argument("--third", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()

    u100_raw = _read(args.u100)
    if not {"inat_taxon_id", "species", "after_observer_cap"}.issubset(u100_raw.columns):
        raise ValueError("U100 input lacks identity/capacity fields")
    u100 = u100_raw[["inat_taxon_id", "species", "after_observer_cap"]].copy()
    u100["inat_taxon_id"] = pd.to_numeric(u100["inat_taxon_id"], errors="raise").astype("int64")
    u100["species"] = u100["species"].astype(str).str.strip()
    u100["after_observer_cap"] = pd.to_numeric(
        u100["after_observer_cap"], errors="raise"
    ).astype("int64")
    u100 = u100.loc[u100["after_observer_cap"] >= 100].drop_duplicates(
        ["inat_taxon_id", "species"]
    )
    if len(u100) != EXPECTED_U100:
        raise ValueError(f"U100 denominator changed: {len(u100)} != {EXPECTED_U100}")
    if u100["inat_taxon_id"].duplicated().any() or u100["species"].duplicated().any():
        raise ValueError("U100 identity is not one-to-one")

    legacy = _identity(_read(args.legacy), "legacy")
    p500 = _identity(_read(args.p500), "p500")
    third = _identity(_read(args.third), "third")

    if len(legacy) != EXPECTED_LEGACY:
        raise ValueError(f"legacy identity count changed: {len(legacy)}")
    if len(p500) != EXPECTED_P500:
        raise ValueError(f"P500 identity count changed: {len(p500)}")
    if len(third) != EXPECTED_THIRD:
        raise ValueError(f"third identity count changed: {len(third)}")

    def keys(frame: pd.DataFrame) -> set[tuple[int, str]]:
        return set(zip(frame["inat_taxon_id"].astype(int), frame["species"].astype(str), strict=True))

    legacy_keys = keys(legacy)
    p500_keys = keys(p500)
    third_keys = keys(third)
    if legacy_keys & p500_keys or legacy_keys & third_keys or p500_keys & third_keys:
        raise ValueError("previous cohort identities overlap unexpectedly")

    used = pd.concat([legacy, p500, third], ignore_index=True)
    if len(used) != EXPECTED_LEGACY + EXPECTED_P500 + EXPECTED_THIRD:
        raise ValueError("used identity denominator drift")

    paired = u100.merge(
        legacy, on=["inat_taxon_id", "species"], how="inner", validate="one_to_one"
    )
    if len(paired) < QUEUE_PER_PANEL:
        raise ValueError("fewer than 300 current-U100 legacy species available for Panel P")

    used_ids = set(used["inat_taxon_id"].astype(int))
    used_names = set(used["species"].astype(str))
    novel = u100.loc[
        ~u100["inat_taxon_id"].isin(used_ids)
        & ~u100["species"].isin(used_names)
    ].copy()
    if len(novel) < QUEUE_PER_PANEL:
        raise ValueError("fewer than 300 novel U100 species available for Panel N")

    paired_q = _queue(paired, "P")
    novel_q = _queue(novel, "N")

    if set(paired_q["inat_taxon_id"]) & set(novel_q["inat_taxon_id"]):
        raise ValueError("Panel P/N taxon overlap")
    if set(paired_q["species"]) & set(novel_q["species"]):
        raise ValueError("Panel P/N species overlap")

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    paired_path = out / "panel_P_candidate_queue.csv"
    novel_path = out / "panel_N_candidate_queue.csv"
    paired_q.to_csv(paired_path, index=False, lineterminator="\n")
    novel_q.to_csv(novel_path, index=False, lineterminator="\n")

    receipt = {
        "schema": "fcp_v2_species_queue_v1",
        "date_jst": "2026-09-23",
        "status": "OUTCOME_BLIND_METADATA_QUEUE_FROZEN",
        "salt": SALT,
        "selection_inputs_used": ["inat_taxon_id", "species", "after_observer_cap"],
        "biological_outcomes_read": False,
        "image_pixels_opened": False,
        "queue_per_panel": QUEUE_PER_PANEL,
        "terminal_target_per_panel_after_fresh_metadata": TARGET_PER_PANEL,
        "input_counts": {
            "U100": int(len(u100)),
            "legacy": int(len(legacy)),
            "P500": int(len(p500)),
            "third": int(len(third)),
            "used_total": int(len(used)),
        },
        "candidate_counts": {
            "paired_current_U100": int(len(paired)),
            "novel_current_U100_after_all_used_exclusions": int(len(novel)),
        },
        "queue_counts": {"P": int(len(paired_q)), "N": int(len(novel_q))},
        "queue_sha256": {
            "P": _csv_sha(paired_q),
            "N": _csv_sha(novel_q),
        },
        "fresh_metadata_rule": (
            "Within each frozen queue, retrieve metadata only in queue order and "
            "freeze the first 200 species with exactly 100 fresh eligible photo IDs "
            "after excluding every previously used photo ID. Stop before pixel opening."
        ),
        "no_rescue": [
            "no biological outcome may enter queue construction",
            "no queue reordering after fresh metadata retrieval",
            "no expansion beyond rank 300 without a new pre-pixel protocol",
            "no species replacement after pixel opening",
        ],
    }
    (out / "result.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
