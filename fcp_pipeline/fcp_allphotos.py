"""Outcome-blind, uncapped photo opportunity and within-species FCP support ledger.

Species contribute equally to downstream comparisons. Every valid unique photo is
used for the photographic distribution, while observation/observer weighting
exposes repeated-image and recorder biases. This is *not* a population-genetics
estimator or a replacement for frozen H2 and H1 null procedures.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sqlite3
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

MORPHS = ("white", "yellow_orange", "red_pink", "blue_purple")
MIN_CLASSIFIABLE = 40  # Frozen paper threshold; not a photograph acquisition cap.
SECOND_GATE = 0.10  # Frozen H2 gate; continuous two-mode gate is a separate step.
REFERENCE_SHA = {
    "discovery": "ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4",
    "validation": "0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6",
    "prospective": "57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186",
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _parse_int(v: object, kind: str) -> int:
    if pd.isna(v) or str(v).strip() == "":
        raise ValueError(f"Missing {kind}")
    try:
        iv = int(str(v))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid {kind}: {v}") from exc
    if iv <= 0 or str(iv) != str(v).strip():
        raise ValueError(f"Invalid {kind}: {v}")
    return iv


def opportunity_frame(path: Path, expected: int = 42111) -> pd.DataFrame:
    """Allow *every* source-frame species; no >=100 restriction or photo cap."""
    df = pd.read_csv(path, dtype={"inat_taxon_id": str}, low_memory=False)
    if not {"inat_taxon_id", "species"} <= set(df):
        raise ValueError("Expected canonical inat_taxon_id and species columns")
    df["inat_taxon_id"] = df["inat_taxon_id"].map(lambda x: _parse_int(x, "taxon id"))
    df["species"] = df["species"].astype("string").str.strip()
    if df["inat_taxon_id"].duplicated().any() or df["species"].isna().any() or not df["species"].str.len().ge(3).all():
        raise ValueError("Duplicate taxon ID or missing species name")
    if df["species"].duplicated().any():
        raise ValueError("Duplicate species names; resolve taxonomy before pool assembly")
    if len(df) != expected:
        raise ValueError(f"Opportunity frame changed from {expected}: {len(df)}")
    if "after_observer_cap" in df:
        x = pd.to_numeric(df["after_observer_cap"], errors="raise")
        if x.isna().any() or (x < 0).any():
            raise ValueError("Invalid observer-cap opportunity metadata")
        # This legacy capacity is a *priority hint*, not a sampling constraint.
        df["legacy_capacity_hint"] = x.astype(int)
    return df.sort_values("inat_taxon_id", kind="stable").reset_index(drop=True)


def opportunity_schedule(frame: pd.DataFrame, shards: int = 32) -> pd.DataFrame:
    if shards <= 0:
        raise ValueError("shards must be positive")
    result = frame[[c for c in ("inat_taxon_id", "species", "legacy_capacity_hint") if c in frame]].copy()
    digests = [hashlib.sha256(f"FCP_ALL_OBSERVABLE_20261014|{ident}|{species}".encode()).hexdigest()
               for ident, species in frame[["inat_taxon_id", "species"]].itertuples(index=False, name=None)]
    result["partition"] = [int(h[:16], 16) % shards for h in digests]
    result["selection_hash"] = digests
    result["acquisition_cap"] = "none"
    result["eligibility_rule"] = "n_classifiable>=40_after_measurement"
    # Legacy photo availability is only a queue priority, never an exclusion.
    result["priority_tier"] = ("unknown" if "legacy_capacity_hint" not in result
                               else np.select([result.legacy_capacity_hint >= 40,
                                               result.legacy_capacity_hint >= 20],
                                              ["available_ge40", "available_20_39"],
                                              default="currently_lt20"))
    return result


def _bool(s: pd.Series) -> pd.Series:
    return s.fillna("").astype(str).str.lower().str.strip().isin(["1", "true", "yes", "y"])


def init_index(db: sqlite3.Connection) -> None:
    db.execute("PRAGMA journal_mode=WAL")
    db.execute("CREATE TABLE IF NOT EXISTS photos (photo_id TEXT PRIMARY KEY, taxon_id INTEGER NOT NULL, observation_id TEXT NOT NULL, observer_id TEXT NOT NULL, morph TEXT, classifiable INTEGER NOT NULL)")
    db.execute("CREATE INDEX IF NOT EXISTS idx_taxon ON photos(taxon_id)")


def ingest_photos(db: sqlite3.Connection, paths: dict[str, Path], frame: pd.DataFrame,
                  chunksize: int = 12000) -> tuple[dict, int]:
    """Stream all measured records; use SQLite to guarantee global photo-ID uniqueness.

    Input chunks may contain *any* number of photos per species or observation;
    only duplicate photo IDs are prohibited. Photo failures remain in denominator.
    """
    reference = dict(zip(frame.inat_taxon_id, frame.species.astype(str)))
    hashes = {}
    n_input = 0
    for label, path in paths.items():
        digest = sha256_file(path)
        if label in REFERENCE_SHA and digest != REFERENCE_SHA[label]:
            raise ValueError(f"Frozen {label} photo bytes differ")
        for chunk in pd.read_csv(path, chunksize=chunksize, dtype={"photo_id": str, "observation_id": str, "observer_id": str}, low_memory=False):
            required = {"inat_taxon_id", "species", "photo_id", "observation_id", "observer_id", "global_classifiable", "morph"}
            if not required <= set(chunk.columns):
                raise ValueError(f"{label}: missing measurement fields {sorted(required-set(chunk.columns))}")
            chunk["inat_taxon_id"] = chunk["inat_taxon_id"].map(lambda x: _parse_int(x, "photo taxon id"))
            if not set(chunk.inat_taxon_id).issubset(reference):
                raise ValueError("Measured taxon absent from pre-frozen opportunity frame")
            if any(str(s).strip() != reference[int(t)] for t, s in chunk[["inat_taxon_id", "species"]].itertuples(index=False, name=None)):
                raise ValueError("Measured taxon name changed; needs explicit taxonomy ledger")
            for field in ("photo_id", "observation_id", "observer_id"):
                if chunk[field].isna().any() or chunk[field].astype(str).str.strip().eq("").any():
                    raise ValueError(f"{label}: missing {field}; cannot audit independence")
            valid = _bool(chunk.global_classifiable) & chunk.morph.isin(MORPHS)
            batch = [(str(pid), int(tax), str(obs), str(user), str(col) if ok else None, int(ok))
                     for pid, tax, obs, user, col, ok in zip(chunk.photo_id, chunk.inat_taxon_id,
                         chunk.observation_id, chunk.observer_id, chunk.morph, valid)]
            try:
                with db:
                    db.executemany("INSERT INTO photos(photo_id,taxon_id,observation_id,observer_id,morph,classifiable) VALUES (?,?,?,?,?,?)", batch)
            except sqlite3.IntegrityError as exc:
                raise ValueError("Duplicate photo_id across input files; no double-counting allowed") from exc
            n_input += len(chunk)
        hashes[label] = {"sha256": digest, "rows": n_input - sum(z["rows"] for z in hashes.values())}
    return hashes, n_input


def _diversity(shares: np.ndarray) -> float:
    if shares.sum() == 0:
        return float("nan")
    p = shares / shares.sum()
    return float(1 - np.dot(p, p))


def species_ledger(db: sqlite3.Connection, frame: pd.DataFrame) -> pd.DataFrame:
    """One species = one comparable unit, regardless of photograph count."""
    base = frame[["inat_taxon_id", "species"]].copy()
    for c in ("legacy_capacity_hint",):
        if c in frame: base[c] = frame[c]
    groups = {}
    for row in db.execute("SELECT taxon_id, COUNT(*), COUNT(DISTINCT observation_id), COUNT(DISTINCT observer_id), SUM(classifiable) FROM photos GROUP BY taxon_id"):
        groups[int(row[0])] = dict(n_photos=int(row[1]), n_observations=int(row[2]), n_observers=int(row[3]), n_classifiable=int(row[4]))
    # Observation-weighted and observer-weighted composition use every photo but
    # do not treat 12 shots by one observer as 12 independent biological records.
    morph_map = {v:i for i,v in enumerate(MORPHS)}
    rows = []
    q = "SELECT observation_id, observer_id, morph, COUNT(*) FROM photos WHERE taxon_id=? AND classifiable=1 GROUP BY observation_id,observer_id,morph"
    for taxon, species in base[["inat_taxon_id", "species"]].itertuples(index=False, name=None):
        entry = groups.get(int(taxon), dict(n_photos=0, n_observations=0, n_observers=0, n_classifiable=0))
        n = entry["n_classifiable"]
        counts = np.zeros(4, dtype=float)
        obs: dict[str, np.ndarray] = defaultdict(lambda: np.zeros(4, dtype=float))
        observers: dict[str, np.ndarray] = defaultdict(lambda: np.zeros(4, dtype=float))
        for obs_id, user_id, morph, count in db.execute(q, (int(taxon),)):
            i = morph_map[morph]; counts[i] += count
            obs[str(obs_id)][i] += count
            observers[str(user_id)][i] += count
        assert abs(sum(counts) - n) < 1e-9
        if n:
            p = counts / n
            n_states = int((counts > 0).sum())
            first, second = sorted(counts, reverse=True)[:2]
            second_fraction = float(second / n)
            obs_p = np.vstack([x / x.sum() for x in obs.values()]).mean(axis=0)
            ob_p = np.vstack([x / x.sum() for x in observers.values()]).mean(axis=0)
            raw_D = _diversity(counts)
            o_D = _diversity(obs_p)
            b_D = _diversity(ob_p)
        else:
            raw_D = o_D = b_D = float("nan")
            second_fraction = float("nan");n_states=0
        # Rank is a detection label, NOT genetic polymorphism diagnosis.
        if entry["n_photos"] == 0: status="not_measured"
        elif n < MIN_CLASSIFIABLE: status="insufficient_classifiable"
        elif n_states == 1: status="single_photo_state_detected"
        else: status="multiple_photo_states_detected"
        eligible = n >= MIN_CLASSIFIABLE
        rows.append({"inat_taxon_id": int(taxon), **entry,
                     "n_classifiable_observations": len(obs),
                     "n_classifiable_observers": len(observers),
                     "n_observed_states": n_states,
                     "D_raw": raw_D if eligible else float("nan"),
                     "D_observation_equal": o_D if eligible else float("nan"),
                     "D_observer_equal": b_D if eligible else float("nan"),
                     "coarse_second_fraction": second_fraction if eligible else float("nan"),
                     "H2_coarse_0_10_pass": bool(eligible and second_fraction >= SECOND_GATE),
                     "H2_continuous_gate_checked": False,
                     "measured_status": status})
    result = base.merge(pd.DataFrame(rows), on="inat_taxon_id", how="left", validate="one_to_one")
    return result


def coverage_receipt(frame: pd.DataFrame, ledger: pd.DataFrame, hashes: dict, n_input: int) -> dict:
    n_eligible = int(ledger.n_classifiable.ge(MIN_CLASSIFIABLE).sum())
    statuses = ledger.measured_status.value_counts().to_dict()
    thresholds = (1, 5, 10, 20, 40, 60, 80, 100)
    return {
        "schema": "fcp_uncapped_photo_opportunity_and_FCP_comparison_v1",
        "status": "MEASUREMENT_LEDGER_ONLY_NOT_H2_TEST",
        "n_source_species": len(frame), "n_unique_measured_photo_ids": n_input,
        "n_D_evaluable_species": n_eligible,
        "n_H2_coarse_gate_species": int(ledger.H2_coarse_0_10_pass.sum()),
        "n_H2_continuous_confirmed_species": None,
        "legacy_opportunity_thresholds_only": {
             str(t): int(frame.legacy_capacity_hint.ge(t).sum()) for t in thresholds
        } if "legacy_capacity_hint" in frame else None,
        "measured_status": {str(k): int(v) for k, v in statuses.items()},
        "source_input_sha256": hashes,
        "all_rows_used_no_per_species_photo_cap": True,
        "species_equal_comparative_weighting": True,
        "source_frame_not_global_flora_or_FCP_prevalence_denominator": True,
        "no_genetic_or_fitness_FCP_claim": True,
        "frozen_original_H2_recomputed": False,
        "environmental_causation_identified": False,
    }


def load_inputs(spec: Iterable[str]) -> dict[str, Path]:
    out = {}
    for item in spec:
        if "=" not in item:
            raise ValueError("--measured expects label=path")
        label, path = item.split("=", 1)
        if not label or label in out or not path:
            raise ValueError("Invalid/duplicate input label")
        out[label] = Path(path)
    return out


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--opportunity", type=Path, required=True)
    p.add_argument("--measured", action="append", default=[], help="label=csv, repeat as needed")
    p.add_argument("--expected-species", type=int, default=42111)
    p.add_argument("--shards", type=int, default=32)
    p.add_argument("--outdir", type=Path, required=True)
    args = p.parse_args(argv)
    frame = opportunity_frame(args.opportunity, args.expected_species)
    args.outdir.mkdir(parents=True, exist_ok=True)
    schedule = opportunity_schedule(frame, args.shards)
    schedule.to_csv(args.outdir / "all_species_acquisition_schedule.csv", index=False)
    paths = load_inputs(args.measured)
    with tempfile.TemporaryDirectory(prefix="fcp_allphotos_") as temp:
        with sqlite3.connect(Path(temp) / "photo_ids.sqlite") as conn:
            init_index(conn)
            hashes, n = ingest_photos(conn, paths, frame)
            ledger = species_ledger(conn, frame)
    ledger.to_csv(args.outdir / "all_photo_FCP_species_ledger.csv", index=False)
    receipt = coverage_receipt(frame, ledger, hashes, n)
    receipt["opportunity_source_sha256"] = sha256_file(args.opportunity)
    (args.outdir / "all_photo_FCP_receipt.json").write_text(json.dumps(receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"source_species": len(frame), "photos": n,
                      "D_evaluable_species": receipt["n_D_evaluable_species"],
                      "outdir": str(args.outdir)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
