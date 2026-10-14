#!/usr/bin/env python3
"""Freeze the historical 4,730 species' outcome-blind identity/status roster.

This script reads only (taxon ID, species name, historical observer-capped
capacity), P100 / P500 membership and the old third-cohort selection; it
never opens photographic colour measurements, other traits or climate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

VERSION = "fcp_u100_within_species_outcome_blind_roster_v1"
FRAME_SHA256 = "5d871bd1f1190c68fd513bb527e6c93de1606dc35ff8cfff5f6187850dc382cc"
P100_SHA256 = "1473aad680fe2fa84903c5e11ee104fd0eceef828957f16c2ef1a35f9dd6993c"
P500_SHA256 = "f54d07fb2338a20a3f92808b58f0ebc948ab0d5af3cb01fb26702f861184b7a4"
THIRD_SHA256 = "16db6a2fc265f0ab20e7dae4e6f9058b907f5c19af3ddab5b6077797dd1def59"
STATUS = {
    "legacy_discovery_or_validation": 1000,
    "earlier_P500_selected_not_new_holdout": 500,
    "third_cohort_selected_500_of_499_measured": 500,
    "still_never_selected_for_high_depth": 2730,
}
FIELDS = ["inat_taxon_id", "species", "after_observer_cap"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def canonical(df: pd.DataFrame, *, require_cap=True) -> pd.DataFrame:
    wanted = FIELDS if require_cap else FIELDS[:2]
    missing = set(wanted) - set(df.columns)
    if missing:
        raise ValueError(f"Missing outcome-blind identity columns: {sorted(missing)}")
    d = df[wanted].copy()
    d["inat_taxon_id"] = pd.to_numeric(d.inat_taxon_id, errors="raise").astype("int64")
    d["species"] = d.species.astype("string").str.strip()
    if d.inat_taxon_id.duplicated().any() or d.species.duplicated().any():
        raise ValueError("Taxon ID and species name are not one-to-one")
    if d.species.isna().any() or (d.species.str.len() == 0).any():
        raise ValueError("Blank species name")
    if require_cap:
        d["after_observer_cap"] = pd.to_numeric(d.after_observer_cap, errors="raise").astype(int)
    return d.sort_values("inat_taxon_id", kind="mergesort").reset_index(drop=True)


def assert_subframe(child: pd.DataFrame, parent: pd.DataFrame, child_name: str) -> None:
    identity_parent = dict(zip(parent.inat_taxon_id, parent.species))
    for tid, name in zip(child.inat_taxon_id, child.species):
        if identity_parent.get(int(tid)) != name:
            raise ValueError(f"{child_name} identity missing/mismatched in parent: {tid} / {name}")


def make_roster(frame: pd.DataFrame, p100: pd.DataFrame, p500: pd.DataFrame, third: pd.DataFrame,
                *, expect_counts: dict[str, int] | None = None):
    # The 42,111 row source may contain other metadata, but the only columns used
    # for the U100 census are three outcome-blind identifiers/capacities.
    f = canonical(frame)
    pool = canonical(p100)
    old500 = canonical(p500)
    new500 = canonical(third, require_cap=False)
    universe = f.loc[f.after_observer_cap >= 100].copy()
    assert_subframe(pool, universe, "P100")
    assert_subframe(old500, pool, "old P500")
    assert_subframe(new500, pool, "third 500")
    if set(old500.inat_taxon_id) & set(new500.inat_taxon_id):
        raise ValueError("P500 and third-cohort selected species overlap")
    # Validate capacity from source versus once-frozen P100 selection; never
    # silently prefer one if historical opportunity changed.
    capacities = universe.set_index("inat_taxon_id").after_observer_cap
    if not (pool.set_index("inat_taxon_id").after_observer_cap == capacities.loc[pool.inat_taxon_id].to_numpy()).all():
        raise ValueError("Historic observer-capped capacity differs between sources")

    ids100 = set(map(int, pool.inat_taxon_id))
    idsP500 = set(map(int, old500.inat_taxon_id))
    idsThird = set(map(int, new500.inat_taxon_id))
    statuses = []
    for tid in universe.inat_taxon_id:
        tid = int(tid)
        if tid not in ids100: status = "legacy_discovery_or_validation"
        elif tid in idsP500: status = "earlier_P500_selected_not_new_holdout"
        elif tid in idsThird: status = "third_cohort_selected_500_of_499_measured"
        else: status = "still_never_selected_for_high_depth"
        statuses.append(status)
    universe["historical_membership"] = statuses
    universe["fresh_outcome_unexposed_candidate"] = universe.historical_membership.eq("still_never_selected_for_high_depth")
    universe["new_colour_measurement_status"] = "NOT_PERFORMED_BY_THIS_FREEZE"
    found = {name: int((universe.historical_membership == name).sum()) for name in STATUS}
    if expect_counts is not None and found != expect_counts:
        raise ValueError(f"unexpected historical cohort attrition {found} != {expect_counts}")
    if len(universe) != sum(found.values()):
        raise ValueError("cohort status coverage failure")
    return universe.reset_index(drop=True), found


def run(frame_path: Path, p100_path: Path, p500_path: Path, selected_path: Path, outdir: Path) -> dict:
    inputs = {"frame": (frame_path, FRAME_SHA256), "P100": (p100_path, P100_SHA256),
              "P500": (p500_path, P500_SHA256), "third_selected": (selected_path, THIRD_SHA256)}
    for name, (path, sha) in inputs.items():
        observed = sha256(path)
        if observed != sha:
            raise RuntimeError(f"{name} must match frozen byte SHA-256: {observed} != {sha}")
    # usecols prohibits accidentally opening picture colour, D, W, local
    # environment, pigments, phylogenies or any other outcome-derived column.
    frame = pd.read_csv(frame_path, usecols=FIELDS)
    p100 = pd.read_csv(p100_path, usecols=FIELDS)
    p500 = pd.read_csv(p500_path, usecols=FIELDS)
    third = pd.read_csv(selected_path, sep="\t", usecols=FIELDS[:2])
    if len(frame) != 42111 or len(p100) != 3730 or len(p500) != 500 or len(third) != 500:
        raise RuntimeError("historical input denominator mismatch")
    roster, found = make_roster(frame, p100, p500, third, expect_counts=STATUS)
    if len(roster) != 4730:
        raise RuntimeError(f"expected all 4730 historical U100 IDs, not {len(roster)}")
    outdir.mkdir(parents=True, exist_ok=True)
    roster_path = outdir / "all_u100_4730_species_status.csv"
    never_path = outdir / "never_selected_2730_preoutcome_species.csv"
    roster.to_csv(roster_path, index=False, lineterminator="\n")
    roster.loc[roster.fresh_outcome_unexposed_candidate, FIELDS].to_csv(never_path, index=False, lineterminator="\n")
    receipt = {
        "schema": VERSION,
        "status": "U100_IDENTITIES_FROZEN_NO_NEW_COLOUR_OR_ENVIRONMENT_OPENED",
        "historical_opportunity_species": int(len(frame)),
        "high_depth_opportunity_species": int(len(roster)),
        "historical_memberships": found,
        "distinct_unselected_hypothesis_test_candidates": found["still_never_selected_for_high_depth"],
        "historical_colour_measurement_observation": "The discovery+validation cohort (1000 species) and third group (499/500) have measured photos; original P500 (500) is historically selected and must be audited separately; never-selected 2730 are not measured here",
        "photo_colours_retrieved_by_this_step": 0,
        "new_W_tests_performed": 0,
        "new_within_species_environment_tests_performed": 0,
        "q_white_changes": False,
        "original_source_SHA256": {k: v[1] for k,v in inputs.items()},
        "output_SHA256": {"all_u100_4730_species_status.csv": sha256(roster_path),
                          "never_selected_2730_preoutcome_species.csv": sha256(never_path)},
        "source_column_allowlist": FIELDS,
        "not_measurement_status": "4730 species satisfy historic photo capacity only, NOT 4730 confirmed biological polymorphs nor 4730 observed two-mode W axes",
    }
    (outdir / "outcome_blind_roster_receipt.json").write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({"status":receipt["status"],"n":found,"sha":receipt["output_SHA256"]},ensure_ascii=False))
    return receipt


def main():
    p=argparse.ArgumentParser()
    for option in ("frame","p100","p500","third-selected"):
        p.add_argument("--"+option,required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    a=p.parse_args()
    run(a.frame,a.p100,a.p500,a.third_selected,a.outdir)

if __name__ == "__main__": main()
