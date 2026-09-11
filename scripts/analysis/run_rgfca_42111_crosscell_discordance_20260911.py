#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/rgfca_42111_taxon_cell_measurement_contract_v2.json"
PAIR_RESULT = ROOT / "results/rgfca_42111_crosscell_observer_preflight_20260911/result.json"
PAIR_TABLE = ROOT / "results/rgfca_42111_crosscell_observer_preflight_20260911/crosscell_observer_disjoint_pairs.csv.gz"
MAP_RESULT = ROOT / "results/rgfca_42111_taxon_cell_measurement_step8g_20260911/result.json"
MAP_TABLE = ROOT / "results/rgfca_42111_taxon_cell_measurement_step8g_20260911/taxon_cell_measured_85337.csv.gz"
OUT = ROOT / "results/rgfca_42111_crosscell_discordance_step8g_20260911"
BIOLOGICAL = {"white", "yellow_orange", "red_pink", "blue_purple"}
EXPECTED_PAIRS = 13416


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def endpoint_join(pairs: pd.DataFrame, atlas: pd.DataFrame, side: int) -> pd.DataFrame:
    map_cols = [
        "inat_taxon_id", "cell_id", "observation_id", "photo_id", "morph",
        "measurement_status", "measurement_origin"
    ]
    m = atlas[map_cols].copy()
    rename = {c: f"map_{c}_{side}" for c in map_cols}
    m = m.rename(columns=rename)
    out = pairs.merge(
        m,
        how="left",
        left_on=["inat_taxon_id", f"cell_id_{side}", f"observation_id_{side}", f"photo_id_{side}"],
        right_on=[
            f"map_inat_taxon_id_{side}", f"map_cell_id_{side}",
            f"map_observation_id_{side}", f"map_photo_id_{side}"
        ],
        validate="many_to_one",
        indicator=f"_merge_{side}",
    )
    if not out[f"_merge_{side}"].eq("both").all():
        missing = int((~out[f"_merge_{side}"].eq("both")).sum())
        raise RuntimeError(f"{missing} frozen pair endpoints failed exact atlas join on side {side}")
    out = out.drop(columns=[f"_merge_{side}"])
    return out


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    pair_result = json.loads(PAIR_RESULT.read_text(encoding="utf-8"))
    map_result = json.loads(MAP_RESULT.read_text(encoding="utf-8"))

    if contract.get("status") != "frozen_before_any_taxon_cell_specific_image_pixel_v2":
        raise RuntimeError("v2 taxon-cell contract is not frozen")
    if pair_result.get("status") != "complete_metadata_only_observer_disjoint_pair_freeze":
        raise RuntimeError("observer-disjoint pair frame is not formally frozen")
    if int(pair_result.get("observer_disjoint_crosscell_pair_species", -1)) != EXPECTED_PAIRS:
        raise RuntimeError("frozen pair denominator drift")
    if pair_result.get("image_pixels_opened") is not False or pair_result.get("flower_colour_used") is not False:
        raise RuntimeError("pair selection was not outcome blind")
    if map_result.get("status") != "complete_85337_taxon_cell_measurement_and_postcomplete_join":
        raise RuntimeError("complete taxon-cell atlas prerequisite absent")
    if int(map_result.get("taxon_cell_rows", -1)) != 85337:
        raise RuntimeError("taxon-cell atlas denominator drift")

    pairs = pd.read_csv(PAIR_TABLE).fillna("")
    atlas = pd.read_csv(MAP_TABLE).fillna("")
    if len(pairs) != EXPECTED_PAIRS or pairs["inat_taxon_id"].nunique() != EXPECTED_PAIRS:
        raise RuntimeError("pair table is not exactly one frozen pair per eligible species")
    if (pd.to_numeric(pairs["cell_id_1"]) == pd.to_numeric(pairs["cell_id_2"])).any():
        raise RuntimeError("pair table contains same-cell endpoint")
    if (pairs["observer_id_1"].astype(str) == pairs["observer_id_2"].astype(str)).any():
        raise RuntimeError("pair table contains same-observer endpoint")
    if len(atlas) != 85337 or atlas[["inat_taxon_id", "cell_id"]].drop_duplicates().shape[0] != 85337:
        raise RuntimeError("atlas is not exactly one row per frozen taxon-cell")

    int_cols_pairs = [
        "inat_taxon_id", "cell_id_1", "cell_id_2", "observation_id_1", "observation_id_2",
        "photo_id_1", "photo_id_2"
    ]
    for c in int_cols_pairs:
        pairs[c] = pd.to_numeric(pairs[c], errors="raise").astype(int)
    for c in ["inat_taxon_id", "cell_id", "observation_id", "photo_id"]:
        atlas[c] = pd.to_numeric(atlas[c], errors="raise").astype(int)

    joined = endpoint_join(pairs, atlas, 1)
    joined = endpoint_join(joined, atlas, 2)
    joined["morph_1"] = joined["map_morph_1"].astype(str)
    joined["morph_2"] = joined["map_morph_2"].astype(str)
    joined["measurement_status_1"] = joined["map_measurement_status_1"].astype(str)
    joined["measurement_status_2"] = joined["map_measurement_status_2"].astype(str)
    joined["measurement_origin_1"] = joined["map_measurement_origin_1"].astype(str)
    joined["measurement_origin_2"] = joined["map_measurement_origin_2"].astype(str)

    c1 = joined["morph_1"].isin(BIOLOGICAL) & joined["measurement_status_1"].eq("classified_four_state_morph")
    c2 = joined["morph_2"].isin(BIOLOGICAL) & joined["measurement_status_2"].eq("classified_four_state_morph")
    both = c1 & c2
    same = both & joined["morph_1"].eq(joined["morph_2"])
    discordant = both & ~joined["morph_1"].eq(joined["morph_2"])
    joined["pair_state"] = "unclassifiable"
    joined.loc[same, "pair_state"] = "same"
    joined.loc[discordant, "pair_state"] = "discordant"
    joined["both_endpoints_classifiable"] = both

    n_both = int(both.sum())
    n_same = int(same.sum())
    n_discordant = int(discordant.sum())
    n_unclassifiable = int(EXPECTED_PAIRS - n_both)
    if n_same + n_discordant + n_unclassifiable != EXPECTED_PAIRS:
        raise RuntimeError("pair-state accounting failed")

    conditional = float(n_discordant / n_both) if n_both else None
    result = {
        "protocol": contract["protocol"],
        "analysis": "rgfca_42111_observer_disjoint_crosscell_discordance",
        "status": "complete_fixed_observer_disjoint_crosscell_coarse_colour_discordance",
        "fixed_pair_species": EXPECTED_PAIRS,
        "both_endpoints_classifiable": n_both,
        "same_coarse_colour": n_same,
        "discordant_coarse_colour": n_discordant,
        "one_or_both_unclassifiable": n_unclassifiable,
        "conditional_discordance_given_both_classifiable": conditional,
        "unconditional_proportions_over_all_fixed_pairs": {
            "same": float(n_same / EXPECTED_PAIRS),
            "discordant": float(n_discordant / EXPECTED_PAIRS),
            "unclassifiable": float(n_unclassifiable / EXPECTED_PAIRS),
        },
        "observer_disjoint_by_construction": True,
        "distinct_equal_area_cells_by_construction": True,
        "post_outcome_distance_bins_used": False,
        "claim_boundary": contract["crosscell_estimand"]["claim_boundary"],
        "hard_nonclaims": contract["hard_nonclaims"],
        "lineage": {
            "contract_sha256": sha256_file(CONTRACT),
            "pair_result_sha256": sha256_file(PAIR_RESULT),
            "pair_table_sha256": sha256_file(PAIR_TABLE),
            "map_result_sha256": sha256_file(MAP_RESULT),
            "map_table_sha256": sha256_file(MAP_TABLE),
        }
    }

    keep = [
        "species", "inat_taxon_id", "pair_hash", "cell_id_1", "cell_id_2",
        "observation_id_1", "observation_id_2", "photo_id_1", "photo_id_2",
        "observer_id_1", "observer_id_2", "cell_centroid_distance_km",
        "morph_1", "morph_2", "measurement_status_1", "measurement_status_2",
        "measurement_origin_1", "measurement_origin_2", "both_endpoints_classifiable", "pair_state"
    ]
    pair_out = joined[keep].sort_values("inat_taxon_id", kind="mergesort").reset_index(drop=True)
    pair_path = OUT / "crosscell_discordance_pairs_13416.csv.gz"
    pair_out.to_csv(pair_path, index=False, compression="gzip", lineterminator="\n")
    result["lineage"]["pair_output_sha256"] = sha256_file(pair_path)
    (OUT / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    cond_text = "NA" if conditional is None else f"{conditional:.3%}"
    (OUT / "RESULT.md").write_text(
        "# RGFCA Step 8G — observer-disjoint cross-cell coarse-colour discordance\n\n"
        f"- frozen pair species: **{EXPECTED_PAIRS:,}**\n"
        f"- both endpoints classifiable: **{n_both:,}**\n"
        f"- same coarse colour: **{n_same:,}**\n"
        f"- discordant coarse colour: **{n_discordant:,}**\n"
        f"- one/both unclassifiable: **{n_unclassifiable:,}**\n"
        f"- conditional discordance among classifiable pairs: **{cond_text}**\n"
        "- claim ceiling: **cross-cell coarse colour discordance between two geographically separated, observer-disjoint photos**.\n"
        "- this is not within-population coexistence, genetic differentiation, adaptation, or flower-colour-polymorphism prevalence.\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
