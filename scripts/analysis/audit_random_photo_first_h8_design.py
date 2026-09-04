#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/random_photo_first_h8_design_audit_contract_v1.json"
H7_MANIFEST = ROOT / "docs/supporting/random_photo_first_h7_fresh_metadata_manifest_v1.json"
H7_AUDIT = ROOT / "data/frozen/random_photo_first_h7_fresh_metadata_target_audit_v1.csv"
OUT_JSON = ROOT / "docs/supporting/random_photo_first_h8_design_audit_result_v1.json"
OUT_SPECIES = ROOT / "data/frozen/random_photo_first_h8_design_species_v1.csv"
OUT_TARGETS = ROOT / "data/frozen/random_photo_first_h8_design_species_cells_v1.csv"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    contract = json.loads(CONTRACT.read_text())
    h7 = json.loads(H7_MANIFEST.read_text())
    if h7["status"] != contract["source"]["required_h7_status"]:
        raise RuntimeError("H7 metadata status mismatch")
    gate = h7["premeasurement_gate"]
    if gate["decision"] != contract["source"]["required_h7_decision"] or gate["pass"] is not False:
        raise RuntimeError("H7 terminal decision mismatch")
    if int(h7["retained_fresh_photos"]) != int(contract["source"]["required_h7_fresh_photos"]):
        raise RuntimeError("H7 fresh-photo count mismatch")
    if bool(h7["candidate_image_pixels_opened"]):
        raise RuntimeError("H7 pixels were opened; metadata-only H8 audit is not authorized")

    tab = pd.read_csv(H7_AUDIT)
    required = {"species", "inat_taxon_id", "cell_id", "target_cell_order", "retained"}
    missing = required - set(tab.columns)
    if missing:
        raise RuntimeError(f"H7 target audit missing columns: {sorted(missing)}")
    if len(tab) != 312 or tab["inat_taxon_id"].nunique() != 52:
        raise RuntimeError("H7 target audit frame mismatch")
    if not (tab.groupby("inat_taxon_id")["cell_id"].nunique() == 6).all():
        raise RuntimeError("H7 target audit must contain six unique cells per species")

    thresholds = [int(x) for x in contract["audit_matrix"]["minimum_retained_fresh_metadata_per_cell"]]
    cell_requirements = [int(x) for x in contract["audit_matrix"]["minimum_qualifying_cells_per_species"]]
    minimum_species = int(contract["audit_matrix"]["minimum_species_required_for_h8"])

    matrix: list[dict[str, int | bool]] = []
    for threshold in thresholds:
        per_species = (
            tab.assign(ok=tab["retained"].astype(int) >= threshold)
            .groupby(["species", "inat_taxon_id"], observed=True)["ok"]
            .sum()
            .reset_index(name="qualifying_cells")
        )
        for cells in cell_requirements:
            n = int((per_species["qualifying_cells"] >= cells).sum())
            matrix.append({
                "minimum_retained_fresh_metadata_per_cell": threshold,
                "minimum_qualifying_cells_per_species": cells,
                "qualifying_species": n,
                "meets_minimum_species": bool(n >= minimum_species),
            })

    eligible = [row for row in matrix if row["meets_minimum_species"]]
    selected_design = None
    selected_species = pd.DataFrame(columns=["species", "inat_taxon_id", "qualifying_cells"])
    selected_targets = pd.DataFrame(
        columns=["species", "inat_taxon_id", "h8_target_order", "cell_id", "h7_retained_fresh_metadata"]
    )
    if eligible:
        selected_design = sorted(
            eligible,
            key=lambda row: (
                int(row["minimum_retained_fresh_metadata_per_cell"]),
                int(row["minimum_qualifying_cells_per_species"]),
            ),
            reverse=True,
        )[0]
        threshold = int(selected_design["minimum_retained_fresh_metadata_per_cell"])
        cells_required = int(selected_design["minimum_qualifying_cells_per_species"])
        support = (
            tab.assign(ok=tab["retained"].astype(int) >= threshold)
            .groupby(["species", "inat_taxon_id"], observed=True)["ok"]
            .sum()
            .reset_index(name="qualifying_cells")
        )
        selected_species = support.loc[support["qualifying_cells"] >= cells_required].copy()
        selected_species = selected_species.sort_values(
            ["inat_taxon_id", "species"], kind="mergesort"
        ).reset_index(drop=True)
        keep_ids = set(selected_species["inat_taxon_id"].astype(int).tolist())
        targets = tab.loc[tab["inat_taxon_id"].astype(int).isin(keep_ids)].copy()
        targets["retained"] = targets["retained"].astype(int)
        targets["cell_id"] = targets["cell_id"].astype(int)
        targets = targets.sort_values(
            ["inat_taxon_id", "retained", "cell_id"],
            ascending=[True, False, True],
            kind="mergesort",
        )
        targets = targets.groupby("inat_taxon_id", sort=False, observed=True).head(5).copy()
        targets["h8_target_order"] = targets.groupby("inat_taxon_id", sort=False).cumcount() + 1
        selected_targets = targets[
            ["species", "inat_taxon_id", "h8_target_order", "cell_id", "retained"]
        ].rename(columns={"retained": "h7_retained_fresh_metadata"})
        if not (selected_targets.groupby("inat_taxon_id")["cell_id"].nunique() == 5).all():
            raise RuntimeError("H8 target-cell selection failed to return five cells per species")

    OUT_SPECIES.parent.mkdir(parents=True, exist_ok=True)
    selected_species.to_csv(OUT_SPECIES, index=False, lineterminator="\n")
    selected_targets.to_csv(OUT_TARGETS, index=False, lineterminator="\n")

    result = {
        "protocol": contract["protocol"],
        "status": "complete_h8_metadata_only_design_audit",
        "h7_terminal_decision_preserved": gate["decision"],
        "h7_pixels_opened": False,
        "matrix": matrix,
        "automatic_selection": {
            "minimum_species_required": minimum_species,
            "selected_design": selected_design,
            "selected_species": int(len(selected_species)),
            "target_cells_per_species": 5 if selected_design is not None else 0,
            "species_cell_targets": int(len(selected_targets)),
            "h8_metadata_design_feasible": bool(selected_design is not None),
        },
        "outcome_firewall": {
            "h7_colour_columns_read": False,
            "h7_measured_table_read": False,
            "h6_h6b_species_effects_read": False,
            "h8_api_queries_opened": False,
            "h8_pixels_opened": False,
        },
        "lineage": {
            "contract_sha256": sha256_file(CONTRACT),
            "h7_manifest_sha256": sha256_file(H7_MANIFEST),
            "h7_target_audit_sha256": sha256_file(H7_AUDIT),
            "selected_species_sha256": sha256_file(OUT_SPECIES),
            "selected_targets_sha256": sha256_file(OUT_TARGETS),
        },
        "claim_ceiling": contract["claim_ceiling"],
    }
    OUT_JSON.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
