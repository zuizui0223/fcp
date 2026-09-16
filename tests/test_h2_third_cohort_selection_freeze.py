import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FREEZE = ROOT / "results/polymorphism_h2_third_cohort_selection_20260916"
RESULT = FREEZE / "result.json"
MANIFEST = FREEZE / "selected_species_manifest.tsv"
EXPECTED_CANONICAL_SELECTION_SHA = "4ad1191f39068e0fb2229f84361b1004803e24793190566d21e5c474aef2002a"
EXPECTED_MANIFEST_SHA = "16db6a2fc265f0ab20e7dae4e6f9058b907f5c19af3ddab5b6077797dd1def59"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_selection_receipt_is_durable_and_outcome_blind():
    receipt = json.loads(RESULT.read_text(encoding="utf-8"))
    assert receipt["status"] == "THIRD_COHORT_SELECTION_FROZEN_PREOUTCOME"
    assert receipt["biological_outcomes_opened"] is False
    assert receipt["lineage"]["third_candidate_species"] == 3230
    assert receipt["lineage"]["third_selected_species"] == 500
    assert receipt["selection"]["salt"] == "FCP_H2_THIRD_COHORT_20260916_V1"
    assert receipt["selection"]["canonical_selected_csv_sha256"] == EXPECTED_CANONICAL_SELECTION_SHA
    assert receipt["selection"]["ordered_identity_manifest_sha256"] == EXPECTED_MANIFEST_SHA
    assert sha256(MANIFEST) == EXPECTED_MANIFEST_SHA
    assert receipt["selection"]["p500_taxon_overlap"] == 0
    assert receipt["selection"]["p500_species_overlap"] == 0
    assert receipt["outcome_firewall"]["P500_measurement_table_read"] is False
    assert receipt["outcome_firewall"]["P500_recovered_H2_read"] is False
    assert receipt["outcome_firewall"]["H2_W_opened"] is False


def test_ordered_identity_manifest_contains_exactly_500_unique_ranked_species():
    with MANIFEST.open(encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 500
    assert [int(r["prospective_rank"]) for r in rows] == list(range(1, 501))
    assert len({r["inat_taxon_id"] for r in rows}) == 500
    assert len({r["species"] for r in rows}) == 500
    assert rows[0] == {
        "prospective_rank": "1",
        "inat_taxon_id": "210141",
        "species": "Ajuga australis",
    }
    assert rows[-1] == {
        "prospective_rank": "500",
        "inat_taxon_id": "566867",
        "species": "Oxalis ciliaris",
    }
