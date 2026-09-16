import csv
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "analysis" / "build_polymorphism_h2_third_cohort_frame_20260916.py"


def load_module():
    spec = importlib.util.spec_from_file_location("third_frame", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_csv(path: Path, rows, fieldnames):
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def test_p500_selected_species_are_removed_even_if_not_measured():
    m = load_module()
    p100 = [
        {"inat_taxon_id": 1, "species": "A a", "after_observer_cap": 100},
        {"inat_taxon_id": 2, "species": "B b", "after_observer_cap": 120},
        {"inat_taxon_id": 3, "species": "C c", "after_observer_cap": 140},
    ]
    p500 = [{"inat_taxon_id": 2, "species": "B b", "after_observer_cap": 120}]
    out = m.derive_candidate_rows(p100, p500)
    assert [r["inat_taxon_id"] for r in out] == [1, 3]


def test_candidate_rows_are_sorted_and_output_only_firewall_columns():
    m = load_module()
    p100 = [
        {"inat_taxon_id": 3, "species": "C c", "after_observer_cap": 140},
        {"inat_taxon_id": 1, "species": "A a", "after_observer_cap": 100},
    ]
    out = m.derive_candidate_rows(p100, [])
    assert [r["inat_taxon_id"] for r in out] == [1, 3]
    assert list(out[0]) == ["inat_taxon_id", "species", "after_observer_cap"]


def test_source_reader_rejects_colour_outcome_columns(tmp_path):
    m = load_module()
    p = tmp_path / "bad.csv"
    write_csv(
        p,
        [{"inat_taxon_id": 1, "species": "A a", "after_observer_cap": 100, "morph": "white"}],
        ["inat_taxon_id", "species", "after_observer_cap", "morph"],
    )
    with pytest.raises(RuntimeError, match="forbidden or unexpected source columns"):
        m.read_outcome_blind_source(p)


def test_source_reader_requires_u100_capacity(tmp_path):
    m = load_module()
    p = tmp_path / "low.csv"
    write_csv(
        p,
        [{"inat_taxon_id": 1, "species": "A a", "after_observer_cap": 99}],
        ["inat_taxon_id", "species", "after_observer_cap"],
    )
    with pytest.raises(RuntimeError, match="below frozen U100 gate"):
        m.read_outcome_blind_source(p)
