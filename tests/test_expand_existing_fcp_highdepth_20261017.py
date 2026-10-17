"""Only synthetic identity/selection guards. Not biological photograph results."""
import csv
import hashlib
from pathlib import Path

import pytest

from scripts.acquisition import expand_existing_fcp_highdepth_20261017 as m


def write_csv(path, rows, names):
    with path.open("w", newline="", encoding="utf8") as f:
        w = csv.DictWriter(f, fieldnames=names, delimiter="\t" if path.suffix == ".tsv" else ",", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def prep(tmp_path, monkeypatch):
    p100 = tmp_path / "p100.csv"
    p500 = tmp_path / "p500.csv"
    old = tmp_path / "third.tsv"
    allrows = [{"inat_taxon_id": i + 10000, "species": f"Species example{i:04d}", "after_observer_cap": 120} for i in range(3730)]
    parents = allrows[:500]
    leftover = allrows[500:]
    oldrows = [{"prospective_rank": rank, "inat_taxon_id": r["inat_taxon_id"], "species": r["species"]} for rank, r in enumerate(leftover[:500], 1)]
    write_csv(p100, allrows, m.SOURCE_COLUMNS)
    write_csv(p500, parents, m.SOURCE_COLUMNS)
    write_csv(old, oldrows, ["prospective_rank", "inat_taxon_id", "species"])
    candidate = "".join([",".join(m.SOURCE_COLUMNS) + "\n"] + [
        f'{r["inat_taxon_id"]},{r["species"]},{r["after_observer_cap"]}\n'
        for r in sorted(leftover, key=lambda x: x["inat_taxon_id"])
    ])
    monkeypatch.setattr(m, "P100_SHA", m.sha256(p100))
    monkeypatch.setattr(m, "P500_SHA", m.sha256(p500))
    monkeypatch.setattr(m, "THIRD_SELECTED_SHA", m.sha256(old))
    monkeypatch.setattr(m, "CANDIDATE_3230_SHA", hashlib.sha256(candidate.encode()).hexdigest())
    return p100, p500, old


def test_repeats_original_500_species_capacity_and_identity_selection(tmp_path, monkeypatch):
    p100, p500, third = prep(tmp_path, monkeypatch)
    out = tmp_path / "out"
    a = m.freeze_selection(p100, p500, third, out)
    b = m.freeze_selection(p100, p500, third, out)
    assert a == b
    assert a["remaining_candidate_species"] == 2730
    assert a["wave4_selected_species"] == 500
    assert a["target_additional_photo_rows"] == 50000
    with (out / "selected_species_manifest.tsv").open() as f:
        data = list(csv.DictReader(f, delimiter="\t"))
    ids = {int(x["inat_taxon_id"]) for x in data}
    assert len(ids) == 500
    assert not (ids & {10000 + i for i in range(1000)})
    assert data[0]["prospective_rank"] == "1"


def test_identity_drift_is_rejected(tmp_path, monkeypatch):
    p100, p500, third = prep(tmp_path, monkeypatch)
    with third.open("a") as f:
        f.write("501\t123\tTestus bogus\n")
    with pytest.raises(RuntimeError, match="SHA drift"):
        m.freeze_selection(p100, p500, third, tmp_path / "out")


def test_parent_file_drift_is_rejected(tmp_path, monkeypatch):
    p100, p500, third = prep(tmp_path, monkeypatch)
    with p500.open("a") as f:
        f.write("99999,Testus bogus,100\n")
    with pytest.raises(RuntimeError, match="SHA drift"):
        m.freeze_selection(p100, p500, third, tmp_path / "out")


def test_cannot_bypass_historical_candidate_hash(tmp_path, monkeypatch):
    p100, p500, third = prep(tmp_path, monkeypatch)
    monkeypatch.setattr(m, "CANDIDATE_3230_SHA", "0"*64)
    with pytest.raises(RuntimeError, match="candidate CSV"):
        m.freeze_selection(p100, p500, third, tmp_path / "out")


def test_exclusion_inputs_required():
    with pytest.raises(RuntimeError, match="five old"):
        m.load_exclusions([])


def test_no_outcome_in_selection_receipt(tmp_path, monkeypatch):
    p100, p500, third = prep(tmp_path, monkeypatch)
    rec = m.freeze_selection(p100, p500, third, tmp_path / "out")
    assert rec["no_image_pixels_opened"] is True
    assert rec["no_flower_colour_used"] is True
    assert "morph" not in rec
    assert rec["overlap_with_old_1000"] == 0
    assert rec["overlap_with_old_500"] == 0
