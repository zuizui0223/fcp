from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "acquisition" / "run_fcp_v2_fresh_metadata_terminalization_20260923.py"


def load_module():
    spec = importlib.util.spec_from_file_location("fcp_v2_metadata", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def queue_frame(panel: str) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "queue_rank": range(1, 301),
            "panel": [panel] * 300,
            "inat_taxon_id": range(1000 if panel == "P" else 2000, 1300 if panel == "P" else 2300),
            "species": [f"{panel}species {i}" for i in range(1, 301)],
            "after_observer_cap": [150] * 300,
            "selection_hash": [f"{i:064x}" for i in range(1, 301)],
        }
    )


def test_validate_queue_requires_exact_ranked_300() -> None:
    m = load_module()
    q = m.validate_queue(queue_frame("P"), "P")
    assert len(q) == 300
    assert q["queue_rank"].tolist() == list(range(1, 301))
    bad = queue_frame("P")
    bad.loc[299, "queue_rank"] = 299
    with pytest.raises(RuntimeError, match="rank"):
        m.validate_queue(bad, "P")


def test_terminalize_panel_takes_first_200_full_species_by_frozen_rank() -> None:
    m = load_module()
    q = m.validate_queue(queue_frame("P"), "P")
    audit = pd.DataFrame(
        {
            "inat_taxon_id": q["inat_taxon_id"],
            "species": q["species"],
            "retained": [100 if i not in {2, 4, 6} and i <= 213 else 80 for i in range(1, 301)],
            "full_fixed_n": [i not in {2, 4, 6} and i <= 213 for i in range(1, 301)],
            "request_error": [""] * 300,
        }
    )
    full_taxa = audit.loc[audit["full_fixed_n"], "inat_taxon_id"].tolist()
    rows = []
    for taxon in full_taxa:
        for j in range(100):
            rows.append(
                {
                    "inat_taxon_id": taxon,
                    "observation_id": taxon * 1000 + j,
                    "photo_id": taxon * 1000 + j + 10_000_000,
                    "h9_selection_order": j + 1,
                }
            )
    observations = pd.DataFrame(rows)
    _, terminal, terminal_rows = m.terminalize_panel(
        panel="P",
        queue=q,
        species_audit=audit,
        observations=observations,
    )
    expected = audit.loc[audit["full_fixed_n"]].merge(
        q[["queue_rank", "inat_taxon_id"]],
        on="inat_taxon_id",
        how="left",
    ).sort_values("queue_rank").head(200)["inat_taxon_id"].tolist()
    assert terminal["inat_taxon_id"].tolist() == expected
    assert terminal["terminal_panel_rank"].tolist() == list(range(1, 201))
    assert len(terminal_rows) == 20_000
    assert terminal_rows["inat_taxon_id"].nunique() == 200


def test_terminalize_panel_does_not_extend_when_underidentified() -> None:
    m = load_module()
    q = m.validate_queue(queue_frame("N"), "N")
    audit = pd.DataFrame(
        {
            "inat_taxon_id": q["inat_taxon_id"],
            "species": q["species"],
            "retained": [100 if i <= 199 else 99 for i in range(1, 301)],
            "full_fixed_n": [i <= 199 for i in range(1, 301)],
            "request_error": [""] * 300,
        }
    )
    observations = pd.DataFrame(
        columns=["inat_taxon_id", "observation_id", "photo_id", "h9_selection_order"]
    )
    _, terminal, terminal_rows = m.terminalize_panel(
        panel="N",
        queue=q,
        species_audit=audit,
        observations=observations,
    )
    assert len(terminal) == 199
    assert terminal_rows.empty
