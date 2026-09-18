import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
FIXTURE = ROOT / "packages" / "disttrait" / "fixtures" / "sf_street_trees_v1.csv"
META = ROOT / "packages" / "disttrait" / "fixtures" / "sf_street_trees_v1_metadata.json"


def test_sf_street_tree_fixture_integrity() -> None:
    frame = pd.read_csv(FIXTURE)
    meta = json.loads(META.read_text(encoding="utf-8"))

    assert len(frame) == 1600
    assert frame["taxon_label"].nunique() == 20
    assert set(frame.groupby("taxon_label").size().tolist()) == {80}
    assert frame["tree_id"].nunique() == len(frame)
    assert (frame["dbh"] > 0).all()
    assert np.isfinite(frame[["dbh", "latitude", "longitude"]].to_numpy(float)).all()

    assert meta["source"]["git_blob_sha"] == "bdc06c1297b7dd88bea0df77de4007eaab30198e"
    assert meta["selection"]["eligible_taxa"] == 52
    assert meta["selection"]["selected_taxa"] == 20
    assert meta["selection"]["selected_individuals_per_taxon"] == 80
    assert meta["selection"]["total_rows"] == 1600
    assert "PDDL" in meta["source"]["upstream_license"]
