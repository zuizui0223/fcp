"""Synthetic-only guards: do not interpret synthetic outputs as biological results."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from fcp_pipeline.u100_comparison import (
    attach_species_metrics, environment_contrasts, main, opportunity,
    photo_tables, registry,
)


def input_files(tmp_path: Path):
    frame = pd.DataFrame([
        {"inat_taxon_id": 11, "species": "Species alba", "family": "Fam1", "after_observer_cap": 100},
        {"inat_taxon_id": 12, "species": "Species rubra", "family": "Fam2", "after_observer_cap": 100},
        {"inat_taxon_id": 13, "species": "Species ambigua", "family": "Fam3", "after_observer_cap": 100},
        {"inat_taxon_id": 14, "species": "Species missing", "family": "Fam4", "after_observer_cap": 100},
    ])
    frame_path = tmp_path / "U100.csv"
    frame.to_csv(frame_path, index=False)
    rows, env = [], []
    for taxon, name, count in ((11, "Species alba", 50), (12, "Species rubra", 50), (13, "Species ambigua", 10)):
        for i in range(count):
            morph = "white" if taxon != 11 or i < 20 else "red_pink"
            rows.append({"inat_taxon_id": taxon, "species": name, "photo_id": f"p{taxon}-{i}",
                         "morph": morph, "global_classifiable": "true", "observer_id": f"o{i%5}",
                         "latitude": 30 + i * 0.1, "longitude": 130})
            env.append({"inat_taxon_id": taxon, "photo_id": f"p{taxon}-{i}",
                        "bio5": 30 + (5 if morph == "white" else 0) + i/100,
                        "bio14": 40 + i, "srad_mean": 200 + i,
                        "near_clip_fraction": 0.1, "high_clip_excluded": i == 0})
    measured_path = tmp_path / "measured.csv"
    env_path = tmp_path / "env.csv"
    pd.DataFrame(rows).to_csv(measured_path, index=False)
    pd.DataFrame(env).to_csv(env_path, index=False)
    return frame_path, measured_path, env_path


def test_full_entrypoint_preserves_unmeasured_and_monochrome_distinction(tmp_path):
    frame_path, measured_path, env_path = input_files(tmp_path)
    out = tmp_path / "output"
    assert main(["--u100", str(frame_path), "--expected-u100", "4", "--shards", "3",
                 "--measured", f"expansion={measured_path}", "--environment", str(env_path),
                 "--outdir", str(out)]) == 0
    d = pd.read_csv(out / "u100_species_registry.csv").set_index("inat_taxon_id")
    assert d.loc[11, "measurement_status"] == "D_evaluable"
    assert d.loc[11, "D"] == pytest.approx(0.48)
    assert d.loc[11, "D_pair"] == pytest.approx(0.48 * 50 / 49)
    assert d.loc[12, "measurement_status"] == "D_evaluable"
    assert d.loc[12, "D"] == 0  # all white in sampled photographs, NOT a natural-population claim
    assert d.loc[13, "measurement_status"] == "insufficient_classifiable"
    assert pd.isna(d.loc[13, "D"])
    assert d.loc[14, "measurement_status"] == "not_measured"
    assert pd.isna(d.loc[14, "D"])
    contrasts = pd.read_csv(out / "white_environment_species_contrasts.csv")
    bio = contrasts[(contrasts.inat_taxon_id == 11) & (contrasts.variable == "bio5")].iloc[0]
    assert bio.white_minus_nonwhite_sd > 0
    assert bio.paired_observers == 5
    assert np.isfinite(bio.paired_observer_delta)
    result = json.loads((out / "comparison_summary.json").read_text())
    assert result["n_u100"] == 4
    assert result["n_measured"] == 3
    assert result["n_D_evaluable"] == 2
    assert result["confirmatory_decisions_changed"] is False
    assert result["cohorts"]["expansion"]["W_alignment_squared"]["status"] == "not_provided"


def test_shards_are_order_independent(tmp_path):
    frame_path, measured_path, _ = input_files(tmp_path)
    frame = opportunity(frame_path, 4)
    photos, _ = photo_tables({"test": measured_path}, frame)
    a = registry(frame, photos, 7).set_index("inat_taxon_id")
    b = registry(frame.sample(frac=1, random_state=13), photos.sample(frac=1, random_state=9), 7).set_index("inat_taxon_id")
    assert a["shard"].to_dict() == b["shard"].to_dict()
    assert a["D"].fillna(-1).to_dict() == b["D"].fillna(-1).to_dict()


def test_disjoint_cohorts_enforced(tmp_path):
    frame_path, measured_path, _ = input_files(tmp_path)
    frame = opportunity(frame_path, 4)
    with pytest.raises(ValueError, match="multiple disjoint"):
        photo_tables({"a": measured_path, "b": measured_path}, frame)


def test_high_depth_capacity_gate_enforced(tmp_path):
    frame_path, *_ = input_files(tmp_path)
    frame = pd.read_csv(frame_path)
    frame.loc[0, "after_observer_cap"] = 99
    frame.to_csv(frame_path, index=False)
    with pytest.raises(ValueError, match="observer-capped"):
        opportunity(frame_path, 4)


def test_w_requires_frozen_specific_method_and_spatial_D_identity(tmp_path):
    frame_path, measured_path, _ = input_files(tmp_path)
    frame = opportunity(frame_path, 4)
    photos, _ = photo_tables({"extension": measured_path}, frame)
    led = registry(frame, photos, 4)
    path = tmp_path / "w.csv"
    pd.DataFrame([{"inat_taxon_id": 11, "cohort": "extension", "w_alignment_sq": .7,
                   "method": "generic_disttrait_refit"}]).to_csv(path, index=False)
    with pytest.raises(ValueError, match="frozen study-specific"):
        attach_species_metrics(led, path, kind="white_axis")
    spatial = tmp_path / "spatial.csv"
    pd.DataFrame([{"inat_taxon_id": 11, "cohort": "extension", "D_pair": .2,
                   "local_pairs_50km": 40, "depletion_50km": .1}]).to_csv(spatial, index=False)
    with pytest.raises(ValueError, match="differs"):
        attach_species_metrics(led, spatial, kind="spatial")


def test_environment_requires_exact_exposure_controls_and_unique_keys(tmp_path):
    frame_path, measured_path, env_path = input_files(tmp_path)
    frame = opportunity(frame_path, 4)
    photos, _ = photo_tables({"extension": measured_path}, frame)
    env = pd.read_csv(env_path)
    env.drop(columns="high_clip_excluded").to_csv(env_path, index=False)
    with pytest.raises(ValueError, match="missing columns"):
        environment_contrasts(photos, env_path)
    env["high_clip_excluded"] = False
    env = pd.concat([env, env.iloc[[0]]], ignore_index=True)
    env.to_csv(env_path, index=False)
    with pytest.raises(ValueError, match="duplicate"):
        environment_contrasts(photos, env_path)
