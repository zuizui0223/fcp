"""Publication evidence guards; never load reserve outcomes or rerun an analysis."""
import copy
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("rgfca_figures", ROOT / "scripts/analysis/make_rgfca_publication_figures.py")
figures = importlib.util.module_from_spec(spec)
spec.loader.exec_module(figures)


@pytest.fixture(scope="module")
def inputs():
    return figures.load_sources()


def test_frozen_inputs_and_reconstruction(inputs):
    result, species, null, measured, hashes = inputs
    pool, checks = figures.validate_inputs(result, species, null, measured)
    assert len(pool) == 21424 and checks["species"] == 369
    assert checks["observed_mean_rho"] == pytest.approx(.02702130374565584, abs=1e-14)
    assert checks["p_upper"] == .001
    assert set(hashes) == set(figures.SOURCES.values())
    assert all("reserve" not in p and "jbi" not in p for p in hashes)


@pytest.mark.parametrize("target", ["species", "null", "measured"])
def test_partial_or_duplicate_inputs_rejected(inputs, target):
    result, species, null, measured, _ = inputs
    frames = dict(species=species, null=null, measured=measured)
    frames[target] = frames[target].iloc[:-1].copy()
    with pytest.raises(ValueError, match="census"):
        figures.validate_inputs(result, **frames)


def test_no_significant_species_selection(inputs):
    result, species, null, measured, _ = inputs
    selected = species.loc[species.bh_detectable_0_05]
    assert len(selected) == 9
    with pytest.raises(ValueError, match="Species census"):
        figures.validate_inputs(result, selected, null, measured)


def test_invalid_colour_is_not_silently_repaired(inputs):
    result, species, null, measured, _ = inputs
    measured = measured.copy()
    idx = measured.loc[measured.global_classifiable & measured.species.isin(species.species)].index[0]
    measured.loc[idx, "colour_white"] = -1
    with pytest.raises(ValueError, match="colour vectors"):
        figures.validate_inputs(result, species, null, measured)


def test_frozen_result_mismatch_is_not_relabelled(inputs):
    result, species, null, measured, _ = inputs
    result = copy.deepcopy(result)
    result["primary"]["p_upper"] = .05
    with pytest.raises(ValueError, match="Recomputed p_upper"):
        figures.validate_inputs(result, species, null, measured)


def test_map_retains_all_points_without_species_colour_labels(inputs):
    result, species, null, measured, _ = inputs
    pool, _ = figures.validate_inputs(result, species, null, measured)
    rgb = figures.display_rgb(pool)
    assert rgb.shape == (21424, 3)
    assert np.isfinite(rgb).all() and (rgb >= 0).all() and (rgb <= 1).all()


def test_two_renders_identical_and_manifest_is_guarded(tmp_path):
    first = figures.build(tmp_path / "a", tmp_path / "a.json")
    second = figures.build(tmp_path / "b", tmp_path / "b.json")
    assert first == second
    assert len(first["outputs"]) == 4
    for row in first["outputs"]:
        assert row["bytes"] > 1000
    for key in ["reserve_outcomes_read", "legacy_results_read", "new_inferential_tests_run",
                "photo_bar_present", "independent_replication_claim_allowed", "shared_boundary_claim_allowed"]:
        assert first[key] is False


def test_manuscript_explicitly_incomplete():
    text = (ROOT / "docs/RGFCA_MANUSCRIPT.md").read_text(encoding="utf-8")
    for required in ["not submission-ready", "0.0270213", "0.001", "independent replication",
                     "not a confidence interval", "not evidence of absence", "photo bar"]:
        assert required in text
    assert "rgfca_figure1_discovery_atlas.png" in text
    assert "rgfca_figure2_discovery_omnibus.png" in text


def test_committed_figure_bytes_and_source_ledger(inputs):
    manifest = json.loads((ROOT / "docs/supporting/rgfca_publication_figure_manifest_v1.json").read_text(encoding="utf-8"))
    assert manifest["input_sha256"] == inputs[-1]
    assert manifest["source_commit"] == figures.SOURCE_COMMIT
    assert len(manifest["outputs"]) == 4
    for row in manifest["outputs"]:
        path = ROOT / "docs/figures" / row["name"]
        assert path.stat().st_size == row["bytes"]
        assert figures.digest(path.read_bytes()) == row["sha256"]
