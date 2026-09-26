import hashlib
import importlib.util
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "analysis" / "make_polymorphism_manuscript_figures.py"
SOURCE_DIR = ROOT / "results" / "polymorphism_publication_figure_source_20260918"
SOURCE_MANIFEST = SOURCE_DIR / "manifest.json"
H2_PROSPECTIVE = ROOT / "results" / "polymorphism_h2_third_cohort_prospective_white_axis_20260917" / "result.json"
H2_PRIMARY_NULL = ROOT / "results" / "polymorphism_h2_third_cohort_prospective_white_axis_20260917" / "primary_0_10_structured_null.csv"
H2_STRICT_NULL = ROOT / "results" / "polymorphism_h2_third_cohort_prospective_white_axis_20260917" / "strict_0_20_structured_null.csv"
SPATIAL = ROOT / "results" / "polymorphism_spatial_organization_clue_20260918" / "result.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_module():
    assert SCRIPT.exists(), "publication figure generator is not implemented"
    spec = importlib.util.spec_from_file_location("polymorphism_figures", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_reporting_source_freeze_matches_artifact_derived_hashes() -> None:
    manifest = json.loads(SOURCE_MANIFEST.read_text(encoding="utf-8"))
    assert manifest["schema"] == "polymorphism_publication_figure_source_v1"
    assert manifest["purpose"].startswith("reporting-only")

    for name, expected in manifest["files"].items():
        path = SOURCE_DIR / name
        assert path.exists()
        assert sha256(path) == expected["sha256"]

    d = pd.read_csv(SOURCE_DIR / "h3a_species_D.csv")
    assert len(d) == 732
    assert d.groupby("cohort").size().to_dict() == {"discovery": 369, "reserve": 363}


def test_prospective_h2_null_arrays_are_complete() -> None:
    result = json.loads(H2_PROSPECTIVE.read_text(encoding="utf-8"))
    primary = pd.read_csv(H2_PRIMARY_NULL)
    strict = pd.read_csv(H2_STRICT_NULL)

    assert len(primary) == 999
    assert len(strict) == 999
    assert list(primary.columns) == ["W"]
    assert list(strict.columns) == ["W"]
    assert result["decision"]["verdict"] == "H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED"
    assert result["thresholds"]["primary_0_10"]["species"] == 158
    assert result["thresholds"]["strict_0_20"]["species"] == 86


def test_spatial_reporting_receipt_is_reporting_only_and_frozen() -> None:
    result = json.loads(SPATIAL.read_text(encoding="utf-8"))
    assert result["schema"] == "polymorphism_spatial_organization_reporting_receipt_v1"
    assert result["new_biological_analysis"] is False
    assert result["reserve"]["span_plus_technical_adjusted_primary"]["partial_rho"] == 0.09928771129095708
    assert result["reserve"]["span_plus_technical_adjusted_flower_minus_background"]["partial_rho"] == 0.1162411363016301


def test_generate_all_publication_figures(tmp_path: Path) -> None:
    module = load_module()
    manifest_path = module.generate_all(ROOT, tmp_path)

    expected_stems = [
        "polymorphism_figure1_measurement_frame",
        "polymorphism_figure2_h1_reproducibility",
        "polymorphism_figure3_h2_target_localization",
        "polymorphism_figure4_prospective_h2",
        "polymorphism_figure5_explanatory_boundaries",
    ]
    for stem in expected_stems:
        png = tmp_path / f"{stem}.png"
        pdf = tmp_path / f"{stem}.pdf"
        assert png.exists() and png.stat().st_size > 10_000
        assert pdf.exists() and pdf.stat().st_size > 5_000

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["schema"] == "polymorphism_manuscript_figure_manifest_v1"
    assert manifest["scientific_claims_changed"] is False
    assert manifest["figures"]["figure4"]["primary"]["species"] == 158
    assert manifest["figures"]["figure4"]["primary"]["p"] == 0.001
    assert manifest["figures"]["figure4"]["strict"]["species"] == 86
    assert manifest["figures"]["figure4"]["strict"]["p"] == 0.001
    assert manifest["figures"]["figure5"]["spatial_organization"]["new_biological_analysis"] is False
    assert manifest["figures"]["figure5"]["spatial_organization"]["reserve_adjusted_partial_rho"] == 0.09928771129095708
    assert manifest["figures"]["figure5"]["spatial_organization"]["reserve_adjusted_p"] == 0.025
    assert manifest["figures"]["figure5"]["spatial_organization"]["reserve_background_partial_rho"] == 0.1162411363016301
    assert manifest["figures"]["figure5"]["spatial_organization"]["reserve_background_p"] == 0.01
    assert manifest["figures"]["figure5"]["h3a"]["verdict"] == "H3A_PHYLOGENETIC_SIGNAL_NOT_SUPPORTED"
    assert manifest["figures"]["figure5"]["h3b"]["verdict"] == "H3B_SAMPLED_SPAN_REPLICATION_NOT_SUPPORTED"

    f1_layout = manifest["figures"]["figure1"]["layout_contract"]
    assert f1_layout["cohort_topology"] == "global_frame_branches_to_original_and_third_cohort"
    assert f1_layout["arrow_direction"] == "top_to_bottom"
    assert f1_layout["secondary_followup"] == "dashed_post_h2_reuse_of_third_cohort"

    f2_layout = manifest["figures"]["figure2"]["layout_contract"]
    assert f2_layout["stress_annotations"] == "offset_no_legend_overlap"

    f3_layout = manifest["figures"]["figure3"]["layout_contract"]
    assert f3_layout["legend"] == "outside_below_axis"
