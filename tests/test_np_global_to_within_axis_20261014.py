"""Guard the maximum currently analyzable H2 scale against invented prospective samples."""
from __future__ import annotations
import importlib.util
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts/analysis/reconcile_np_global_to_within_axis_20261014.py"
RECEIPT=ROOT/"results/np_unified_20261014/global_to_within_white_axis_coverage_and_reconciliation.json"
spec=importlib.util.spec_from_file_location("fcp_global_to_within_audit",SCRIPT)
audit=importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def test_all_3_original_cohorts_contribute_to_descriptive_axis_only():
    got=audit.calculate()
    stored=json.loads(RECEIPT.read_text(encoding="utf8"))
    assert got==stored
    primary=got["white_axis_descriptive_all_cohort_summary"]["primary"]
    strict=got["white_axis_descriptive_all_cohort_summary"]["strict"]
    assert primary["eligible_species_total"]==152+129+158==439
    assert strict["eligible_species_total"]==75+65+86==226
    assert math.isclose(primary["species_equal_weighted_observed_W"],0.5155565569772226,abs_tol=1e-12)
    assert math.isclose(primary["species_equal_weighted_cohort_null_MEANS"],0.45036227000622664,abs_tol=1e-12)
    assert primary["no_pooled_null_worlds_constructed"]
    assert primary["no_pooled_p_value_computed"]
    assert primary["do_not_relabel_post_audit_cohorts_as_prospective"]
    rows=got["white_axis_already_evaluable_three_cohorts"]
    assert [x["inferential_role"] for x in rows]==[
        "retrospective target localization",
        "retrospective target localization",
        "untouched prospective confirmation",
    ]
    assert [x["primary"]["source_original_null_p"] for x in rows]==[.001,.001,.001]


def test_global_original_measured_colour_not_genetic_polymorphism_prevalence():
    got=audit.calculate()
    b=got["worldwide_four_colour_coverage"]
    assert (b["one_photo_nominal_species"],b["one_photo_colour_classifiable_species"])==(42111,18457)
    assert sum(b["one_photo_original_four_colour_counts"].values())==18457
    assert b["taxon_geographic_cell_photo_records"]==85337
    assert b["taxon_cell_colour_classifiable"]==39075
    assert (b["merged_unique_photo_ids"],b["merged_photo_four_colour_classifiable"])==(100543,45953)
    assert not b["global_colour_polymorphism_prevalence_estimable"]
    p=got["cross_geographic_cell_photo_pairs"]
    assert p["n_fixed_observer_disjoint_pairs"]==13416
    assert p["both_classifiable_pairs"]==3196
    assert p["discordant_four_colour_pairs"]==797
    assert p["unknown_hue_in_either"]==10220
    assert math.isclose(p["discordance_among_both_classified"],797/3196,abs_tol=1e-15)
    assert p["no_assumption_full_frame_lower"]<p["discordance_among_both_classified"]<p["no_assumption_full_frame_upper"]
    assert not p["distinct_individual_plants_or_genotypes_identified"]


def test_existing_high_depth_limit_and_next2000_unopened():
    got=audit.calculate()
    h=got["high_depth_sample_coverage"]
    assert h["n_historical_high_depth_photographed_species"]==1499
    assert h["n_at_least_40_classifiable_species_all_cohorts"]==369+363+377==1109
    assert h["n_prospective_classifiable_rows"]==25788
    f=got["new_high_depth_species_not_yet_measured"]
    assert f["historical_u100_metadata_opportunity"]==4730
    assert f["never_allocated_capacity_species"]==2730
    assert f["old_selected_taxon_ids_excluded"]==2000
    assert f["preselected_2000_species_pixels_unopened"]==2000
    assert f["reserve_730_species_unopened"]==730
    assert f["max_possible_if_all_new_2000_measured"]==3499
    assert f["new_colour_axis_confirmation_available"] is False


def test_paper_global_to_within_to_tradeoff_has_provenance_fence():
    path=ROOT/"docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST_UNIFIED_20261014.md"
    s=path.read_text(encoding="utf8")
    r=s[s.index("## Results"):s.index("## Discussion")]
    d=s[s.index("## Discussion"):s.index("## Acknowledgements")]
    assert r.index("### Global photographed flower-colour diversity")<r.index("### Pattern: achromatic–chromatic variation")
    assert "439" in r
    assert "158" in r
    assert "not a pooled" in r.lower() or "no pooled" in r.lower()
    assert "relative fitness" in d
    assert "P = 0.214" in r
    si=(ROOT/"docs/POLYMORPHISM_SUPPORTING_INFORMATION_UNIFIED_20261014.md").read_text(encoding="utf8")
    assert "## S20." in si
    assert "439" in si and "226" in si
    assert "2,000" in si and "not newly measured" in si.lower()
