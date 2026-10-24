"""Recheck source-derived colour-label richness; do not infer genetic morph clustering."""
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "scripts/analysis/fcp_observed_colour_cluster_qa_20261023.py"
spec = importlib.util.spec_from_file_location("fcp_colour_cluster_source_qa", MODULE)
qa = importlib.util.module_from_spec(spec)
spec.loader.exec_module(qa)


def test_original_source_photo_denominators_and_exact_richness():
    rows = qa.rows_from_csv(qa.SOURCE)
    assert len(rows) == 4475
    assert sum(r["n"] for r in rows) == 227107
    assert sum(r["n"]+r["missing"] for r in rows) == 447500
    from collections import Counter
    assert dict(Counter(r["k"] for r in rows)) == {
        1: 286, 2: 1743, 3: 1632, 4: 814
    }
    assert sum(r["k"] >= 2 for r in rows) == 4189
    assert len({r["taxon_id"] for r in rows}) == 4475


def test_minority_colour_support_collapses_apparent_three_four_label_cases():
    r = qa.analyse(qa.rows_from_csv(qa.SOURCE))
    supported = r["each_detected_colour_ge5_and_atleast40photos"]
    assert [supported[str(i)]["species"] for i in range(1,5)] == [
        268, 488, 134, 11
    ]
    t = r["nested_photo_support_sensitivity"]
    assert t["at_least_two_photo_classes_ge5_and_nclass_ge40"]["species"] == 1592
    assert t["at_least_three_photo_classes_ge5_and_nclass_ge40"]["species"] == 309
    assert t["four_photo_classes_ge10_and_nclass_ge40"]["species"] == 0
    assert r["minor_colour_threshold_is_exploratory_not_biological_validation"]
    assert r["photo_colour_classes_are_not_unsupervised_9_palette_modes"]
    assert r["all_independent_observers_per_colour_unknown_in_aggregate"]


def test_original_spatial_gate_is_untouched_and_subgroup_null_not_claimed():
    r = qa.analyse(qa.rows_from_csv(qa.SOURCE))
    n = sum(r["observed_colour_label_richness_1to4"][str(k)][
        "species_strict_50km_ge40_ge30pairs"] for k in range(1,5))
    assert n == 1761
    assert r["observed_colour_label_richness_1to4"]["1"][
        "strict_50km_mean_global_minus_local_mismatch"] == 0
    assert r["observed_colour_label_richness_1to4"]["4"][
        "species_strict_50km_ge40_ge30pairs"] == 251
    assert r["each_detected_colour_ge5_and_atleast40photos"]["4"][
        "species_strict_50km_ge40_ge30pairs"] == 3
    assert r["only_observed_geographic_delta_no_new_stratum_conditioned_spatial_null"]
    assert r["original_matched_source_50km_null_not_retested_by_this_summary"]
    assert r["matched_genotype_or_petal_stage_fitness_evidence"] is False
    assert r["strict_5km_source_null_p_ge5pairs"] == 0.3


def test_synthetic_cluster_support_edge_cases():
    qa.self_test()

def test_original_photo_mode_quality_geo_result_is_real_and_not_genetic():
    import json
    path = ROOT / "results/fcp_original_photo_palette_mode_geo_quality_audit_20261023.json"
    d = json.loads(path.read_text(encoding="utf-8"))
    assert d["authentic_photo_sha256"] == qa.PHOTO_SHA
    assert d["source_photo_rows"] == 447500
    assert d["classifiable_photo_rows"] == 227107
    assert d["verified_ci_run_id"] == 39880111954
    assert d["verified_ci_conclusion"] == "success"
    assert d["full_anonymized_json_artifact_id"] == 12469982240
    assert d["full_artifact_zip_sha256"] == "f337b7dfe2b44a7dd8e943922593734acad8db32d64382af22a415c174e39d7a"
    q = d["photo_quality_observed_class_support"]
    assert (q["q65_ge2_supported_colour_groups"],
            q["q65_ge3_supported_colour_groups"],
            q["q65_ge4_supported_colour_groups"]) == (895, 81, 0)
    assert (q["q75_ge2_supported_colour_groups"],
            q["q75_ge3_supported_colour_groups"],
            q["q75_ge4_supported_colour_groups"]) == (445, 26, 0)
    results = d["quality_conditioned_50km_species_matched_label_null"]
    g = results["all_supported_two_or_more"]
    assert g["species"] == 306
    assert abs(g["mean_source_photo_global_minus_local"]-0.04614868076348339) < 1e-11
    assert abs(g["distinct_observer_pair_mean_depletion"]-0.04557112043546877) < 1e-11
    assert g["matched_199_world_upper_p"] == 0.005
    assert g["distinct_observer_pair_matched_199_world_upper_p"] == 0.005
    assert results["exactly_2_supported_classes"]["species"] == 279
    assert results["exactly_3_supported_classes"]["species"] == 27
    assert results["exactly_4_supported_classes"]["species"] == 0
    pilot = d["photographic_PCA2_GMM_BIC_pilot"]
    assert pilot["n_species_tested"] == 160
    assert pilot["n_species_with_photo_mixture_candidate_and_observer_stability"] == 95
    assert pilot["bics_selected_components_1to4"] == {"2": 53, "3": 87, "4": 20}
    assert d["original_P2_frozen_outcome_not_modified"]
    assert d["no_individual_coordinates_observer_or_photo_ids_exported"]
    assert "not a census" in d["note"].lower() or "not a census" in d["nulls"].lower() or "genetically inherited morphs" in d["note"].lower()
    # Past same-year-month 5km null remains negative, not overridden by 50km P.
    assert d["strict_5km_month_matched_source_result_unaltered"]["p"] == 0.3
