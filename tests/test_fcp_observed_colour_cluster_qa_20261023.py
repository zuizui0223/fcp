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
