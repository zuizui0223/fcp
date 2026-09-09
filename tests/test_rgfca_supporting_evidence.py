"""Read-only publication reconstruction of completed discovery evidence.

All inputs are literal immutable Git objects. Never glob result directories,
open reserve outcomes, fetch images or rerun an inferential permutation.
"""
import csv
from functools import lru_cache
import hashlib
import io
import json
from pathlib import Path
import statistics
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]
DISCOVERY = "730612e4af6c64178b1c079d23bb6dc802e1dac0"
HETEROGENEITY = "f2f9c58e1d857a5d4b5b35a78adab8e8324b25eb"
SOURCES = {
    "measurement": ("global_monte_carlo_measurement_result_v1.json", "864a13f22228e7feff65bbac2fbbd276cc54ec35"),
    "omnibus": ("global_rgfca_within_species_spatial_omnibus_result_v1.json", "a0f0f87d41bdc57365b52c664072cf610fb3220f"),
    "g1": ("global_rgfca_g1_result_v1.json", "eb2acdf8f642cd13197426d2aed58ef902621f05"),
    "robustness": ("global_rgfca_prespecified_robustness_result_v1.json", "881fe4abaf91200fa85ca851d77b1d865bfdcbe9"),
    "commonness": ("global_rgfca_species_disjoint_commonness_result_v1.json", "71ecaf582a78d5830210574c6eafc808e7bcf391"),
    "environment": ("global_rgfca_expanded_environmental_panel_final_result_v1.json", "88f9771fd49ccfdd0f82cfa5cc47f5216d7e7d6a"),
    "interactions": ("global_rgfca_environmental_interaction_family_verification_v1.json", "c2ad1e29dacf63e7a2ef2eb9994f3281d43292d3"),
    "climate_qualification": ("hypervolume_real_climate_synthetic_qualification_result_v1.json", "166859f75254d7a1809d419e9369f5ef070ba682"),
    "sharedness_qualification": ("global_rgfca_sharedness_specific_predictive_result_v1.json", "1c92d11f82101634b8903ec2a2a17d1f5546f82a"),
}
HET_SOURCES = {
    "result": ("docs/supporting/global_rgfca_environmental_species_heterogeneity_inference_result_v1.json", "c9b8e29c6a9b8f8588417031c9394309be973a1c17aca97c677ce89f91defd36"),
    "null": ("data/derived/global_rgfca_environmental_species_heterogeneity_inference_null_v1.csv", "810bb23a89f651fa6c184ccaa60ee5352f27f6429db134ec729c4448413600fd"),
    "species": ("data/derived/global_rgfca_environmental_species_heterogeneity_inference_species_v1.csv", "ecb95ed7e92145df44957cf8877b602e9d648278bf5042d59e158083f620b138"),
}


@lru_cache(None)
def source(key):
    name, expected = SOURCES[key]
    raw = subprocess.check_output(["git", "show", f"{DISCOVERY}:docs/supporting/{name}"], cwd=ROOT)
    oid = hashlib.sha1(f"blob {len(raw)}\0".encode() + raw).hexdigest()
    assert oid == expected, f"Committed source differs: {name}"
    return json.loads(raw)


@lru_cache(None)
def heterogeneity(key):
    path, expected = HET_SOURCES[key]
    raw = subprocess.check_output(["git", "show", f"{HETEROGENEITY}:{path}"], cwd=ROOT)
    assert hashlib.sha256(raw).hexdigest() == expected
    return json.loads(raw) if key == "result" else list(csv.DictReader(io.StringIO(raw.decode("utf-8"))))


def holm(values):
    adjusted = [None] * len(values)
    running = 0.0
    for rank, index in enumerate(sorted(range(len(values)), key=values.__getitem__)):
        running = max(running, min(1.0, (len(values) - rank) * values[index]))
        adjusted[index] = running
    return adjusted


def test_all_literal_discovery_sources_and_population_counts():
    for key in SOURCES:
        assert source(key)["status"].startswith("complete")
    measured, omnibus = source("measurement"), source("omnibus")
    assert measured["terminal_result_rows"] == measured["classifiable_rows"] + measured["mixed_uncertain_rows"] == 50000
    assert measured["classifiable_rows"] == 25377
    assert measured["postmeasurement_gate"]["evaluable_species"] == omnibus["eligible_species"] == 369
    assert omnibus["classifiable_rows"] == 21424
    assert omnibus["primary"]["observed_mean_rho"] == pytest.approx(.02702130374565584)
    assert not omnibus["six_species_used"] and not omnibus["thirty_four_species_used"]


def test_non_support_and_all_spatial_sensitivities_retained():
    g1, robust, common = source("g1"), source("robustness"), source("commonness")
    assert g1["p_upper"] == (1 + g1["null_exceed_or_equal_count"]) / (1 + g1["null_permutations"]) == .07
    assert g1["g1_supported"] is False and robust["primary_g1"]["may_be_rescued_by_this_suite"] is False
    for config, expected in {"coarse": .043, "fine": .006, "support10": .26}.items():
        assert robust["spatial_support"][config]["descriptive_p_upper_not_primary_gate"] == expected
    assert robust["full_strong_stability_claim_not_yet_evaluable"] is True
    values = [r["observed_train_vs_heldout_weighted_r"] for r in common["fold_results"]]
    assert statistics.median(values) == common["observed_global_median_fold_r"]
    assert sum(v > 0 for v in values) == common["observed_positive_fold_count"] == 2
    assert common["p_upper"] == .856 and not common["common_across_species_supported"]


def test_full_environment_and_interaction_corrections():
    env = source("environment")
    rows = env["blocks"]
    assert len(rows) == env["family_size"] == 5
    assert holm([r["raw_p_upper"] for r in rows]) == pytest.approx([r["p_holm"] for r in rows])
    assert not any(r["supported_after_five_block_holm"] for r in rows)
    thermal = next(r for r in rows if r["block"] == "thermal_regime")
    assert thermal["p_holm"] == env["alpha"] == .05  # Strict threshold: not support.
    soil = next(r for r in rows if r["block"] == "edaphic_regime")
    assert soil["block_edge_occurrence_coverage_fraction"] == pytest.approx(.5845144555634617)
    interactions = source("interactions")["pairs"]
    assert len(interactions) == 10
    assert holm([r["raw_p_two_sided"] for r in interactions]) == pytest.approx([r["p_holm"] for r in interactions])
    assert not any(r["supported"] for r in interactions)


@pytest.mark.parametrize("family,prefix,observed_key", [
    ("main_effect_variance_family", "variance__main__", "observed_variance"),
    ("interaction_variance_family", "variance__interaction__", "observed_variance"),
    ("coordinated_syndrome_omnibus", "pc1_fraction__", "observed_pc1_fraction"),
])
def test_parent_branch_heterogeneity_all_tail_values_and_families(family, prefix, observed_key):
    result, null = heterogeneity("result"), heterogeneity("null")
    assert len(null) == result["null_permutations"] == 999
    assert sorted(int(r["null_index"]) for r in null) == list(range(999))
    rows = result[family]
    assert len(rows) == {"main_effect_variance_family": 5, "interaction_variance_family": 10, "coordinated_syndrome_omnibus": 2}[family]
    raw_p = []
    for row in rows:
        suffix = row["feature"] if "feature" in row else row["statistic"].removesuffix("_pc1_fraction")
        observed = row[observed_key]
        values = [float(r[prefix + suffix]) for r in null]
        p = (1 + sum(v >= observed for v in values)) / 1000
        assert p == row["raw_p_upper"]
        assert statistics.mean(values) == pytest.approx(row.get("null_mean_variance", row.get("null_mean")), abs=1e-14)
        raw_p.append(p)
    assert holm(raw_p) == pytest.approx([r["p_holm"] for r in rows])
    assert all(r["p_holm"] >= .05 for r in rows)


def test_parent_branch_species_variances_and_no_rescue():
    result, species = heterogeneity("result"), heterogeneity("species")
    assert len(species) == len({r["species"] for r in species}) == result["n_species"] == 369
    for family, prefix in (("main_effect_variance_family", "main__"), ("interaction_variance_family", "interaction__")):
        for row in result[family]:
            values = [float(r[prefix + row["feature"]]) for r in species]
            assert statistics.variance(values) == pytest.approx(row["observed_variance"], abs=1e-14)
    for key in ("supported_main_variance_features", "supported_interaction_variance_features", "supported_syndrome_statistics"):
        assert result[key] == []
    assert result["parent_mean_effects_reclassified"] is False


def test_failed_synthetic_qualifications_do_not_become_ecological_absence():
    climate = source("climate_qualification")["qualification_gates"]
    assert climate["moderate_full_sharing_recovery_by_block"] == {"thermal_regime": 0.0, "water_balance": .032, "atmospheric_energy_dryness": .068}
    assert climate["synthetic_qualification_pass"] is False
    shared = source("sharedness_qualification")["qualification_gates"]
    assert shared["amp1_shared_fraction_1_recovery"] == .008
    assert shared["synthetic_qualification_pass"] is False


def test_supporting_index_retains_source_locations_and_claim_ceiling():
    text = (ROOT / "docs/RGFCA_SUPPORTING_EVIDENCE.md").read_text(encoding="utf-8")
    for name, _ in SOURCES.values():
        assert name in text
    for path, sha in HET_SOURCES.values():
        assert f"/blob/{HETEROGENEITY}/{path}" in text and sha in text
    for phrase in ("not submission-ready", "not an effect", "not estimate prevalence",
                   "not a biological negative", "58.45%", "No reserve outcome"):
        assert phrase in text
