"""Fail-closed checks for the new *photographic*, not genetic, 447500-photo audit."""
import json
from pathlib import Path
import math

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "results/fcp_photo_label_error_audit_20261023/aggregate_results.json"
LEDGER = ROOT / "results/fcp_np_all4475_unified_species_ledger_receipt_20261021.json"
MANUSCRIPT = ROOT / "docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md"
NOTE = ROOT / "docs/FCP_PHOTOGRAPHIC_MULTICOLOUR_ERROR_AUDIT_20261023.md"

def test_original_photo_source_and_arithmetic_richness():
    a = json.loads(AUDIT.read_text())
    l = json.loads(LEDGER.read_text())
    c = a["species_counts"]
    assert a["source_photo_sha256"] == l["source_actual_photo_sha256"]
    assert a["source_ledger_sha256"] == l["merged_species_ledger_sha256"]
    assert a["execution_receipt"]["source_run_id"] == 39879387172
    assert a["execution_receipt"]["raw_photo_identifiers_or_exact_coordinates_exported"] is False
    assert c["n_species"] == 4475
    assert c["n_original_photo_rows"] == 447500
    assert c["n_four_colour_classifiable_photo_rows"] == 227107
    assert c["n_unclassified_photo_rows"] == 220393
    assert c["observed_n_colour_classes"] == {
        "1":286,"2":1743,"3":1632,"4":814}
    tier=c["colour_photo_count_tiers"]
    assert tier["8"]["at_least_2_colours_each_ge_5"] == 2400
    assert tier["8"]["at_least_3_colours_each_ge_5"] == 435
    assert tier["8"]["at_least_4_colours_each_ge_5"] == 19
    assert tier["40"]["at_least_2_colours_each_ge_5"] == 1592
    assert tier["40"]["at_least_3_colours_each_ge_5"] == 309
    assert tier["40"]["at_least_4_colours_each_ge_5"] == 11
    assert tier["8"]["at_least_4_colours_each_ge_10"] == 0
    assert c["minimum_misclassified_photos_required_to_be_monochrome"]["5"] == 2470

def test_original_observer_and_spatial_null_are_not_conflated_with_morph_genetics():
    a=json.loads(AUDIT.read_text())
    v=a["original_photo_and_observer_check"]
    assert v["max_original_photos_per_observer_per_taxon"] == 2
    assert v["n_species_two_colours_each_ge5_photos_and_3_observers_and_3_observations"] == 2400
    assert v["n_species_three_colours_each_ge5_photos_and_3_observers_and_3_observations"] == 435
    assert v["n_species_four_colours_each_ge5_photos_and_3_observers_and_3_observations"] == 19
    # There are no *comparable* repeated photos from the same observation.
    assert v["same_observation_multi_photo_comparisons"] == 0
    assert v["same_observation_conflicting_colour_label_groups"] == 0
    assert "uninformative" in a["execution_receipt"]["within_observation_diagnostic_note"]
    t=v["species_conditioned_spatial_null_sensitivity"]
    checks={
      "full_original":(3600,0.016115883798023662),
      "one_per_observer":(3592,0.01590008209481584),
      "photo_5_each_n40_pairs30":(813,0.028549525382587742),
      "one_per_observer_photo_5_each_n40_pairs30":(727,0.026544914587817294),
      "photo_10_each_n40_pairs30":(462,0.03743967347859275),
      "one_per_observer_photo_10_each_n40_pairs30":(406,0.036319269463526455),
      "photo_5_each_3types_n40_pairs30":(139,0.031084283455232835)
    }
    for name,(n,mean) in checks.items():
        o=t[name]
        assert o["n_species"] == n
        assert math.isclose(o["mean_delta"],mean,rel_tol=0,abs_tol=1e-12)
        assert o["upper_p_199"] == 0.005
        assert o["permutations"] == 199
        assert o["mean_delta"] > o["mean_null"]
    doc=NOTE.read_text()
    for token in ("0.015900","0.026545","0.036319","0.300","0.275"):
        assert token in doc
    assert "mathematically redundant" in doc
    assert "not established" in doc.lower() or "not establish" in doc.lower()
    paper=MANUSCRIPT.read_text()
    assert "not a prevalence survey" in paper
    assert "No native-range restriction or explicit captive/wild filter" in paper
    assert "genetic" in paper.lower()

def test_new_phytologist_manuscript_preserves_photo_error_caveats():
    paper=MANUSCRIPT.read_text(encoding="utf-8")
    methods=paper.split("## Results\n",1)[1].split("## Discussion\n",1)[0]
    discussion=paper.split("## Discussion\n",1)[1].split("## Acknowledgements\n",1)[0]
    for token in ("2,400 species","435","19","3,592","+0.015900","813",
                  "+0.028550","727","+0.026545","FCP_PHOTOGRAPHIC_MULTICOLOUR_ERROR_AUDIT_20261023.md"):
        assert token in methods, token
    for token in ("No same-observation repeat-photo", "remains unidentified"):
        assert token in discussion, token
    assert "within-population polymorphism" in paper
    assert "not a prevalence survey" in paper
