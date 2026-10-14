"""Source-receipt checks for global photographed colour diversity; no reclassification."""
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/np_unified_20261014/source_pr135_global_photographic_colour_geography.json"
MANUSCRIPT = ROOT / "docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST_UNIFIED_20261014.md"
SUPPORTING = ROOT / "docs/POLYMORPHISM_SUPPORTING_INFORMATION_UNIFIED_20261014.md"


def data():
    return json.loads(SOURCE.read_text(encoding="utf-8"))


def test_original_global_breadth_four_hue_totals_and_missingness():
    q = data()
    b = q["breadth_single_original_photo_per_species"]
    hues = b["colour_counts"]
    assert q["original_species_denominator"] == 42111
    assert b["species_source_rows"] == 42111
    assert b["species_four_colour_classifiable"] == 18457
    assert sum(hues.values()) == 18457
    assert list(hues.values()) == [8597, 6323, 2001, 1536]
    assert b["unclassifiable"] == 42111 - 18457 == 23654
    for key, value in hues.items():
        assert math.isclose(
            b["conditional_colour_fractions"][key],
            value / b["species_four_colour_classifiable"],
            rel_tol=0, abs_tol=1e-14
        )


def test_original_taxon_cell_and_crosscell_pairs_not_biological_genetics():
    q = data()
    r = q["regional_taxon_cell_measurement"]
    assert r["all_taxon_cell_rows"] == 85337
    assert r["classified_taxon_cell_rows"] == 39075
    p = q["observer_disjoint_crosscell_morph_comparison"]
    assert p["n_fixed_pairs"] == 13416
    assert p["n_identified_pairs"] == 3196
    assert p["n_discordant_pairs"] == 797
    assert p["n_fixed_pairs"] - p["n_identified_pairs"] == 10220
    world = next(x for x in p["regional_and_global_descriptive"] if x["region"] == "ALL")
    assert world["same_photo_colour"] + world["discordant_photo_colour"] == world["both_photos_classified"]
    assert math.isclose(world["discordance_fraction_among_both_classified"],797/3196,abs_tol=1e-13)
    assert world["species_bootstrap_95CI_conditional"][0] < 797/3196 < world["species_bootstrap_95CI_conditional"][1]
    assert world["no_assumption_full_fixed_pair_discordance_bounds"][0] < 797/3196 < world["no_assumption_full_fixed_pair_discordance_bounds"][1]
    assert any("not genetic" in x or "not genetic or biochemical" in x for x in q["hard_nonclaims"])


def test_latitude_pair_counts_have_unassignable_midpoints():
    p = data()["observer_disjoint_crosscell_morph_comparison"]
    sectors = [x for x in p["regional_and_global_descriptive"] if x["region"] != "ALL"]
    assert len(sectors) == 3
    assert sum(x["fixed_observer_disjoint_species_pairs"] for x in sectors) + p["n_not_assignable_pair_midpoints"] == p["n_fixed_pairs"]
    assert sum(x["both_photos_classified"] for x in sectors) == 1168+1646+374


def test_manuscript_and_s19_retain_correct_scope_and_limit():
    manuscript=MANUSCRIPT.read_text(encoding="utf-8")
    s=SUPPORTING.read_text(encoding="utf-8")
    assert "### Global photographed flower-colour diversity and cross-cell variation" in manuscript
    assert manuscript.index("### Global photographed flower-colour diversity and cross-cell variation") < manuscript.index("### Ecology: worldwide flower-colour correspondence")
    for token in ("42,111", "18,457", "85,337", "39,075", "13,416", "3,196", "797", "10,220"):
        assert token in manuscript,token
    assert "### S19.1." in s and "### S19.2." in s
    assert "global polymorphic-species prevalence" in manuscript
    assert "genetic polymorphism incidence" in s
    assert "globally unbiased estimate" in s and "not" in s.lower()
