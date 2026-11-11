"""Editorial checks for one strictly intraspecific New Phytologist FCP paper.

The frozen original H1/H2 results and manuscript are never modified by these tests.
"""
import json
import math
import re
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PAPER=ROOT/"docs/POLYMORPHISM_MANUSCRIPT_NP_INTRASPECIFIC_20261014.md"
ORIGINAL=ROOT/"docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md"
PLAN=ROOT/"docs/POLYMORPHISM_NP_INTRA_PAPER_AND_EVIDENCE_20261014.md"
COVER=ROOT/"docs/POLYMORPHISM_NP_INTRA_COVER_LETTER_20261014.md"
H2=ROOT/"results/polymorphism_h2_third_cohort_prospective_white_axis_20260917/result.json"
H2COHORTS=ROOT/"results/polymorphism_white_axis_three_cohorts_descriptive_20261014.json"
SPATIAL=ROOT/"results/polymorphism_distributed_polymorphism_posthoc_20261007/result.json"


def _section(s,a,b):
    i=s.index(a)+len(a)
    return s[i:s.index(b,i)]


def test_title_and_original_submission_result_unchanged():
    s=PAPER.read_text(encoding="utf-8")
    assert s.startswith("# Recurrent achromatic–chromatic variation and geographic allocation of flower colours within species")
    assert ORIGINAL.read_text(encoding="utf-8").startswith(
        "# Species-wide flower-colour variation is geographically partitioned across plant species"
    )
    assert "149,900" in s and "1,499" in s
    assert "genetically verified" in s
    assert "global between-species moisture analyses are outside this study" in s


def test_four_bullet_abstract_and_one_full_paper_figure_plan():
    s=PAPER.read_text(encoding="utf-8")
    summary=_section(s,"## Summary","## Introduction").split("**Keywords")[0].strip()
    assert len([row for row in summary.splitlines() if row.startswith("- ")])==4
    assert 130<=len(summary.split())<=200,len(summary.split())
    assert len(re.findall(r"^\*\*Figure [1-5]\.",s,flags=re.M))==5
    assert len(re.findall(r"^## Table 1\.",s,flags=re.M))==1
    assert len(re.findall(r"^\*\*Keywords",s,flags=re.M))==1


def test_submission_length_and_methods_pattern_mechanism():
    s=PAPER.read_text(encoding="utf-8")
    sections=[
        _section(s,"## Introduction","## Materials and Methods"),
        _section(s,"## Materials and Methods","## Table 1."),
        _section(s,"## Results","## Discussion"),
        _section(s,"## Discussion","## Acknowledgements"),
    ]
    n=sum(len(x.strip().split()) for x in sections)
    assert 6000<=n<=7500,n
    assert f"Main text (Introduction through Discussion): {n} words" in s
    assert "### Method: reproducible distributional colour measurement" in s
    assert "### Pattern: a frozen achromatic–chromatic axis is prospectively confirmed" in s
    assert "### Biological hypotheses for recurrent achromatic–chromatic variation" in s
    assert "### An ecological experiment for the water-balance hypothesis" in s


def test_actual_prospective_h2_is_158_not_439_confirmations():
    s=PAPER.read_text(encoding="utf-8")
    data=json.loads(H2.read_text())
    cohorts=json.loads(H2COHORTS.read_text())
    assert data["decision"]["verdict"]=="H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED"
    primary=next(x for x in cohorts["cohorts"] if x["tier"]=="primary")
    assert primary["n_eligible"]==439
    assert [q["n"] for q in primary["cohorts"]]==[152,129,158]
    assert primary["pooled_p_value"] is None
    assert math.isclose(primary["cohorts"][2]["W"],0.5172457461053418,abs_tol=1e-12)
    assert "439 species" in s
    assert "158 newly sampled species" in s
    assert "descriptive" in s
    assert "not prospectively confirmed by these same cohorts" in s


def test_spatial_photographic_pattern_not_population_fitness():
    s=PAPER.read_text(encoding="utf-8")
    j=json.loads(SPATIAL.read_text())
    assert j["primary_radius_km"]==50
    for k,n,approx in (("discovery",166,0.02053),("validation",181,0.01867),("third",204,0.01468)):
        v=j["primary"][k]
        assert v["n_species"]==n
        assert abs(v["mean_local_depletion"]-approx)<1e-5
        assert v["p_upper"]==0.005
    assert "0.02053" in s and "0.01867" in s and "0.01468" in s
    assert "genetic" in s.lower() and "fitness" in s.lower()
    assert "not a local allele-frequency" in PLAN.read_text() or "not a local" in PLAN.read_text()


def test_macroscopic_inter_specific_rainfall_does_not_fake_within_species_morphs():
    s=PAPER.read_text(encoding="utf-8")
    text_body=_section(s,"## Introduction","## Acknowledgements")
    assert "20,694" not in text_body
    assert "+2.255" not in text_body
    assert "20,786" not in text_body
    assert "not an effect estimated from the separate large between-species rainfall atlas" in text_body
    assert "within-species" in text_body
    p=PLAN.read_text()
    assert "zuizui0223/mois" in p
    assert "not a genetic" not in p.lower() or "genetic" in p.lower()


def test_technical_white_measurement_boundary_unchanged():
    s=PAPER.read_text(encoding="utf-8")
    assert "OR 1.444" in s
    assert "INDETERMINATE" in s
    assert "89.9%" in s
    assert "not a single conserved pigment mutation" in s or "not evidence" in s


def test_one_cover_letter_asks_three_research_questions():
    s=COVER.read_text(encoding="utf-8")
    assert "within-species" in s.lower()
    assert "158-species" in s
    assert "20,694" not in s
    for num in (1,2,3):
        title="**"+str(num)+". "
        assert title in s
    assert "[CORRESPONDING AUTHOR" in s
    assert "No submission has occurred" in s



def test_table1_and_figure1_keep_endpoint_specific_species_denominators():
    """Only 158 tested the frozen axis prospectively; 439 is not a new pooled test."""
    s=PAPER.read_text(encoding="utf-8")
    table=_section(s,"## Table 1.","## Results")
    legend=_section(s,"**Figure 1.","**Figure 2.")
    high=json.loads(H2.read_text())
    dat=json.loads(H2COHORTS.read_text())
    spatial=json.loads(SPATIAL.read_text())
    p=next(x for x in dat["cohorts"] if x["tier"]=="primary")
    assert high["measurement_support"]["evaluable_species"]==377
    assert p["n_eligible"]==439
    assert sum(x["n"] for x in p["cohorts"][:2])==281
    assert p["cohorts"][2]["n"]==158
    assert p["pooled_p_value"] is None
    assert "369 discovery + 363 validation + 377 prospective = 1,109" in table
    assert "152 discovery + 129 validation + 158 prospective = 439" in table
    assert "pre-frozen confirmatory test n=158" in table
    assert "166 discovery / 181 validation / 204 prospective" in table
    assert "not a sequential 1,499" in table
    assert all(spatial["primary"][key]["n_species"]==n for key,n in
               (("discovery",166),("validation",181),("third",204)))
    assert "439 eligible W vectors" in legend
    assert "No source-wide polymorphism prevalence is estimated" in legend



def test_submission_readiness_identifies_real_limits_and_requires_human_authority():
    note=ROOT/"docs/POLYMORPHISM_NP_INTRA_SUBMISSION_READINESS_20261014.md"
    assert note.exists()
    t=note.read_text(encoding="utf-8")
    for phrase in ("1,499", "1,109", "439", "158", "166", "181", "204", "149,900", "4,730"):
        assert phrase in t,phrase
    assert "not 439 preregistered confirmed species" in t
    assert "genetically inherited" not in t or "genetically inherited" in t.lower()
    assert "INDETERMINATE" in t
    assert "not a direct measurement" in t or "unmeasured" in t
    assert "No manuscript has been submitted" in t
