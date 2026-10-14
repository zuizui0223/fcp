"""One-paper New Phytologist editorial/claim guard; never changes frozen science."""
from pathlib import Path
import json
import math
import re

ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST_UNIFIED_20261014.md"
ORIGINAL = ROOT / "docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md"
CLAIMS = ROOT / "docs/POLYMORPHISM_UNIFIED_EVIDENCE_AND_CLAIM_LEDGER_20261014.md"
COVER = ROOT / "docs/POLYMORPHISM_NEW_PHYTOLOGIST_UNIFIED_COVER_LETTER_20261014.md"
MACRO = ROOT / "results/fcp_macro_direction_20261013/decision_receipt.json"
SPATIAL = ROOT / "results/fcp_m1_matched_spatial_support_20261014/decision_receipt.json"
EXTERNAL = ROOT / "results/np_unified_20261014/nhm_source_verified_colour_summary.json"


def sections():
    s=MANUSCRIPT.read_text(encoding="utf8")
    names=("Summary","Introduction","Materials and Methods","Results","Discussion")
    chunks={}
    for i,name in enumerate(names):
        start=s.index("## "+name)+len("## "+name)
        stop=s.index("## "+names[i+1],start) if i+1<len(names) else s.index("## Acknowledgements",start)
        chunks[name]=s[start:stop]
    return s,chunks


def test_one_manuscript_primary_title_and_source_untouched():
    assert ORIGINAL.exists()
    orig=ORIGINAL.read_text(encoding="utf8")
    current,chunks=sections()
    assert orig.startswith("# Species-wide flower-colour variation is geographically partitioned across plant species")
    assert current.startswith("# Geographic partitioning of flower-colour variation is widespread, while climatic associations are scale dependent")
    assert "149,900" in current and "100,543" in current
    assert "not a fourth" in current
    assert "not an additional strong causal or confirmatory paper" not in chunks["Discussion"]


def test_journal_four_bullet_summary_under_200_words_and_seven_display_items():
    full,sections_=sections()
    summary=sections_["Summary"].split("**Keywords")[0].strip()
    bullets=[line for line in summary.splitlines() if line.startswith("- ")]
    assert len(bullets)==4
    assert 0 < len(summary.split()) <= 200
    assert len(re.findall(r"^\*\*Figure [1-6]\\.",full,flags=re.M))==6
    assert len(re.findall(r"^## Table [1-9]\\.",full,flags=re.M))==1
    assert "Figure 6" in full and "### Global colour–environment prediction" in full
    assert len(re.findall(r"\*\*Fig\\. S1[2-5]\\.",full))==4


def test_np_word_target_and_discussion_proportion():
    full,c=sections()
    main=sum(len(c[x].strip().split()) for x in ("Introduction","Materials and Methods","Results","Discussion"))
    assert 6500<=main<=7500,main
    assert len(c["Discussion"].split())/main<=0.30
    assert "Main text (Introduction through Discussion): "+str(main)+" words" in full


def test_original_high_depth_primary_effects_and_h2_kept():
    full,c=sections()
    for token in ("0.02053","0.01867","0.01468","0.517","0.457","1.444","INDETERMINATE"):
        assert token in full,token
    assert "not a within-population" in full or "not a measure of coexistence within a natural population" in full
    assert "genetic maintenance" in c["Summary"]


def test_macro_numeric_claims_backed_by_original_frozen_receipts():
    full,c=sections()
    photo=json.loads(MACRO.read_text())
    spatial=json.loads(SPATIAL.read_text())
    assert photo["n_macro_scope_cohorts"]==20
    assert photo["directional_estimates_and_sensitivities"]==1280
    q=[x for x in photo["key_white_directional_observations"]
       if (x["mode"],x["quality"],x["region"],x["adjustment"],x["axis"],x["hue"]) ==
       ("climate","all","world","genus_geo","wc_bio12","white")]
    assert len(q)==1
    assert math.isclose(q[0]["point_pp"],2.254545235123587,abs_tol=1e-10)
    assert q[0]["family_95_interval_pp"][0]>0
    assert spatial["n_environmental_direction_rows"]==576
    east=spatial["BIO12_white_by_region_and_weight"]["E90_E180|genus_equal"]
    assert east["n_matched_species"]==2886
    assert not east["difference_excludes_zero_after_24_maxT"]
    assert east["paired_24_hue_axis_maxT_95_interval_pp"][0]<0<east["paired_24_hue_axis_maxT_95_interval_pp"][1]
    assert "−0.258" in c["Results"]


def test_external_botanical_colour_not_claimed_as_independent_geography():
    s,c=sections()
    nhm=json.loads(EXTERNAL.read_text())
    assert nhm["original_overlap"]["botanically_labeled_species"]==2476
    x=nhm["cases"]["min3_all_photo_opportunity"]["equal_genera"]
    assert x["east"]["support"]["n_genus"]==58
    assert x["east"]["support"]["named_species"]==670
    assert math.isclose(x["east_result"]["point"],0.035550472136345004,abs_tol=1e-12)
    assert x["east_result"]["CI95"][0]<0<x["east_result"]["CI95"][1]
    assert "photo" in c["Results"].lower() and "site" in c["Materials and Methods"].lower()
    assert "not an independent high-depth" in c["Materials and Methods"] or "not a fourth" in c["Materials and Methods"]


def test_hard_selection_and_genetic_boundaries():
    s,c=sections()
    assert "not a causal spatial effect" in c["Materials and Methods"]
    assert "not measured" in c["Materials and Methods"] or "Neither data system estimates" in c["Materials and Methods"]
    assert "drift" in c["Discussion"] and "fitness" in c["Discussion"]
    assert "universal climatic" in s or "universal climatic" in c["Discussion"]
    ledger=CLAIMS.read_text()
    for word in ("42,111","149,900","45,953","INDETERMINATE","NOT_ESTABLISHED"):
        assert word in ledger,word


def test_cover_letter_is_one_submission_and_50_word_questions():
    s=COVER.read_text()
    assert "single Full Paper" in s
    for number in (1,2,3):
        m=re.search(r"^\*\*"+str(number)+r"\\. .+?\*\*\s*\n(.+?)(?=\n\n)",s,flags=re.M|re.S)
        assert m is not None,number
        assert len(m.group(1).split())<=50,(number,len(m.group(1).split()))
    assert "[CORRESPONDING AUTHOR" in s
