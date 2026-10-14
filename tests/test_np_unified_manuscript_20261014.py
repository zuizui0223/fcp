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
        stop=(s.index("## Table 1.",start) if name=="Materials and Methods" else s.index("## "+names[i+1],start)) if i+1<len(names) else s.index("## Acknowledgements",start)
        chunks[name]=s[start:stop]
    return s,chunks


def test_one_manuscript_primary_title_and_source_untouched():
    assert ORIGINAL.exists()
    orig=ORIGINAL.read_text(encoding="utf8")
    current,chunks=sections()
    assert orig.startswith("# Species-wide flower-colour variation is geographically partitioned across plant species")
    assert current.startswith("# Moisture-associated floral colours and recurrent achromatic–chromatic variation across plant species")
    assert "149,900" in current and "100,543" in current
    assert "not a fourth" in current
    assert "not an additional strong causal or confirmatory paper" not in chunks["Discussion"]


def test_journal_four_bullet_summary_under_200_words_and_seven_display_items():
    full,sections_=sections()
    summary=sections_["Summary"].split("**Keywords")[0].strip()
    bullets=[line for line in summary.splitlines() if line.startswith("- ")]
    assert len(bullets)==4
    assert 0 < len(summary.split()) <= 200
    assert len(re.findall(r"^\*\*Figure [1-6]\.",full,flags=re.M))==6
    assert len(re.findall(r"^## Table [1-9]\.",full,flags=re.M))==1
    assert "Figure 6" in full and "### Global colour–environment prediction" in full
    assert len(re.findall(r"\*\*Fig\. S1[2-5]\.",full))==4


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
    summary=c["Summary"].split("**Keywords")[0]
    for fragment in ("distributional photographic framework", "achromatic–chromatic direction", "annual precipitation", "water-vapour", "ecological trade-offs"):
        assert fragment in summary, fragment
    assert "universal climatic selection" not in summary.lower()
    assert "not robustly supported" not in summary.lower()


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
    assert "benefit–cost hypothesis" in c["Discussion"]
    assert "pigment synthesis" in c["Discussion"] and "do **not** directly measure" in c["Discussion"]
    assert "moisture" in c["Discussion"].lower()
    ledger=CLAIMS.read_text()
    for word in ("42,111","149,900","45,953","INDETERMINATE"):
        assert word in ledger,word
    assert "pigment costs and morph-specific fitness were not measured" in ledger.lower()


def test_cover_letter_is_one_submission_and_50_word_questions():
    s=COVER.read_text()
    assert "single Full Paper" in s
    for number in (1,2,3):
        m=re.search(r"^\*\*"+str(number)+r"\. .+?\*\*\s*\n(.+?)(?=\n\n)",s,flags=re.M|re.S)
        assert m is not None,number
        assert len(m.group(1).split())<=50,(number,len(m.group(1).split()))
    assert "[CORRESPONDING AUTHOR" in s



def test_figure6_vector_and_data_source_match_frozen_macro_results():
    data_path=ROOT/"results/np_unified_20261014/figure6_source_data.json"
    svg=ROOT/"docs/figures/polymorphism_20261014/figure6_unified_environmental_support.svg"
    render=ROOT/"scripts/analysis/render_np_unified_figure6_20261014.py"
    assert svg.exists() and render.exists() and svg.stat().st_size>100000
    svg_text=svg.read_text(encoding="utf8")
    assert "<svg " in svg_text and "</svg>" in svg_text
    ledger=json.loads(data_path.read_text(encoding="utf8"))
    assert ledger["schema"]=="fcp_np_unified_figure6_source_20261014_v1"
    assert set(ledger["panels"])=={"a","b","c","d"}
    assert sum(len(v["items"]) for v in ledger["panels"].values())==9
    for panel in ledger["panels"].values():
        for item in panel["items"]:
            assert item["lo"]<=item["point"]<=item["hi"]
    env=json.loads(MACRO.read_text())
    def env_row(quality,cell):
        row=[x for x in env["key_white_directional_observations"]
             if x["mode"]=="climate" and x["region"]=="world" and
                x["quality"]==quality and x["adjustment"]==("genus_geocell" if cell else "genus_geo") and
                x["axis"]=="wc_bio12" and x["hue"]=="white"]
        assert len(row)==1
        return row[0]
    for i,(quality,cell) in enumerate([("all",False),("strict",False),("all",True)]):
        r=env_row(quality,cell);p=ledger["panels"]["b"]["items"][i]
        assert math.isclose(p["point"],r["point_pp"],abs_tol=1e-9)
        assert all(math.isclose(p[k],r["family_95_interval_pp"][j],abs_tol=1e-9) for j,k in enumerate(("lo","hi")))
    sp=json.loads(SPATIAL.read_text())
    for i,key in enumerate(["world|genus_equal","E90_E180|genus_equal"]):
        r=sp["BIO12_white_by_region_and_weight"][key]
        p=ledger["panels"]["c"]["items"][i]
        assert math.isclose(p["point"],r["paired_adjustment_difference_pp"],abs_tol=1e-9)
        assert all(math.isclose(p[k],r["paired_24_hue_axis_maxT_95_interval_pp"][j],abs_tol=1e-9) for j,k in enumerate(("lo","hi")))
    nhm=json.loads(EXTERNAL.read_text())
    for i,key in enumerate(["equal_species","equal_genera"]):
        r=nhm["cases"]["min3_all_photo_opportunity"][key]["east_result"]
        p=ledger["panels"]["d"]["items"][i]
        assert math.isclose(p["point"],100*r["point"],abs_tol=1e-9)
        assert all(math.isclose(p[k],100*r["CI95"][j],abs_tol=1e-9) for j,k in enumerate(("lo","hi")))


def test_single_unified_supporting_information_includes_all_four_macro_layers():
    target=ROOT/"docs/POLYMORPHISM_SUPPORTING_INFORMATION_UNIFIED_20261014.md"
    assert target.exists()
    si=target.read_text(encoding="utf8")
    for heading in ("## S13.","## S14.","## S15.","## S16."):
        assert heading in si,heading
    manuscript=MANUSCRIPT.read_text(encoding="utf8")
    assert target.name in manuscript



def test_climate26_vs_all36_species_counts_are_distinct_and_source_verified():
    """Do not subtract 20786 - 20694 as if one filter removed 92 taxa."""
    r=json.loads(MACRO.read_text())
    observations=r["key_white_directional_observations"]
    def matching(mode):
        q=[x for x in observations if
           x["mode"]==mode and x["quality"]=="all" and x["region"]=="world"
           and x["adjustment"]=="genus_geo" and x["axis"]=="wc_bio12" and x["hue"]=="white"]
        assert len(q)==1
        return q[0]
    climate=matching("climate")
    all36=matching("soil")
    assert (climate["n_photo_species"],climate["n_genera"])==(20694,1934)
    assert (all36["n_photo_species"],all36["n_genera"])==(16801,1665)
    assert 25028-20694==4334
    assert 20786-16801==3985
    ms=MANUSCRIPT.read_text(encoding="utf8")
    memo=(ROOT/"docs/POLYMORPHISM_20K_SCALING_SUPPORT_20261014.md").read_text(encoding="utf8")
    for count in ("25,028","20,694","20,786","16,801","1,934","1,665"):
        assert count in ms and count in memo,count
    assert "20,786 − 20,694 = 92 is *not an attrition step*" in memo


def test_ecology_first_story_and_noncausal_separate_sampling_grains():
    s, c = sections()
    r, d = c["Results"], c["Discussion"]
    assert r.index("### Ecology: worldwide flower-colour correspondence") < r.index("### Pattern: achromatic–chromatic variation")
    assert r.index("### Pattern: achromatic–chromatic variation") < r.index("### Geographic context: local mixing")
    assert d.index("### Moisture-associated ecology") < d.index("### Recurring achromatic–chromatic patterns")
    assert d.index("### A water-dependent benefit–cost hypothesis") < d.index("### Conclusion")
    assert "nonnested sampling designs" in c["Introduction"]
    assert "not a nested" in c["Summary"]
    assert "strict" in c["Summary"] and "spatial" in c["Summary"]
    assert "relative fitness" in d and "do **not** directly measure" in d
    assert "figures still retain original 1–6 numbering" in s
