#!/usr/bin/env python3
"""Fail-closed New Phytologist global FCP manuscript claim and format audit.

This is an editorial/source audit of pre-existing authentic empirical receipts,
not new species, photo, genotype, fitness or transcriptomic inference. The
manuscript has not been submitted and graphics are not re-rendered by this test.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
MANUSCRIPT=ROOT/"docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md"
PHOTO_SOURCE_SHA="8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6"


def require(condition,reason):
    if not condition:
        raise RuntimeError("FCP_NEW_PHYTOLOGIST_GLOBAL_AUDIT_FAILED: "+reason)


def readjson(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def readcsv(path):
    with Path(path).open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))


def approx(x,y,tol=5e-5):
    return abs(float(x)-y)<=tol


def between(s,start,end):
    a=s.find(start);b=s.find(end,a+len(start))
    require(a>=0 and b>a,"missing or reordered manuscript section: "+start)
    return s[a:b]


def audit_manuscript(md):
    lines=md.splitlines()
    require(lines[0].startswith("# "),"main title absent")
    title=lines[0][2:]
    require(len(title)<=130,"Full Paper title exceeds 130 characters")
    top=[
        "## Summary","## Introduction","## Materials and Methods",
        "## Results","## Discussion","## Acknowledgements",
        "## Competing interests","## Author contributions",
        "## Data availability","## References","## Figure legends",
        "## Supporting Information"
    ]
    location=[md.find(h+"\n") for h in top]
    require(all(p>=0 for p in location) and location==sorted(location),
            "New Phytologist Full Paper required section ordering")
    summ=between(md,"## Summary","## Introduction")
    bullets=[s for s in summ.splitlines() if s.startswith("- ")]
    require(len(bullets)==4,"Full Paper Summary must contain exactly four research bullets")
    summary_words=len(" ".join(bullets).split())
    require(summary_words<=200,f"four-bullet summary exceeds 200 words: {summary_words}")
    keyline=next((x for x in lines if x.startswith("**Keywords (alphabetical):** ")),None)
    require(keyline is not None,"keywords missing")
    words=[x.strip() for x in keyline.split("** ",1)[1].split(";")]
    require(5<=len(words)<=8 and words==sorted(words,key=str.casefold),
            "keywords must be 5–8 alphabetical phrases")
    wc={title:len(between(md,"## "+title,"## "+after).split())
        for title,after in [
          ("Introduction","Materials and Methods"),
          ("Materials and Methods","Table 1."),
          ("Results","Discussion"),
          ("Discussion","Acknowledgements"),
        ]}
    main_words=sum(wc.values())
    require(6500<=main_words<=7500,
            f"New Phytologist Full Paper main text outside 6500–7500 working target: {main_words}")
    figs=between(md,"## Figure legends","## Supporting Information")
    require(len(re.findall(r"^\*\*Figure [1-5]\.",figs,re.M))==5,
            "planned primary five figure captions not intact")
    require(md.count("## Table 1.")==1,"one source/candidate Table 1 required")
    # Never let separate hypothesis-specific n be described as an attrition
    # funnel or as three independently collected data providers.
    src_table=between(md,"## Table 1.","## Results")
    for raw in ("4,475 species / 447,500 photographs",
                "499 new species → 377 photo-classifiable → 185 coarse-second-morph-eligible → 158 continuous W vectors",
                "1,761 species","1,574 species / 86,114 photographs"):
        require(raw in src_table,"missing distinct-source-panel eligibility: "+raw)
    require("distinct eligibility subsets" in src_table
            and "not three consecutive attrition stages" in src_table,
            "New Phytologist source panel mistakenly implied nested design")
    for heading in ["### Worldwide photographic variation","### Prospective H2: the frozen white axis is confirmed",
                    "### Geographic sorting across all 447,500",
                    "### Weak global climatic direction does not predict novel regions reliably",
                    "### H4 and H5: local maintenance and pigment modularity are separate tests of persistence"]:
        require(heading in md,"evolutionary evidence narrative missing "+heading)
    for caveat in ["not a global genetic-polymorphism census",
                   "INDETERMINATE","not disproportionately",
                   "p = **0.300**","p = **0.220**",
                   "p = **0.075**","**not** a direct measurement",
                   "unmeasured causal mechanisms"]:
        require(caveat in md,"non-support / inferential boundary missing: "+caveat)
    require("**Authors:** [AUTHOR LIST TO CONFIRM]" in md,
            "unconfirmed author list may not be fabricated")
    return {"title":title,"title_characters":len(title),"four_summary_bullets":len(bullets),
            "summary_words":summary_words,"key_word_count":len(words),"body_word_counts":wc,
            "body_total_words":main_words,"planned_main_figures":5,"main_tables":1,
            "illustrations_still_require_regeneration":True,
            "authors_and_affiliations_not_finalized":True}


def audit_source(h2,locality,white,clim,transfer,atlas,world_white,nested):
    axis=h2["thresholds"]["primary_0_10"]
    require(h2["decision"]["verdict"]=="H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED","H1 old frozen axis changed")
    require(axis["species"]==158 and approx(axis["observed_W"],.5172457461,1e-9)
            and approx(axis["structured_null_summary"]["q50"],.457142815,1e-9)
            and axis["structured_null_upper_p"]==.001,
            "prospective axis receipt inconsistent")
    pooled=white["pooled_effect_decomposition"]
    require(approx(pooled["total"],.0176149021685411,1e-10)
            and approx(pooled["white_nonwhite_over_positive_pooled_total"],.8381108938332189,1e-9),
            "actual pooled 50-km matched source changed")
    cohorts={r["group"]:r for r in locality if r["statistic"]=="mean"}
    require(len(cohorts)==6 and "pooled" in cohorts,"5 original cohorts + pooled 50-km receipt not found")
    require(int(cohorts["pooled"]["n_species"])==1761
            and approx(cohorts["pooled"]["observed"],.0176149022,1e-9),
            "original 1761 species locality source drift")
    native=next((r for r in clim if r["group"]=="pooled"
                 and r["variable"]=="BIO5" and r["model"]=="primary"),None)
    require(native is not None and int(native["n_species"])==1574
            and approx(native["equal_species_mean"],.0323286448565001,1e-10)
            and approx(native["genus_cluster_percentile_low"],.017855650726,1e-9)
            and approx(native["genus_cluster_percentile_high"],.04727667440,1e-9),
            "exploratory signed BIO5 association differs")
    spatial={r["comparison"]:r for r in transfer if float(r["buffer_km"])==0}
    require(approx(spatial["climate_vs_intercept"]["logloss_gain"],-.05992424800536707,1e-10)
            and approx(spatial["geoclimate_vs_geography"]["logloss_gain"],-.14977710772409594,1e-10),
            "spatial holdout outcome drift")
    world=atlas["world_photo_atlas"]
    require(atlas["source_original_photo_csv_sha256"]==PHOTO_SOURCE_SHA,
            "external world photo atlas source SHA mismatch")
    require(world["original_measured_rows"]==447500
            and world["original_measured_species"]==4475
            and world["original_classifiable_four_state_rows"]==227107
            and world["global_D_evaluable_species_min40photos_3observers"]==3297
            and world["global_species_with_2plus_detected_photo_colour_labels"]==3029,
            "worldwide photo eligibility census drift")
    require(world["n_cells_with_eligible_taxon_units"]==65,"equal-area eligible cell count drift")
    require(len([c for c in world["equal_area_grid_cells"] if c["n_species_eligible"]>=20])==28,
            "world geographic reporting denominator drift")
    require(world["paired_trop_extratrop_same_species"]["n_spanning_both_zones"]==246,
            "trop-extratrop paired comparison denominator changed")
    wp=world_white["result"]["matched_same_species_same_cells_high_purity"]
    require(wp["n_species"]==1282
            and approx(wp["white_share_of_mean_excess_turnover"],.8946024398624991,1e-9),
            "strict source composition site matching changed")
    normalized=world_white.get("baseline_normalized_result_qualification",{})
    require(normalized and "not supported" in str(normalized).lower(),
            "explicit white-specific relative regional non-support missing")
    n=nested["original_result_from_production_job"]["strict"]["all_eligible_5deg_month_strata"]
    temp=n["mean_signed_white_minus_nonwhite_C"]
    require(n["n_species"]==327
            and approx(temp["total_5deg_month"]["BIO5"],.3298491039,1e-8)
            and approx(temp["between_1deg_month"]["BIO5"],.2715821283,1e-8)
            and approx(temp["within_1deg_month"]["BIO5"],.0582669756,1e-8)
            and n["one_degree_nested_conditioned_within_component_upper_p"]["BIO5"]==.075,
            "nested photo BIO5 thermal test drift")
    return {"photographs":447500,"species":4475,"classifiable":227107,
        "W_158":axis["observed_W"],"W_null_median":axis["structured_null_summary"]["q50"],
        "geo_mean_1761":float(cohorts["pooled"]["observed"]),
        "BIO5_equal_mean_log_odds_per_C":float(native["equal_species_mean"]),
        "spatial_holdout_climate_gain":float(spatial["climate_vs_intercept"]["logloss_gain"]),
        "strict_nested_taxa":n["n_species"],
        "strict_nested_within_one_degree_BIO5_p":.075,
        "source_boundaries":"photographic patterns; no genetic history, fitness or tissue intervention"}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--atlas",type=Path,required=True)
    ap.add_argument("--world-white",type=Path,required=True)
    ap.add_argument("--nested",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    man=MANUSCRIPT.read_text(encoding="utf-8")
    try:
        review=audit_manuscript(man)
        values=audit_source(
            readjson(ROOT/"results/polymorphism_h2_third_cohort_prospective_white_axis_20260917/result.json"),
            readcsv(ROOT/"results/fcp_real_447500_robust_matched_null_20261018/group_statistic_summary.csv"),
            readjson(ROOT/"results/fcp_real_447500_white_axis_decomposition_20261018/result.json"),
            readcsv(ROOT/"results/fcp_environment_native_units_20261019/primary_summary.csv"),
            readcsv(ROOT/"results/fcp_environment_transport_20261019/spatial_summary.csv"),
            readjson(args.atlas),readjson(args.world_white),readjson(args.nested))
    except Exception as e:
        raise RuntimeError("Full Paper/source QA halted before success") from e
    obj={
        "schema":"fcp_np_global_evolutionary_manuscript_verified_20261021_v1",
        "manuscript_sha256":hashlib.sha256(man.encode("utf-8")).hexdigest(),
        "verified_main_authoring_sections":review,
        "verified_original_data_receipts":values,
        "submission_status":"WORKING_FULL_PAPER_NOT_SUBMITTED_FIGURES_AND_AUTHOR_DETAILS_PENDING",
        "frozen_original_analysis_modified":False,
        "current_original_peer_review_claims_not_strengthened_by_editorial_qa":True,
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+"\n")
    print("FCP_NEW_PHYTOLOGIST_GLOBAL_MAIN_PAPER_SOURCE_AND_FORMAT_QA_PASS",
          json.dumps({"words":review["body_total_words"],"summary_words":review["summary_words"],
                      "figures":review["planned_main_figures"],"real_source_photos":values["photographs"],
                      "status":obj["submission_status"]},sort_keys=True))


if __name__=="__main__":
    main()
