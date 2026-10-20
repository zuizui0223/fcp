#!/usr/bin/env python3
"""Fail-closed H1..H5 minimal evidence audit; no synthetic successes or new causal claims.

Reads independently frozen actual measured-photo analysis receipts and
literature/molecular metadata; outputs an evidence ledger. References are
pinned through the runner workflow to the main image source and two branches.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path

def require(ok,msg):
    if not ok: raise RuntimeError("FCP_H1H5_MINIMUM_EVIDENCE_GATE_FAILED: "+msg)

def read(p):
    return json.loads(Path(p).read_text())

def check(h1,clip,ibd,coex,cs,pal,molecular,nested):
    # H1 is a signed frozen inference, different from the conditional
    # categorical white frequency association used in H2/H3.
    require(h1["decision"]["verdict"]=="H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED",
            "prospective white/chromatic W no longer supported")
    x=h1["thresholds"]["primary_0_10"]
    require(x["species"]==158 and x["observed_W"]>x["structured_null_summary"]["q50"]
            and x["structured_null_upper_p"]==.001,"original H1 evidence changed")
    require(clip["decision"]["state"]=="INDETERMINATE"
            and clip["h2_high_clip_sensitivity"]["vector_retention"]<.9,
            "unresolved image-exposure gate must be shown")

    # Phenotypic distances are not neutral SNP distances. Geography is not history.
    c=ibd["cohorts"]
    require(all(c[k]["IBD_like"]["mean_partial_rho"]>c[k]["IBE_like"]["mean_partial_rho"]
                for k in ("discovery","validation","third")),
            "old three-cohort IBD-like/IBE-like order changed")
    require("not genetic isolation by environment" in ibd["hard_nonclaims"],
            "genetic-IBD/IBE boundary missing")

    # H3 is conditional on photo quality and long-term BIO5, not a
    # temporal or local-fitness model. Algebraically exact nested split.
    require(nested["source_run"].endswith("/39430364175"),
            "nested climate result is not source-authenticated")
    nr=nested["original_result_from_production_job"]["strict"]["all_eligible_5deg_month_strata"]
    t=nr["mean_signed_white_minus_nonwhite_C"]
    require(nr["n_species"]==327
            and abs(t["total_5deg_month"]["BIO5"]-
                    t["within_1deg_month"]["BIO5"]-
                    t["between_1deg_month"]["BIO5"])<1e-8,
            "real nested BIO5 decomposition not exact")
    require(nr["one_degree_nested_conditioned_within_component_upper_p"]["BIO5"]==.075,
            "strict within1 signal status changed")

    # H4 is local populations with inherited morphs, not regional photos.
    q=coex["matched_month_null"]["min_5_pairs"]
    require(q["n_species"]==194 and q["upper_p"]==.30,
            "same-site month null changed")
    require(cs["totals"]["spatial_segregation_only_spatial_matches"]==0,
            "C/S bridge availability changed, refit required")

    # H5 molecular evidence is literature-derived, without actual photo-genotypes.
    require(pal["PAL"]["instances"]==13 and pal["WAL"]["instances_total"]==13,
            "PAL/WAL literature ascertainment denominator changed")
    require(molecular["union_curated_species"]==14
            and molecular["union_curated_contrasts"]==24
            and molecular["cross_species_mechanism_environment_fit"] is False,
            "molecular comparative estimability changed")

    strict_total=t["total_5deg_month"]["BIO5"]
    between=t["between_1deg_month"]["BIO5"]
    outputs=[
      {
       "hypothesis":"H1_common_colour_axis",
       "reported_test":"Frozen species-disjoint W against construction-preserving nine-colour palette null",
       "n_taxa":x["species"],"estimate":x["observed_W"],
       "null_reference":x["structured_null_summary"]["q50"],
       "permutation_p":x["structured_null_upper_p"],
       "status":"SUPPORTED_WITHIN_ONE_PHOTO_SYSTEM",
       "technical_caveat":"near-clip image coupling; high-clip retention 142/158 = 89.87%, prospectively selected support gate INDETERMINATE",
       "minimum_new_evidence":"independent-photographer, colour-calibrated full-flower reflectance in 20 species; selection and null rules frozen before colours opened",
       "causal_history_claim":False
      },
      {
       "hypothesis":"H2_historical_range_filter",
       "reported_test":"Phenotypic partial rho(colour, distance | BIO5) versus rho(colour, BIO5 | distance) across three cohorts",
       "n_taxa":sum(c[k]["n_species"] for k in ("discovery","validation","third")),
       "ibd_like_by_cohort":[c[k]["IBD_like"]["mean_partial_rho"] for k in ("discovery","validation","third")],
       "ibe_like_by_cohort":[c[k]["IBE_like"]["mean_partial_rho"] for k in ("discovery","validation","third")],
       "status":"HISTORICAL_CAUSE_NOT_IDENTIFIED",
       "minimum_new_evidence":"neutral SNP ancestry and independently informed historic connectivity in the SAME populations; matched current environment and distance",
       "causal_history_claim":False
      },
      {
       "hypothesis":"H3_environment_sorting",
       "reported_test":"Real original-source same-month nested 5 degree BIO5 white/photo contrast, split within/between 1 degree",
       "n_taxa":nr["n_species"],"photo_BIO5_total_C":strict_total,
       "photo_BIO5_between_1deg_C":between,
       "photo_BIO5_within_1deg_C":t["within_1deg_month"]["BIO5"],
       "between_fraction_of_signed_photo_BIO5_association_descriptive":between/strict_total,
       "between_1deg_permutation_p":nr["five_degree_vertex_null_upper_p"]["between_1deg_month"]["BIO5"],
       "within_1deg_conditional_permutation_p":nr["one_degree_nested_conditioned_within_component_upper_p"]["BIO5"],
       "status":"REGIONAL_PHENOTYPIC_ASSOCIATION_FINE_LOCAL_NOT_SUPPORTED",
       "minimum_new_evidence":"same 2 focal species and population-genotype backgrounds, matched local microclimate/soil and morph-by-treatment seed fitness",
       "causal_environmental_selection_claim":False
      },
      {
       "hypothesis":"H4_local_morph_maintenance",
       "reported_test":"5-km same-year-month independent-photo matched null and strict literature C/S to photograph spatial crosswalk",
       "n_taxa":q["n_species"],"local_photo_month_p":q["upper_p"],
       "strict_literature_c_s_species":cs["totals"]["literature_species"],
       "literature_spatial_photo_overlap":cs["totals"]["spatial_50km"],
       "spatial_segregation_only_photo_overlap":cs["totals"]["spatial_segregation_only_spatial_matches"],
       "status":"UNTESTABLE_CURRENT_INHERITANCE_AND_FITNESS_DATA",
       "minimum_new_evidence":"mark individual plants for 2 years, genotype morph and measure morph-specific flowering, pollination and seed fitness within delimited interbreeding populations",
       "genetic_local_polymorphism_confirmed":False
      },
      {
       "hypothesis":"H5_tissue_specific_molecular_reachability",
       "reported_test":"Published PAL/WAL natural frequency and curated gene/tissue contrast evidence",
       "PAL_literature_cases":pal["PAL"]["instances"],
       "WAL_literature_cases":pal["WAL"]["instances_total"],
       "PAL_lower_median_percent":pal["PAL"]["lower_bound_median_percent"],
       "WAL_numeric_upper_median_percent":pal["WAL"]["numeric_upper_bound_median_percent"],
       "curated_distinct_species":molecular["union_curated_species"],
       "curated_contrasts_not_independent_origins":molecular["union_curated_contrasts"],
       "status":"LITERATURE_CANDIDATE_NOT_CAUSALLY_GENERALIZED",
       "minimum_new_evidence":"same field-genotyped individuals: calibrated petal+leaf anthocyanin LC-MS, targeted F3H/MYB/ANS expression, common-garden and maternal seed fitness",
       "global_white_molecular_cause_established":False
      }
    ]
    return {"schema":"fcp_h1_h5_minimum_discriminating_verification_v1",
      "date":"2026-10-20","source":"existing source-locked evidence, not a new H1-H5 causal experiment",
      "n_photos_in_original_measurement":447500,
      "user_goal":"mechanism-led global within-species floral colour diversity, rather than mere spatial aggregation",
      "tests":outputs,
      "decision":"H1_photo_space_supported; H2_history_unknown; H3_regional_current_climate_association_only; H4_populational_maintenance_unknown; H5_molecular_candidate_only",
      "do_not_combine_into_universal_selection_claim":True
    }

def main():
    p=argparse.ArgumentParser()
    for name in ("h1","clip","ibd","coex","cs","pal","molecular","nested"):
        p.add_argument("--"+name,required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    args={k:read(getattr(a,k)) for k in ("h1","clip","ibd","coex","cs","pal","molecular","nested")}
    result=check(**args)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("FCP_H1H5_MINIMAL_AUDIT_VERIFIED",json.dumps({
      "status":[(r["hypothesis"],r["status"]) for r in result["tests"]],
      "H3_between_fraction":result["tests"][2]["between_fraction_of_signed_photo_BIO5_association_descriptive"],
      "H4_S_only_bridge":result["tests"][3]["spatial_segregation_only_photo_overlap"],
      "H5_curated_species":result["tests"][4]["curated_distinct_species"]},sort_keys=True,allow_nan=False))

if __name__=="__main__":main()
