#!/usr/bin/env python3
"""Design feasibility from authentic FCP source rosters + published PAL/WAL frequencies.

No causal ecological conclusions. No genomic data or independent spectroscopy are
imputed. Historical white photo labels are NOT manually verified flower morphs.
Candidate taxa for independent floral-organ QA are fixed deterministically from
preexisting family and photographed colour-opportunity metadata; selection is
ascertained and cannot estimate prevalence.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

SEED="FCP_MINIMUM_PANEL_20261020_PREPUBLISHED_FAMILY_AND_OPPORTUNITY"
COHORTS=("discovery","validation")
PILOT_N=120
TARGET_FAMILIES=20
N_PER_COHORT=10
MIN_CLASSIFIED=40
MIN_EACH_LABEL=10
EXPECTED_ROSTER=4475
EXPECTED_CLASSIFIABLE=227107
EXPECTED_LITERATURE=34
EXPECTED_SPATIAL_OVERLAP=5
MAX_WAL=0.0086

def require(x,msg):
    if not x:
        raise RuntimeError("FCP_MINIMUM_PILOT_SOURCE_OR_DESIGN_ERROR: "+msg)

def read_csv(path):
    with Path(path).open("r",encoding="utf-8",newline="") as f:
        return list(csv.DictReader(f))

def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))

def rank_species(name):
    return hashlib.sha256(f"{SEED}|{name}".encode("utf-8")).hexdigest()

def p_ge_one(n,p):
    return -math.expm1(n*math.log1p(-p))

def p_ge_two(n,p):
    return max(0.,p_ge_one(n,p)-n*p*math.exp((n-1)*math.log1p(-p)))

def minimum_n(p,k=1,threshold=.95):
    fn=p_ge_one if k==1 else p_ge_two
    lo,hi=0,1
    while fn(hi,p)<threshold:
        hi*=2
        require(hi<=10**8,"required sample grows beyond scientific range")
    while hi-lo>1:
        mid=(hi+lo)//2
        if fn(mid,p)>=threshold:hi=mid
        else:lo=mid
    return hi

def power_scenarios(silene):
    require(silene["source"].find("Del Valle et al. 2019")>=0,"published frequency source changed")
    pal=silene["phenotypes"]["PAL"]
    wal=silene["phenotypes"]["WAL"]
    require(silene["surveyed_populations"]==21,"published populations changed")
    scenarios=[
       ("WAL_min_positive",wal["positive_frequency_min_percent"]),
       ("WAL_median_positive",wal["positive_frequency_median_percent"]),
       ("WAL_max_positive",wal["positive_frequency_max_percent"]),
       ("PAL_min_positive",pal["positive_frequency_min_percent"]),
       ("PAL_median_positive",pal["positive_frequency_median_percent"]),
       ("PAL_max_positive",pal["positive_frequency_max_percent"]),
    ]
    result=[]
    for label,percentage in scenarios:
        p=float(percentage)/100
        require(0<p<1,"published frequencies invalid")
        result.append({
           "name":label,
           "positive_population_year_frequency_percent":float(percentage),
           "n_random_individuals_in_one_species":PILOT_N,
           "p_at_least_one_if_independent_binomial":p_ge_one(PILOT_N,p),
           "p_at_least_two_if_independent_binomial":p_ge_two(PILOT_N,p),
           "n_to_have_95pct_at_least_one_binomial":minimum_n(p,1),
           "n_to_have_95pct_at_least_two_binomial":minimum_n(p,2),
        })
    med=next(r for r in result if r["name"]=="WAL_median_positive")
    require(.2<med["p_at_least_one_if_independent_binomial"]<.25
            and med["n_to_have_95pct_at_least_one_binomial"]>1000,
            "unexpected WAL ascertainment feasibility")
    return result

def select_independent_organ_qc_panel(family,photos):
    require(len(photos)==EXPECTED_ROSTER,"photo-opportunity roster denominator changed")
    unique={x["inat_taxon_id"] for x in photos}
    require(len(unique)==len(photos),"photo-opportunity roster contains duplicate taxon IDs")
    require(sum(int(x["n_classifiable"]) for x in photos)==EXPECTED_CLASSIFIABLE,
            "source photo classification denominator changed")
    family_by_species={}
    for x in family:
        if str(x["family_resolved"]).upper()!="TRUE" or not x["family"]:continue
        name=x["species"]
        require(name not in family_by_species or family_by_species[name]==x["family"],
                "same species assigned incompatible families")
        family_by_species[name]=x["family"]
    check=[]
    for p in photos:
        fam=family_by_species.get(p["species"])
        if not fam or p["source_cohort"] not in COHORTS:continue
        nc,nw,nv=(int(p[k]) for k in ("n_classifiable","n_white","n_nonwhite"))
        if nc<MIN_CLASSIFIED or nw<MIN_EACH_LABEL or nv<MIN_EACH_LABEL:continue
        require(nw+nv==nc,"white and nonwhite photographic counts inconsistent")
        check.append({
           "inat_taxon_id":int(p["inat_taxon_id"]),
           "species":p["species"],
           "family":fam,
           "source_cohort":p["source_cohort"],
           "n_classifiable_photo_labels":nc,
           "n_white_photo_labels":nw,
           "n_nonwhite_photo_labels":nv,
           "rank_hash":rank_species(p["species"]),
        })
    chosen=[];used_family=set()
    for cohort in COHORTS:
        cohort_cases=sorted((r for r in check if r["source_cohort"]==cohort),
                            key=lambda r:(r["rank_hash"],r["species"]))
        for r in cohort_cases:
            if r["family"] in used_family:continue
            chosen.append(r);used_family.add(r["family"])
            if sum(x["source_cohort"]==cohort for x in chosen)==N_PER_COHORT:
                break
        require(sum(x["source_cohort"]==cohort for x in chosen)==N_PER_COHORT,
                "insufficient family-disjoint taxa for cohort-balanced panel")
    require(len(chosen)==TARGET_FAMILIES and len(used_family)==TARGET_FAMILIES,
            "20-family independent measurement panel could not be selected")
    for r in chosen:
        r["pilot_status"]="SOURCE_PHOTO_LABEL_CANDIDATE_ONLY_NOT_TRUE_WILD_MORPH"
        r["minimum_independent_organ_measurements"]="two independently confirmed real morphs x at least 8 tagged individual flowers when available; no morphology imputations"
        r["should_replace_after_negative_review"]=False
    return sorted(chosen,key=lambda r:(r["source_cohort"],r["rank_hash"])),{
       "n_source_taxonomy_family_rows":len(family),
       "n_photo_opportunity_taxa":len(photos),
       "n_candidate_photo_labels_with_family":len(check),
       "n_distinct_candidate_families":len({r["family"] for r in check}),
       "n_preselected_independent_floral_organ_QA_taxa":len(chosen),
       "distinct_families_selected":len(used_family),
       "cohort_counts":dict(Counter(r["source_cohort"] for r in chosen)),
       "no_outcome_dependent_replacements":True,
       "selected_by_observed_white_photo_opportunity":True,
       "source_photo_labels_not_independent_flower_validation":True,
       "not_prospective_source_unbiased_global_white_morph_prevalence":True,
       "family_frame_limited_to_existing_two_historical_cohorts":True,
    }

def assess_population_organization(cs,photos,bridge):
    require(len(cs)==EXPECTED_LITERATURE,"strict literature positive-evidence freeze changed")
    require(bridge["totals"]["literature_species"]==EXPECTED_LITERATURE
       and bridge["totals"]["spatial_50km"]==EXPECTED_SPATIAL_OVERLAP,
       "strict literature x photo spatial denominator changed")
    require(bridge["totals"]["spatial_segregation_only_spatial_matches"]==0,
            "S-only spatial species now available; rerun design")
    named={r["species"]:r for r in photos}
    out=[]
    for c in cs:
        name=c["canonical_name"];r=named.get(name)
        if r is None:continue
        white,ncol,n=int(r["n_white"]),int(r["n_nonwhite"]),int(r["n_classifiable"])
        out.append({
           "species":name,"family":c["family"],
           "literature_state":c["organization_state"],
           "photo_classifiable":n,
           "photo_white_labels":white,"photo_nonwhite_labels":ncol,
           "photo_min5_each":white>=5 and ncol>=5,
           "photo_min10_each":white>=10 and ncol>=10,
           "field_colour_genotypes_verified":False,
           "native_population_interbreeding_status_unverified":True,
        })
    require(len(out)==bridge["totals"]["photo_roster"],
            "actual photo x C/S taxonomy crosswalk different from verified original")
    return out,{
      "n_strict_literature_species":len(cs),
      "n_present_in_photo_roster":len(out),
      "n_both_photo_labels_min10":sum(r["photo_min10_each"] for r in out),
      "n_S_only_spatial_photo_overlap":bridge["totals"]["spatial_segregation_only_spatial_matches"],
      "can_test_C_vs_S_from_current_photo_data":False,
      "not_documented_state_is_not_an_absence":True,
    }

def synthetic_tests():
    require(abs(p_ge_one(1,.0021)-.0021)<1e-12,"binomial at least one")
    require(abs(p_ge_two(1,.0021))<1e-12,"no two with one sampled plant")
    require(abs(p_ge_two(2,.5)-.25)<1e-12,"two positives with two draws")
    for p in (.0005,.0021,.0086,.08):
        a=minimum_n(p,1);b=minimum_n(p,2)
        require(p_ge_one(a,p)>=.95 and p_ge_one(a-1,p)<.95,
                "n95 one minimizer incorrect")
        require(p_ge_two(b,p)>=.95 and p_ge_two(b-1,p)<.95,
                "n95 two minimizer incorrect")
    print("FCP_H1H5_SOURCE_PANEL_AND_RARE_WAL_POWER_SYNTHETIC_PASS")

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--family",type=Path)
    parser.add_argument("--photo-roster",type=Path)
    parser.add_argument("--silene-frequencies",type=Path)
    parser.add_argument("--cs",type=Path)
    parser.add_argument("--cs-bridge",type=Path)
    parser.add_argument("--out",type=Path)
    parser.add_argument("--self-test",action="store_true")
    a=parser.parse_args()
    synthetic_tests()
    if a.self_test:return
    require(all((a.family,a.photo_roster,a.silene_frequencies,a.cs,a.cs_bridge,a.out)),
            "source inputs mandatory")
    photo=read_csv(a.photo_roster)
    panel,summary=select_independent_organ_qc_panel(read_csv(a.family),photo)
    source_coverage,coverage=assess_population_organization(read_csv(a.cs),photo,
                             read_json(a.cs_bridge))
    pow=power_scenarios(read_json(a.silene_frequencies))
    a.out.mkdir(parents=True,exist_ok=True)
    record={
       "schema":"fcp_minimal_h1_h5_actual_photo_panel_and_rare_pigment_morph_power_v1",
       "source_type":"verified_4475_species_photo_roster_and_published_Silene_positive_population_year_frequencies",
       "not_an_ecological_selection_or_molecular_expression_result":True,
       "taxa_panel":summary,
       "fcp_literature_population_organization_overlap":coverage,
       "rare_pigment_morph_sampling":pow,
       "rare_WAL_result_interpretation":"A simple unselected 120-plant sampling pilot is unlikely to find two WAL plants at 0.21pct even within a known WAL-positive population-year. Ascertain WAL-positive plants deliberately for functional assays, separate from unbiased population-frequency estimation. These binomial calculations assume independent sampling and may overstate discovery for geographically clustered morphs.",
       "genetic_neutral_SNP_history_present":False,
       "same_individual_petal_leaf_gene_expression_and_fitness_present":False,
       "independent_spectrophotometry_performed":False,
       "field_sampling_power_for_selection_guaranteed":False,
    }
    (a.out/"h1_h5_field_feasibility_actual_source.json").write_text(
        json.dumps(record,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    with (a.out/"h1_independent_organ_qa_candidate_panel_20families.csv").open("w",newline="",encoding="utf-8") as f:
        cols=list(panel[0])
        w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(panel)
    with (a.out/"h4_crosswalk_literature_photos_NO_GENOTYPES.csv").open("w",newline="",encoding="utf-8") as f:
        cols=list(source_coverage[0])
        w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(source_coverage)
    print("FCP_H1H5_REAL_SOURCE_PILOT_FEASIBILITY_COMPLETE",json.dumps({
      "h1_candidate_species":summary["n_preselected_independent_floral_organ_QA_taxa"],
      "h1_distinct_families":summary["distinct_families_selected"],
      "h4_C_S_photo_roster_overlap":coverage["n_present_in_photo_roster"],
      "h4_S_only_spatial_overlap":coverage["n_S_only_spatial_photo_overlap"],
      "wal_median_frequency_percent":next(r["positive_population_year_frequency_percent"] for r in pow if r["name"]=="WAL_median_positive"),
      "wal_median_p_n120_at_least_two":next(r["p_at_least_two_if_independent_binomial"] for r in pow if r["name"]=="WAL_median_positive"),
      "wal_median_n95_at_least_one":next(r["n_to_have_95pct_at_least_one_binomial"] for r in pow if r["name"]=="WAL_median_positive"),
    },sort_keys=True,allow_nan=False))

if __name__=="__main__":
    main()
