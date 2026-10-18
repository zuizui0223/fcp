#!/usr/bin/env python3
"""Source-authenticated FCP photo overlap with independent pathway evidence.

Species name matching only. It is NOT allele association, proof of local
flower-colour polymorphism, a genotype/phenotype join, or evidence that
published pathway variants occur in the photographed populations.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

PHOTO_SHA="8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6"
N_PHOTOS=447500
N_TAXA=4475
COLOURS=("white","yellow_orange","red_pink","blue_purple")
SOURCES=[
 {
  "system":"Silene littorea",
  "fcp_accepted_name":None,
  "family":"Caryophyllaceae",
  "study":"Casimiro-Soriguer et al. 2016",
  "doi":"10.3389/fpls.2016.00204",
  "molecular_level":"bud-petal RNAseq and HPLC; F3h1 pigmented/white fold changes 42.2/49.0; Myb1a regulatory candidate",
  "evidence_grade":"biochemical_expression_candidate_not_causal_variant",
  "genetic_variant_causally_identified":False,
  "field_morph_coexistence_documented":True,
  "experimental_morph_fitness_estimated":False,
  "reproductive_isolation_demonstrated":False,
  "target_trait":"petal anthocyanin pink-white polymorphism"
 },
 {
  "system":"Ipomoea purpurea",
  "fcp_accepted_name":"Ipomoea purpurea",
  "family":"Convolvulaceae",
  "study":"Zufall & Rausher 2003; Coberly & Rausher 2008",
  "doi":"10.1093/jhered/esg098",
  "extra_doi":"10.1111/j.1558-5646.2008.00355.x",
  "molecular_level":"F3'H coding loss-of-function at P locus drives purple/pink; white A allele tested for selfing transmission vs viability costs",
  "evidence_grade":"causal_structural_allele_and_separate_fitness_locus",
  "genetic_variant_causally_identified":True,
  "field_morph_coexistence_documented":True,
  "experimental_morph_fitness_estimated":True,
  "reproductive_isolation_demonstrated":False,
  "target_trait":"petal purple-pink at P; distinct white allele at A; separate loci"
 },
 {
  "system":"Mimulus guttatus",
  "fcp_accepted_name":"Erythranthe guttata",
  "family":"Phrymaceae",
  "study":"M. guttatus PLA1/MYB cluster mapping 2012; unpigmented natural morph 2018",
  "doi":"10.1086/694853",
  "extra_url":"https://pmc.ncbi.nlm.nih.gov/articles/PMC11061637/",
  "molecular_level":"MYB5 low expression associated with anthocyanin loss; PLA1 region contains MYB regulatory cluster for corolla/calyx/leaf pigmentation patterns",
  "evidence_grade":"genetic_mapping_and_MYB_expression_association_not_fixed_causal_sequence",
  "genetic_variant_causally_identified":False,
  "field_morph_coexistence_documented":True,
  "experimental_morph_fitness_estimated":True,
  "reproductive_isolation_demonstrated":False,
  "target_trait":"anthocyanin nectar-guide spots and vegetative/calyx pigmentation; NOT whole yellow-vs-white flower"
 },
 {
  "system":"Mimulus aurantiacus",
  "fcp_accepted_name":"Diplacus aurantiacus",
  "family":"Phrymaceae",
  "study":"Streisfeld et al. 2013; premating isolation study 2014",
  "doi":"10.1371/journal.pgen.1003385",
  "extra_url":"https://pubmed.ncbi.nlm.nih.gov/25545789/",
  "molecular_level":"MaMyb2 cis-regulatory mutation supported by crossing, QTL, expression, knockdown and allele imbalance",
  "evidence_grade":"causal_regulatory_divergence_with_pollinator_premating_isolation",
  "genetic_variant_causally_identified":True,
  "field_morph_coexistence_documented":False,
  "experimental_morph_fitness_estimated":True,
  "reproductive_isolation_demonstrated":True,
  "target_trait":"red vs yellow pollinator-associated ecotypes; species-wide coarse photo labels do NOT identify them"
 },
 {
  "system":"Silene flower-only versus whole-plant anthocyanin loss literature",
  "fcp_accepted_name":None,
  "family":"multiple_lineages",
  "study":"Del Valle et al. 2019",
  "doi":"10.1186/s12870-019-2082-6",
  "molecular_level":"HPLC-defined petal-only anthocyanin loss PAL vs systemic WAL; descriptive natural-frequency tables of 13 PAL and 13 WAL reported systems",
  "evidence_grade":"descriptive_literature_biochemistry_not_unbiased_causal_comparison",
  "genetic_variant_causally_identified":False,
  "field_morph_coexistence_documented":True,
  "experimental_morph_fitness_estimated":False,
  "reproductive_isolation_demonstrated":False,
  "target_trait":"tissue restriction of anthocyanin loss, not presence of visually white photos alone"
 },
]


def verify(condition,msg):
    if not condition:
        raise ValueError("FCP_PATHWAY_BRIDGE_INVALID: "+msg)


def checksum(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as fd:
        for ch in iter(lambda:fd.read(1<<20),b""):
            h.update(ch)
    return h.hexdigest()


def make_registry(photo):
    verify(len(photo)==N_PHOTOS and photo.inat_taxon_id.nunique()==N_TAXA
           and photo.photo_id.nunique()==N_PHOTOS,"full original photo ID/table failed")
    cv=photo.loc[photo.global_classifiable.fillna(False).astype(str).str.lower().isin(
        ("true","yes","1"))&photo.morph.isin(COLOURS)]
    verify(len(cv)==227107,"original coarse-colour classifier denominator drift")
    species_summary={}
    for name,g in photo.groupby("species",sort=False):
        q=cv.loc[cv.species==name]
        species_summary[name]={
            "n_real_original_photos":len(g),
            "n_fourstate_classified":len(q),
            "source_cohort":str(g.fcp_source_cohort.iloc[0]),
            "photographed_coarse_state_counts":{
                morph:int((q.morph==morph).sum()) for morph in COLOURS
            },
        }
    rows=[]
    for source in SOURCES:
        item={k:v for k,v in source.items()}
        match=item["fcp_accepted_name"]
        item["exists_in_fivecohort_measured_photo_frame"]=bool(match in species_summary)
        item["photo_phenotype"]=species_summary.get(match)
        item["external_molecular_genotype_join_established"]=False
        item["original_field_population_same_as_study"]=False
        item["coarse_photo_state_proves_original_molecular_trait"]=False
        item["color_to_gene_cross_species_validation_established"]=False
        rows.append(item)
    verify(len([x for x in rows if x["exists_in_fivecohort_measured_photo_frame"]])==3,
           "real molecular-study taxon/photo alias overlap changed")
    verify(next(x for x in rows if x["system"]=="Silene littorea")["photo_phenotype"] is None,
           "single-species Silene molecular anchor incorrectly joined to photos")
    tgt=next(x for x in rows if x["system"]=="Mimulus guttatus")
    verify(tgt["photo_phenotype"]["photographed_coarse_state_counts"]["yellow_orange"]==60
           and tgt["photo_phenotype"]["n_fourstate_classified"]==60,
           "nectar-guide spots are not encoded by whole-flower coarse-state FCP photo labels")
    im=next(x for x in rows if x["system"]=="Ipomoea purpurea")
    verify(im["photo_phenotype"]["n_fourstate_classified"]==47,
           "Ipomoea matching photo phenotype support changed")
    return {
        "schema":"fcp_fivecohort_external_molecular_pathway_photo_overlap_audit_v1",
        "study_count_literature_items":len(SOURCES),
        "unique_direct_species_photo_overlaps":3,
        "original_real_photo_rows":447500,
        "original_species":4475,
        "original_classifiable_fourstate_photos":227107,
        "same_provider_photo_and_external_molecular_study_are_independent_sources_but_not_joined_at_individual_genotype_level":True,
        "all_pathway_interpretations_must_be_trait_organ_specific":True,
        "do_not_pool_source_locus_pvalues_across_studies":True,
        "not_a_systematic_molecular_meta_analysis":True,
        "molecular_records":rows,
        "conclusion":"Same conserved anthocyanin pathway family is repeatedly implicated, but distinct structural/regulatory loci and pleiotropy/selection/assortative mating outcomes are nonexchangeable; the actual four-colour photographed states cannot identify causal variant or morph maintenance."
    }


def main():
    a=argparse.ArgumentParser()
    a.add_argument("--photos",required=True,type=Path)
    a.add_argument("--out",required=True,type=Path)
    args=a.parse_args()
    verify(args.photos.is_file() and checksum(args.photos)==PHOTO_SHA,
           "original 447500 measured photo archive SHA256 mismatch")
    df=pd.read_csv(args.photos,usecols=[
        "species","inat_taxon_id","photo_id","fcp_source_cohort",
        "global_classifiable","morph"],low_memory=False)
    result=make_registry(df)
    result["source_original_photo_CSV_sha256"]=PHOTO_SHA
    args.out.mkdir(parents=True,exist_ok=True)
    path=args.out/"fcp_molecular_pathway_external_study_photo_bridge.json"
    path.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")
    print("FCP_EXTERNAL_MOLECULAR_PATHWAY_FIVECOHORT_BRIDGE_AUTHENTICATED",
          json.dumps({x["system"]:x["photo_phenotype"] for x in result["molecular_records"]
                      if x["exists_in_fivecohort_measured_photo_frame"]},sort_keys=True))


if __name__=="__main__":
    main()
