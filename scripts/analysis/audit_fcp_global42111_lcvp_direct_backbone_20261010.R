#!/usr/bin/env Rscript
# Source-exact 1761 FCP local congener taxon matches against LCVP mega-phylogeny.
# Genus-only / family-only placement is NOT a measured species phylogenetic tip.
# Output a direct-backbone-tip induced time-scaled phylogeny when possible,
# but do not fit PHOTO COLOUR vs CLIMATE or infer genetic evolution.
suppressPackageStartupMessages({
  library(ape)
  library(jsonlite)
  library(V.PhyloMaker2)
})
a <- commandArgs(trailingOnly=TRUE)
if(length(a)!=2L) stop("Usage: Rscript ... original_local_species_family_genus.csv outdir")
infile<-a[[1]];outdir<-a[[2]]
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
d<-read.csv(infile,stringsAsFactors=FALSE,check.names=FALSE)
need<-c("inat_taxon_id","original_scientific_name","original_genus",
        "source_genus_cell_id","in_fixed_250km_original_source",
        "current_family","current_genus","status")
if(!all(need %in% names(d))) stop("Original source exact family audit columns missing")
if(nrow(d)!=1761L || anyDuplicated(d$inat_taxon_id)) stop("Original 1761 source taxa changed")
if(sum(as.logical(d$in_fixed_250km_original_source))!=872L) stop("Original nested 872 taxa changed")
data("tips.info.LCVP",package="V.PhyloMaker2")
data("GBOTB.extended.LCVP",package="V.PhyloMaker2")
tips<-tips.info.LCVP
if(!all(c("species","genus","family")%in%names(tips))) stop("LCVP backbone taxonomy schema drifted")
if(length(GBOTB.extended.LCVP$tip.label)<70000) stop("Wrong GBOTB backbone object")
clean<-function(s) tolower(gsub(" +"," ",trimws(gsub("_"," ",as.character(s)))))
source_name<-clean(d$original_scientific_name)
backbone_name<-clean(tips$species)
if(anyDuplicated(backbone_name)) stop("LCVP backbone ambiguous duplicate normalized species")
if(anyDuplicated(source_name)) stop("Original source species names cannot map uniquely to tips")
lookup<-match(source_name,backbone_name)
taxonomy_ok<-d$status=="EXACT_SPECIES_FAMILY_GENUS_READY" &
   !is.na(d$current_family) & nzchar(d$current_family) &
   !is.na(d$current_genus) & nzchar(d$current_genus)
direct<-taxonomy_ok & !is.na(lookup)
original_genus<-clean(d$original_genus)
genus_match<-rep(FALSE,length(d));family_match<-genus_match
genus_match[direct]<-clean(tips$genus[lookup[direct]])==clean(d$current_genus[direct]) &
  clean(d$current_genus[direct])==original_genus[direct]
family_match[direct]<-clean(tips$family[lookup[direct]])==clean(d$current_family[direct])
usable_direct<-direct & genus_match & family_match
# A genus/family membership match is readiness for a SYNTHETIC backbone graft,
# not sufficient evidence for the order/length of congeneric branches.
genera<-split(clean(tips$family),clean(tips$genus))
family_for_genus<-lapply(genera,unique)
genus_family_supported<-vapply(seq_len(nrow(d)),function(i) {
  g<-original_genus[i];f<-clean(d$current_family[i])
  if(!taxonomy_ok[i] || is.null(family_for_genus[[g]])) return(FALSE)
  length(family_for_genus[[g]])==1L && identical(family_for_genus[[g]],f)
},logical(1))
status<-rep("CURRENT_SOURCE_TAXONOMY_UNRESOLVED",nrow(d))
status[taxonomy_ok]<-"FAMILY_VERIFIED_GENUS_ABSENT_FROM_LCVP"
status[taxonomy_ok & genus_family_supported]<-"GENUS_FAMILY_IN_LCVP_SYNTHETIC_NOT_DIRECT"
status[direct & (!genus_match | !family_match)]<-"EXACT_SPECIES_TIP_GENUS_FAMILY_TAXONOMIC_CONFLICT"
status[usable_direct]<-"EXACT_SPECIES_TIP_IN_LCVP_BACKBONE"
result<-d
result$phylogeny_placement_status<-status
result$direct_lcvp_tip_label<-NA_character_
result$direct_lcvp_tip_label[usable_direct]<-as.character(tips$species[lookup[usable_direct]])
result$direct_lcvp_backbone_tip<-usable_direct
result$synthetic_graft_required<-taxonomy_ok & !usable_direct
write.csv(result,file.path(outdir,"original_local_1761_LCVP_direct_tip_ledger.csv"),row.names=FALSE,na="")

summarize<-function(x) {
  yes<-as.logical(x$direct_lcvp_backbone_tip)
  direct_names<-x$original_genus[yes]
  tab<-table(tolower(direct_names))
  loc<-table(x$source_genus_cell_id[yes])
  list(n_original_species=nrow(x),
       n_current_genus_family_taxonomy_verified=sum(x$status=="EXACT_SPECIES_FAMILY_GENUS_READY"),
       n_exact_direct_lcvp_tip_species=sum(yes),
       n_direct_lcvp_genera=length(tab),
       n_direct_lcvp_genera_with_at_least_two_source_species=sum(tab>=2L),
       n_local_genus_cell_groups_with_at_least_two_direct_tips=sum(loc>=2L),
       n_taxa_genus_family_supported_but_synthetic_tip_required=sum(x$phylogeny_placement_status=="GENUS_FAMILY_IN_LCVP_SYNTHETIC_NOT_DIRECT"),
       n_taxa_with_missing_or_conflicted_family_genus=sum(!x$status=="EXACT_SPECIES_FAMILY_GENUS_READY"),
       direct_tip_fraction=mean(yes),
       original_500km_cohort_coverage_gate_pass=FALSE)
}
allstats<-summarize(result)
inner<-summarize(result[as.logical(result$in_fixed_250km_original_source),,drop=FALSE])
# Source-partition up-front suitability: same 300/50%/30 congeneric genera
gate<-function(s) s$n_exact_direct_lcvp_tip_species>=300L &&
  s$direct_tip_fraction>=.5 && s$n_direct_lcvp_genera_with_at_least_two_source_species>=30L
allstats$original_500km_cohort_coverage_gate_pass<-gate(allstats)
inner$original_250km_cohort_coverage_gate_pass<-gate(inner)
for(cap in c("500","250")) {
  selected<-if(cap=="500") result else result[as.logical(result$in_fixed_250km_original_source),,drop=FALSE]
  tipskeep<-unique(selected$direct_lcvp_tip_label[as.logical(selected$direct_lcvp_backbone_tip)])
  tipskeep<-tipskeep[!is.na(tipskeep)]
  if(length(tipskeep)>=2) {
    subset<-ape::keep.tip(GBOTB.extended.LCVP,tipskeep)
    if(!all(tipskeep %in% subset$tip.label)) stop("Induced LCVP tree omitted original direct species")
    if(is.null(subset$edge.length) || any(!is.finite(subset$edge.length))) stop("Direct backbone was not dated")
    ape::write.tree(subset,file=file.path(outdir,paste0("original_local_",cap,"km_direct_LCVP_tree.tre")))
  }
}
json<-list(schema="fcp_global42111_local_congener_LCVP_direct_backbone_coverage_v1",
 date_jst="2026-10-10",
 status="DIRECT_MEGA_TREE_EXACT_SPECIES_TIP_READINESS_NOT_COLOUR_EFFECT",
 original_source_1761_taxa=1761L, original_source_872_nested_taxa=872L,
 backbone="GBOTB.extended.LCVP", source_tree_n_tips=length(GBOTB.extended.LCVP$tip.label),
 500km=allstats, 250km=inner,
 source_taxon_colour_labels_read=FALSE, source_photo_pixels_read=FALSE,
 source_42111_global_population_unmodified=TRUE,
 n_source_photo_labels_remeasured=0L,
 decision=if(gate(allstats)&&gate(inner)) "DIRECT_LCVP_PHYLOGENY_EXPLORATORY_COVERAGE_SUFFICIENT" else "HOLD_DIRECT_LCVP_SUBGENUS_COVERAGE_FOR_FULL_250_500",
 limitations=c("Grafting unsampled congeneric taxa by genus/family does not resolve within-genus evolutionary history",
 "Old iNaturalist current taxonomy can conflict with historical original photo taxon names",
 "This output matches original source species to dated backbone only, no rainfall–photo-colour phylogenetic coefficient"))
writeLines(jsonlite::toJSON(json,auto_unbox=TRUE,pretty=TRUE,na="null"),
 file.path(outdir,"result.json"))
cat(jsonlite::toJSON(json,auto_unbox=TRUE,pretty=TRUE,na="null"),"\n")
