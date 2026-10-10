# FCP 42,111 source species: can 342/649 TRUE dated phylogeny tips form local matched congeneric comparisons?

## Current status

- The original 42,111-species photographed flower-colour atlas has 18,457 one-photo four-state classifiable species. A 50–100km locally conditioned original photo-colour-label permutation did not support an independent rainfall signal, unlike an overly broad same-genus×coarse-cell permutation. Both thresholds were explored after the original outcomes were seen.
- Earlier project trees (34-species JBI, 709-tip FCP H3a) cover only 1–20 exact photo species in each of the original 250/500km congeneric samples (1.1%).
- Exact original family/genus metadata was retrieved for all **1,761/1,761** original 500km congeneric source taxa and all nested **872/872** 250km taxa (59 public taxon-ID batches, no synonymous species substitution).
- Source-pinned dated `GBOTB.extended.LCVP` backbone genuinely includes **342/872 (39.2%)** original 250km and **649/1,761 (36.9%)** original 500km photo source taxa as **direct** matching tips. Direct-only induced dated trees were archived at `results/fcp_global42111_original_local_LCVP_backbone_20261010`. The entire original-source ≥50% direct-tip coverage rule fails on both source cohorts, so full-sample tree-supported climate/colour association remains HOLD. Synthetic genus or family grafts are not interpreted as measured phylogenetic distances.

## Pre-model local phylogenetic opportunity audit

Before attempting any source original photographed floral-colour vs rainfall model, the script `scripts/analysis/audit_fcp_global42111_direct_phylo_nearby_photo_pairs_20261010.py` joins the original source 42,111 species public photo positions only by original numerical species identity to the exact DIRECT original LCVP backbone tip ledger; no photo outcomes are opened for local-pair selection.

On BOTH unchanged original 250km and 500km source congeneric samples, examine each source photograph's true geographic coordinates. Require original same-genus, same original 162-cell and ≤50km or ≤100km **complete-link** geographically local groups. Count:
- number of original direct-tip LCVP species with at least one **different original same-genus source species** also a true LCVP tip in the same local source photo group;
- number of separately represented original genera, original genus×cell groups, and true local tipped source pairs;
- distribution of genuine ORIGINAL dated-LCVP phylogenetic path lengths (not family/genus pseudodistance).
- all original source denominators: 872 / 1,761 one-photo congeneric cohorts, 342 / 649 direct species-tip subsets and excluded direct-tip-unmatched taxa.

Fixed exploratory minimum: **≥100 original direct tipped source species**, **≥20 distinct source genera** and **≥30 local multi-species phylogenetic groups** in the same fixed distance cohort. If any fails, status `HOLD_INSUFFICIENT_ORIGINAL_DIRECT_PHYLOGENETIC_LOCAL_COMPARISONS`. If the gate passes, this ONLY authorizes a carefully bounded exploratory *direct-tree-subset* phylogenetic sensitivity, not reinstatement of a full-sample 872/1,761 species rainfall-selection claim.

No new photos, colour labels, genotype data, independent prospective 2,000+730 source taxa, or original 1,499 deep-photo H1/H2 manuscript are touched. The dated LCVP backbone remains a reference phylogenetic hypothesis, and its direct tip subset can be nonrepresentative of the original photo database.

Source-only CI `.github/workflows/fcp-global42111-direct-phylo-photo-pairs-20261010.yml` fails closed on original taxon/photograph identity drift, missing direct tips, or insufficient comparison coverage. Its evidence is written to `results/fcp_global42111_direct_LCVP_local_photo_pairs_20261010/result.json` once verified. **No phylogenetic effect coefficient has yet been computed.**
