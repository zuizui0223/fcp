# RGFCA Step 8A — 12,064-species photo-first atlas

Date: 2026-09-10 JST
Status: frozen before any Step-8A image pixels are opened.

## Target

Construct an outcome-blind atlas frame of exactly 12,064 plant species from the existing metadata-only global capacity universe. The number 12,064 is a user-selected atlas target; it must not be described as the number of species in the historical literature corpus. The historical v2.2 literature completeness surface contains 12,064 works, not 12,064 species.

## Source universe

Use `data/frozen/global_monte_carlo_capacity_scan_species_audit_v3.csv`, generated from the recovered metadata-only scan of 42,111 discovered species. No flower-colour outcome or image pixel may be used for eligibility or ordering.

## Maximum common raw-photo depth

Before acquiring any new pixels, count metadata-eligible species at the following fixed candidate raw-photo depths, in this exact descending order:

`100, 80, 60, 50, 40, 30, 20, 10` raw photographs per species.

The Step-8A common depth is the largest candidate depth for which at least 12,064 species have `after_observer_cap >= depth`. There is no post-result relaxation outside this fixed grid.

If no candidate depth supports 12,064 species, Step 8A is not evaluable under this design and stops before pixels.

## Species selection

Among species meeting the selected common depth, order species deterministically by SHA256(`20260910-step8a|inat_taxon_id|species`) and take exactly the first 12,064. No literature evidence, flower colour, climate, family, geography, C*/S*, D, or previous effect estimate enters the ordering.

Persist the complete 12,064-species frame and its source-row lineage before pixels.

## Nested checkpoints

The 12,064 frame is primary. For saturation curves use fixed nested prefixes from the same ordering:

- 1,000
- 2,000
- 4,000
- 8,000
- 12,064 species

Within-species photo-depth rarefaction uses only depths not exceeding the selected common depth.

## Existing opened observations

The prior 1,000-species discovery/reserve measurements are outcome-opened and cannot be relabelled as independent validation. Any overlap with the Step-8A frame must be explicitly marked. New species/photo records remain distinct from the opened cohorts.

## Purpose

Estimate how global flower-colour summaries and polymorphism organization stabilize as the species universe expands to 12,064 species. The literature polymorphism registry remains an independent validation/annotation layer and must not select Step-8A species.
