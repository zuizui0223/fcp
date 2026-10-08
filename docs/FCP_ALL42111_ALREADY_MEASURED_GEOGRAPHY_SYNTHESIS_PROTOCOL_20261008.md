# FCP 42,111 species — integrate ALREADY measured global and spatial flower colours
## 2026-10-08 retrospective geographic readout contract

### Corrected status of the broad FCP research programme

Contrary to the earlier assumption that only 1,499 species had ever been measured, an independent historical FCP branch `analysis/rgfca-42111-tiered-measurement-step8c` completed a globe-wide measurement programme on 2026-09-11/12. Its terminal commit is **2b5390ea74f8d196012499d2d35dd477f9795938**, and its source model is ROI-v4 with a four-category visible floral palette. No new image pixels are needed to re-read these original terminal results. The completed programme produced:

1. **42,111 species-wide breadth denominators**, using one original photograph per nominal taxon. 18,457 photo flower colours were classified (43.829%); 23,654 could not be classified. Colour counts among eligible photographed labels: white 8,597, yellow_orange 6,323, red_pink 2,001, blue_purple 1,536. This is *photographed state prevalence conditional on ROI-classifiability*, NOT genetic FCP prevalence or even the modal colour of every species.
2. **85,337 species×equal-area-cell anchor records**, with 58,431 newly processed source images, 26,904 reused breadth terminals and 2 unresolved source transports. 39,075 taxon×cell anchors (45.789%) were classified. 128 of 162 global equal-area cells have a usable recorded taxon-cell anchor. Unlike the breadth layer, a widespread species appears repeatedly; never use all 85,337 records as equal-species global colour prevalence.
3. **13,416 fixed species-specific cross-cell observer-disjoint photo pairs**, with 3,196 pairs classifiable on both sides (23.82% of pairs). 797 differ in coarse flower colour (24.94% of doubly classifiable pairs), 2,399 match; 10,220 fixed pairs have one or both endpoints not classifiable. The unconditional discordance fraction across all pairs is bounded by 797/13,416 (5.94%) to (797+10,220)/13,416 (82.12%) without assumptions about unknown photos. No genotype, mating population or selective agent is measured.

### Present integrated scientific target

A single-source, reproducible global FCP evidence table should retain these **three distinct statistical universes**, while adding the geography of their limitations: regional differences in classifiability, colour distributions in latitude zones and regional conditional two-cell discordance among matched same-species observer-disjoint pairs. The output is a **descriptive synthesis using source-measured colours**, not a fresh prospective hypothesis test or an independent replication of the old 1,499 high-density model.

- Global species-equal colour composition uses the 42,111 one-photo-per-species breadth source ONLY.
- Spatial regional distribution uses the 85,337 taxon-cell anchor source, assigning its original equal-area cell centres (18 longitude sectors × 9 equal-width sin(latitude) rows) to absolute-latitude 0–30°, 30–60° and 60–90°.
- Pair regional assignment uses great-circle midpoint between the original two distinct equal-area cell centroids and retains all 13,416 fixed species regardless of photo colour classifiability. For each region report total matched pairs, both-classified, same, different, missing, classified-conditional observed mismatch fraction, bootstrap uncertainty and no-assumption missingness bounds.
- A region with 100 or fewer classifiable pairs is marked coverage-limited. **No inferential claim that temperate photo colours differ causally from tropical photo colours** follows from these unbalanced, different-species, cell-scale aggregates.
- Once species metadata are joined, source observer IDs and taxon-cell identity must be kept in the auditable CSV, not uploaded photograph pixels or guessed photos.

### Source byte authority

**Use ONLY source files at commit 2b5390ea74f8d196012499d2d35dd477f9795938**, verify the complete original bytes before computing any derived geographic summary:
- `results/rgfca_42111_breadth_measurement_step8f_20260911/rgfca_42111_species_breadth_measured.csv.gz`, SHA256 `38aa42123b4e9b05753020ff1a3b050f4d14dbd3de3194557ead90b75c0cc605`.
- `results/rgfca_42111_taxon_cell_measurement_step8g_20260911/taxon_cell_measured_85337.csv.gz`, SHA256 `7fddcf3449fbcdbaba63d857e4028636e30d2b330a6b39c01c67d42b0470140e`.
- `results/rgfca_42111_crosscell_discordance_step8g_20260911/crosscell_discordance_pairs_13416.csv.gz`, SHA256 `a5d8e677a87fc3b0771413bbb5b2aa40a1b2242acb4156813ae986518ffa7e4f`.

The source design precedes the current October post hoc white-benefit-cost discussions. The October comparison is itself a RETROSPECTIVE geographical description and must never relabel source-colour states, change prior experimental outcomes, or be presented as confirmatory biological fitness.

### Global opportunity and URL availability in parallel

A separate current [PR #134](https://github.com/zuizui0223/fcp/pull/134) includes the full photo-opportunity ledger and all 42,111 unique historical image/observation IDs. Its corrected original URL query recovered licence-valid original photo URLs for 42,058 taxa and left 53 with taxon mismatch, deleted or missing original image, or nonpermitted licensing. That was SOURCE-METADATA ONLY, not a new image outcome. Those URL-validation totals are not the same denominator as the September historical ROI outcomes. Do not silently replace historical image measurements with re-resolved photo URL states.

Original September broad ROI result already exists and should be used now instead of paying to repeat 42,111 identical image classifications. Extra independent photo classifications can be planned separately and must protect the prospective 2,000+730 unexposed species cohorts from outcome leakage.

### Ecological limitation

The geographical distribution of photo colour variants can be structured by lineage composition, population history, pollinators, vegetation/abiotic environments and citizen-science sampling. Nothing in the source directly identifies heritable colours, opposing net reproductive-fitness benefits/costs or balancing selection. One photo per species cannot distinguish naturally monomorphic plants from rare morphs, and the same taxon may appear in multiple cells in the 85,337 geographic atlas. Lack of classifiable flower pixels is not evidence of white or no FCP.

A complete *descriptive* atlas with all 42,111 taxa does not imply that all 42,111 species are phenotypically measured. The status of each layer, including its unclassifiable observations, is essential to trustworthy macroecology.
