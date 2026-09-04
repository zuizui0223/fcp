# Random photo-first flower-colour boundary atlas — terminal inference state

## Purpose

This is a **new prospective experiment**, separate from both the frozen six-species Chapter-1 analysis and PR #21's terminal 200-species / 60,000-photo experiment.

Its preregistered core question was:

> When georeferenced flowering photographs are repeatedly resampled without fixing a focal species set, do the same flower-colour transition edges recur more strongly than expected under a species-conditioned null, and—if so—do those recurrent edges track macroclimate contrast?

The core sequence was frozen as:

`fresh metadata-only candidate pool -> location-blind photo measurement -> repeated photo-first H1 persistence -> hierarchical H2 environmental concordance`

No step rescued, subsetted, relabelled, or reopened the terminal PR #21 records.

## Terminal state in one sentence

The fresh random photo-first atlas is **fully evaluable but does not support a universal recurrent global flower-colour transition geography or a common macroclimate explanation**. Post-H1 exploratory follow-ups likewise do not provide a robust general substitute: a spatial signal appears only when low-information species are given equal species weight, and it disappears in the reliable species frame and under information weighting.

## Fresh candidate pool

The metadata-only pool was frozen before candidate image pixels opened.

- source: iNaturalist API v1;
- root taxon: Angiospermae (`taxon_id=47125`);
- Research Grade;
- flowering annotation (`term_id=12`, `term_value_id=13`);
- exact species-rank observations with photos and coordinates;
- positional accuracy ≤5 km;
- allowed Creative Commons photo licences only;
- one random page per frozen equal-area query cell;
- no fixed or targeted species list;
- no outcome-based replacement or favourable rerun.

The frozen pool contains **20,845 photographs from 8,989 species in 128 occupied H1 cells**. Under the frozen species cap of two photographs per species per cell, premeasurement sampling capacity was **17,813**, exceeding the fixed 10,000-photo H1 replicate requirement.

## Fresh location-blind measurement

All **20,845** candidate rows reached exactly one terminal measurement state and were joined back to frozen metadata only after measurement completion.

The location-blind measurement path was:

`candidate image -> ROI-v4 flower mask -> masked RGB pixels -> fixed CIELAB nearest-palette assignment -> four coarse biological colour groups`

The four biological groups were `white`, `yellow_orange`, `red_pink`, and `blue_purple`. `mixed_uncertain` remained structural measurement missingness and was never treated as a fifth biological colour state.

Terminal measurement counts:

- classifiable four-state flower colour: **10,103**;
- `mixed_uncertain` / non-evaluable for biological colour: **10,742**;
- white: **4,705**;
- yellow-orange: **3,038**;
- red-pink: **1,358**;
- blue-purple: **1,002**.

The measured table also retains the prospective nine continuous flower-palette fractions, which were used only in post-H1 exploratory analyses after the primary H1/H2 outcomes were fixed.

Measurement table SHA-256:

`5d17ba931bbf712def34f3e5f47341c3c5bd436ac754978c6e34978cf8675098`

## Preregistered H1 — recurrent global transition geography

Primary H1 used the frozen design:

- **18 × 9** equal-area longitude × sin(latitude) grid;
- **10,000 photographs per replicate × 200 replicates**;
- cell-first colour-blind sampling;
- maximum **2 photographs per species per cell per replicate**;
- minimum **5 classifiable photographs per cell** for edge evaluation;
- Jensen-Shannon divergence of four-state cell compositions;
- exact top 10% of evaluable edge intensities per replicate;
- opportunity-denominator edge persistence;
- global statistic: opportunity-weighted concentration of edge persistence;
- **999 within-species morph-label permutations** preserving species geography and sampling structure.

Primary result:

- observed persistence concentration: **0.0557863**;
- upper-tail permutation **p = 0.413**;
- decision: `no_support_excess_recurrent_boundary_concentration`.

Therefore the atlas does **not** support the hypothesis that a common set of global flower-colour transition edges recurs more strongly than expected after conditioning on species geography.

One coarse-grid sensitivity (12 × 6, species cap 2) was nominally positive (`p = 0.018`), but the primary and the remaining preregistered grid/cap sensitivities were not. It is retained as scale-sensitive secondary evidence, not as a positive H1 conclusion.

## Predeclared regional H1 audit

The preregistered reporting plan included global, continent-within and biome-within views. The frozen environmental table contained biome and realm but no continent column, so continent-within was left `not_evaluable` rather than being reconstructed after the outcome.

Using the **same observed H1 persistence map and the same 999 H1 null maps**:

- within-biome concentration: **p = 0.971**;
- within-realm supplemental diagnostic: **p = 0.867**.

Thus the global null result is not rescued by restricting the analysis to shared biomes or realms.

## Hierarchical H2 — macroclimate concordance

H2 was frozen before the fresh measurement outcome and was hierarchically subordinate to H1. Its primary macroclimate contrast used the frozen 250-km WorldClim-derived H1-cell summaries for standardized BIO1, BIO4, BIO12 and BIO15.

Primary H2 diagnostic result:

- opportunity-weighted correlation between H1 edge persistence and multivariate climate contrast: **r = 0.15465**;
- upper-tail matched-null **p = 0.134**;
- supported edges: **175**;
- hierarchical decision: `diagnostic_only_h1_not_supported_no_climate_mechanism_claim`.

Single-variable secondary diagnostics were also non-significant after Holm correction; BIO15 was the largest candidate (`r = 0.19044`, raw `p = 0.107`, Holm `p = 0.428`). Within-biome, within-realm, 100-km and 500-km H2 sensitivities were all non-significant.

Accordingly, the random atlas does not support a common macroclimate explanation for recurrent flower-colour transition geography.

## Post-H1 exploratory H3 — species turnover

H3 was introduced only after the preregistered H1/H2 outcomes. It asked whether colour-persistence geography was especially aligned with turnover in species composition.

- observed species-turnover versus colour-persistence correlation: **r = 0.164**;
- matched null mean correlation: **r ≈ 0.193**;
- upper-tail **p = 0.826**.

Species turnover therefore does not provide an additional explanation for the observed colour-persistence surface beyond the geography already represented in the species-conditioned null.

## Post-H1 exploratory H4a — within-species climate-colour divergence

A capacity audit was run before opening a new within-species climate test. The reliable frame was frozen as species with at least 10 classifiable photographs, at least 5 H1 cells and at least 2 coarse morphs: **74 species / 1,557 photographs**.

H4a used continuous soft four-group colour composition derived from the retained palette fractions. For each species it correlated pairwise colour JSD with pairwise multivariate climate distance, then averaged species correlations equally. The null permuted colour-cell labels within species.

Primary H4a:

- evaluable species: **74**;
- mean species Spearman rho: **0.02605**;
- upper-tail **p = 0.157**;
- decision: `no_support_within_species_climate_colour_divergence_do_not_open_h4b`.

The predeclared H4b directional colour/climate-driver decomposition was therefore **not opened**.

Looser three-cell sensitivities were nominally positive (`p = 0.013–0.014`) whereas the five-cell primary and stricter frames were not. That threshold dependence motivated H5 but did not alter H4a's primary conclusion.

## Post-H4a exploratory H5 — spatial-scale moderation

H5 was explicitly frozen after seeing the H4a threshold sensitivity and is therefore exploratory. Rather than selecting a favourable cell-count cutoff, it directly tested whether species-level climate-colour coupling weakened continuously as occupied geographic span increased.

Primary H5:

- frame: **79 species** with at least 10 classifiable photographs, at least 3 H1 cells and at least 2 morphs;
- moderator: log1p(mean pairwise H1-cell centroid great-circle distance);
- observed scale-moderation Spearman rho: **-0.17549**;
- lower-tail permutation **p = 0.091**;
- decision: `no_support_scale_dependent_weakening`.

Max-span (`p = 0.028`) and occupied-cell-count (`p = 0.006`) diagnostics were nominally positive but were prospectively labelled non-primary sensitivities that cannot rescue H5.

## Post-H5 exploratory H6 — species-specific spatial colour structure

H6 asked a different question from H1: even if species do not share the same global boundary location, does colour composition generally diverge with geographic distance **within species**?

The reliable primary frame reused the previously established ≥10-photo, ≥5-cell, ≥2-morph rule.

Primary H6:

- evaluable species: **74**;
- equal-species mean geographic-distance versus colour-JSD Spearman rho: **0.02775**;
- upper-tail permutation **p = 0.113**;
- decision: `no_support_general_species_specific_spatial_structuring`.

A much looser sensitivity (≥5 photographs, ≥3 cells, ≥2 morphs) contained **188 species** and was nominally positive:

- mean rho: **0.10792**;
- cell-label-null **p = 0.001**.

The signal was concentrated among the sparsest three- and four-cell species and therefore triggered H6b rather than being promoted to a biological conclusion.

## H6b robustness diagnostic — photo-count-preserving null

H6b was frozen after H6 and targeted only the nominally positive 188-species sparse sensitivity. It replaced the H6 cell-mean exchangeability null with a stricter photo-level null that:

- permuted individual classifiable photo colour vectors only within species;
- kept every photo's H1-cell assignment fixed;
- preserved the exact number of classifiable photographs in every species × cell;
- rebuilt cell mean colours before recomputing spatial statistics.

Results:

- 188-species equal-species mean: rho **0.10792**, photo-level-null **p = 0.002**;
- pair-count-weighted mean: rho **0.02805**, **p = 0.148**;
- reliable 74-species frame: rho **0.02775**, **p = 0.314**.

Diagnostic decision:

`sparse_equal_species_signal_survives_photo_count_preserving_null_but_is_not_information_weight_robust`

Thus the loose-frame signal is not merely an artefact of the original cell-label null, but it is **strongly dependent on giving low-information species equal weight**. It is hypothesis-generating only and does not overturn the non-significant reliable H6 primary result.

## Canonical evidence hierarchy

The inferential hierarchy must remain explicit:

1. **Preregistered random-atlas H1:** no common recurrent global flower-colour transition geography (`p = 0.413`).
2. **Preregistered hierarchical H2:** no supported common macroclimate concordance (`p = 0.134`), and H1 was already unsupported.
3. **Predeclared regional H1 audit:** no biome/realm rescue.
4. **Post-outcome exploratory H3/H4a/H5/H6:** no robust general replacement mechanism or spatial rule in their primary frames.
5. **H6b diagnostic:** a sparse equal-species-weight spatial signal survives a stricter photo-level null but disappears under information weighting and in the reliable frame.

No exploratory or diagnostic result may be used to rewrite H1/H2 as positive.

## Ecological conclusion

The random photo-first experiment changes the interpretation of the broader programme in a useful way.

The data do **not** support a single global flower-colour transition template repeated across angiosperms, nor a universal macroclimate gradient that locates such a template. The negative result is also not cleanly explained by species-turnover geography. When the question is moved inside species, the fresh multi-species atlas still does not show a robust general climate-colour or geographic-colour relationship in the information-rich frame.

The strongest defensible contrast is therefore:

> **Spatial organization of flower colour can be real within well-measured focal species, but its geography is not globally synchronized across species.**

That interpretation is consistent with the retained six-species held-out Stage A result (`p = 0.0113`) together with the six-species Stage B failure to confirm a universal shared boundary (`p = 0.0906`) and the new random-atlas H1/H2 results.

The new atlas therefore supports **heterogeneous, species- or context-specific spatial organization** as the next prospective target, not continued post hoc searching for one universal global boundary in the current data.

## What is closed on the current dataset

The current dataset is closed for further outcome-driven threshold searching. In particular:

- H1 is not to be redefined around the one positive coarse-grid sensitivity;
- H2 cannot be reopened as a mechanism test after H1 failure;
- H4b remains closed;
- H5 nominal sensitivities do not become primary;
- the H6 188-species sensitivity does not become the canonical species-level result;
- no further cutoff search over photo count, cell count, spatial span, palette grouping, or weighting is justified on this same outcome set.

Any attempt to test the H6/H6b low-information signal as a biological hypothesis should be a **new prospective experiment** with its sampling frame, minimum per-cell replication, weighting rule and null frozen before new outcomes are measured.

## Legacy evidence retained unchanged

The random atlas does not alter these completed results:

- six-species Stage A: within-species continuous colour is spatially organized (`p = 0.0113`);
- six-species Stage B: one universal shared global transition boundary was not confirmed at the primary scale (`p = 0.0906`);
- PR #21 terminal scale-out: 60,000 location-blind measurements completed, but the experiment was `not_evaluable` under its prospectively frozen species/cohort measurement-completeness gate.

Those results motivated the random photo-first persistence question but were not input data for it.
