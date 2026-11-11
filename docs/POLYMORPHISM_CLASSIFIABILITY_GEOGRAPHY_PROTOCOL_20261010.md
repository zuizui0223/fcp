# Geographic photograph-classifiability audit (2026-10-10)

## Research question

The New Phytologist FCP paper finds that nearby conspecific photographed flower colours are more homogeneous than composition-preserving within-species nulls. This audit asks whether **the photographic/classification opportunity itself has a comparable spatial structure**. It tests the original binary source status (globally classifiable vs nonclassifiable flower-colour photograph), not a biological morph.

Primary: does technical status exhibit 50km geographic depletion of pairwise discordance?
Secondary, descriptive only: do species with stronger technical local depletion also tend to show stronger four-colour local depletion?

**Important:** technical classifiability is NOT a validated biological negative control. Failure may reflect flowering stage, developmental state, plant morphology, image exposure, occlusion, phenology and observer behaviour. Its spatial structure alone cannot explain the original colour-specific effect; absence of structure cannot rule out morph-dependent missingness.

## Source firewall

Frozen source photos, no image reclassification, replacement, downloads or new environmental data:

- Discovery and validation: immutable source commit 5142f7951af0dde5364bb047a566d67e8c479e51
- Third: immutable source commit 7e538e5c51c05a7cc47b2fcf53eea92634c8a863
- Three original file SHA256 values embedded in analysis script and workflow
- High-depth eligibility >=40 globally classifiable colour photos per species, expected historical eligible species 369/363/377
- Technical classifiability analysis includes ALL selected original photo positions for these eligible species, including nonclassifiable outcomes
- Companion four-colour analysis retains the original four states only (white; yellow-orange; red-pink; blue-purple) and classifiable photo positions only
- All geographic distances refer to actual original photographed sites; no source equal-area cell centre substitutes as plant location

## Statistical contract

At 50km and minimum 30 local conspecific photograph pairs, species i has technical global photo-status pair discordance T_global_i and local status discordance T_local_i. The diagnostic statistic is mean_i(T_global_i - T_local_i), with one equal vote per species. Four-colour local depletion is computed separately, from classifiable original photos only.

The null keeps each species' exact binary classification-status counts and all original photo positions fixed. There are 199 whole-species label shuffles and 199 calendar-month-restricted within-species shuffles; photos without dates form a distinct missing-month stratum. Bootstrap 95% interval resamples species 999 times. A status-constant or month-nonexchangeable species contributes its exact conditional null value; it is never taken to show biological lack of variation.

The species-level correlation between technical-status depletion and colour depletion is descriptive only. No causal mediation, selection mechanism, loss-of-function allele or percent artefact attributable to technical missingness is identified.

## Hard decisions

1. Four synthetic tests check deliberate technical spatial clustering, all-status-success no-signal, month-only sorting and alternating unclustered technical statuses.
2. Full three-cohort empirical numbers are not a result until a successful SHA-verified source replay and stored output are available.
3. Positive technical clustering motivates further *colour-conditional* photo quality calibration, not downgrading all genuine biological ITV to artefact.
4. A technical null does not exclude geographically correlated **morph-specific** segmentation/detection error; a technical positive does not establish such an error.
5. No source H1/H2, manuscript result or prospective confirmation status changes. The distinct 42,111-species global opportunity frame and unopened 2,000+730 species remain out of scope.

Status: retrospective analysis implementation; no empirical positive/negative claim before a completed source workflow.
