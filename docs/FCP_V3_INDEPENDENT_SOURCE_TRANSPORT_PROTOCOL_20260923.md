# FCP v3 — independent measurement-system transport protocol

Date: 2026-09-23 JST  
Status: **PROSPECTIVE SOURCE-QUALIFICATION PROTOCOL; NO SOURCE SELECTED**

## Purpose

FCP v3 asks whether the same flower-colour estimands transport across an image-generation and measurement system that is independent of the iNaturalist pipeline used in FCP v1/v2.

The target is not “another set of species”. The target is:

> **same estimands × different measurement system**

Species overlap is therefore useful and is not treated as leakage.

## Admission prerequisite

v3 cannot open until FCP v2 has a terminal measurement-validity report.

A v2 failure does not prohibit v3, but it changes the purpose: v3 would then test whether a better-controlled measurement system removes the v2 technical sensitivity rather than merely replicating the same measurement pipeline.

## What counts as an independent source

Preferred order:

1. **paired calibrated field photography** of the same flower/individual under:
   - locked/manual exposure with a colour reference;
   - ordinary auto-exposure JPEG;
2. a curated external field-image dataset with documented acquisition, exposure and colour calibration;
3. another community-science source only when its capture/storage/processing pipeline is demonstrably independent of iNaturalist.

The following do **not** qualify:

- another iNaturalist species subset;
- a new iNaturalist date window;
- iNaturalist images redistributed through GBIF or another aggregator;
- a source that has inherited the same image bytes or transformations.

## Minimum source metadata

Before any colour outcome is opened, the independent source must document:

- source organization/provider;
- acquisition device class;
- image file format;
- whether RAW is available;
- exposure mode;
- white-balance mode;
- colour-reference availability;
- compression/resizing pipeline;
- georeferencing precision when spatial analyses are planned;
- licensing/reuse status;
- unique image/individual identifiers;
- whether multiple images come from the same flower/individual.

## Strongest preferred design: paired same-flower measurement

For every sampled flower/individual:

1. calibrated reference image:
   - colour reference included in frame;
   - locked/manual exposure;
   - locked white balance where possible;
   - RAW retained when available;
2. ordinary auto-exposure image:
   - same flower;
   - ordinary field framing;
   - JPEG processing allowed.

This produces a direct within-flower measurement-system contrast. Biological colour is fixed by design.

## Frozen FCP estimands carried into v3

No axis or species-level statistic is rediscovered from the new source.

### D

Four-state diversity:

`D = 1 - sum_k p_k^2`

with the same biological coarse-state definitions.

Report cross-system:

- species-level D difference;
- Spearman rank correlation;
- Lin concordance;
- calibration slope/intercept;
- error as a function of exposure/background/ROI diagnostics.

### q_white / W

Use the existing fixed q_white contrast without rotation or refit.

Report:

- W in the calibrated/reference system;
- W in the ordinary image system;
- per-species direction cosine agreement;
- excess alignment relative to the same construction-preserving null architecture.

v3 does not use the independent source to search for a new best axis.

### Spatial organization

Only if the independent source has adequate georeferencing:

- species-specific geographic colour organization;
- matched flower-minus-background response;
- cross-system agreement of the species-level spatial statistic.

If georeferencing is insufficient, spatial transport is marked unavailable rather than approximated.

## Primary measurement-system contrasts

The first v3 claims concern measurement transport, not ecology.

### V3-M1: image-level agreement

For paired images of the same flower:

- coarse-state agreement;
- continuous palette distance;
- Lab delta-E;
- classifiability agreement.

### V3-M2: species-level D agreement

For species with adequate paired individuals:

- D reference versus D ordinary;
- rank/concordance;
- bias and heteroscedasticity.

### V3-M3: H2 geometry transport

With q_white fixed:

- per-species delta-vector cosine;
- W difference between systems;
- structured-null comparison in each system.

### V3-M4: technical mediation

Ask whether cross-system disagreement is explained by technical channels measured before biological opening:

- clipping;
- luminance;
- background contrast;
- ROI stability;
- compression/resolution.

This is a measurement model, not an adaptive ecological model.

## Outcome firewall

Source qualification, cohort construction and all technical diagnostics are frozen before opening:

- morph;
- D;
- q_white projection;
- W;
- spatial organization.

The independent-source measurement code may know that images are paired by individual, but must not receive the previous FCP species-level outcomes during technical qualification.

## No-rescue rules

After independent-source outcomes are opened:

- no q_white refit;
- no new source substitution;
- no species replacement because transport is poor;
- no new exposure threshold;
- no new colour-reference correction method;
- no switch between RAW/JPEG branches based on H2 success;
- no dropping species because D disagreement is large;
- no ecological predictor is added to explain a measurement-system failure.

## Interpretation states

v3 should not collapse all evidence into one pass/fail label. It reports at least:

- **image transport**;
- **species-level D transport**;
- **H2 geometry transport**;
- **spatial transport**, if available.

This prevents one successful component from masking another failed component.

## Relationship to the current paper

v3 is not required to finish the current New Phytologist submission.

Its role is to answer the strongest unresolved external-validity question left by FCP v1/v2:

> Does the apparent biological structure persist when the image-generation process itself changes?
