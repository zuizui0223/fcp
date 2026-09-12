# H2 evidence ledger — final validation-cohort decision

Date: 2026-09-12 JST

## Final decision

The broad claim that flower-colour polymorphism follows a recurrent arbitrary hue direction is **not supported** after construction controls.

The supported, narrower pattern is an **achromatic–chromatic axis: white versus the equal mean of the eight non-white palette coordinates**.

In the existing high-depth cohorts this axis is a targeted post-audit decomposition, because the candidate was motivated by already opened H2 geometry. It is therefore retrospective evidence and is frozen as a prospective target for future high-depth expansion.

## Evidence sequence

### A. Label-free continuous H2 geometry

H1 admission used the frozen four-state coarse rule, but H2 modes were reconstructed without those labels from the normalized nine-colour flower palette using Hellinger two-means.

Primary 10% gate:

- discovery: N=152, lambda1=0.541412, isotropic p=0.00009999;
- reserve: N=129, lambda1=0.535045, isotropic p=0.00009999;
- reserve mean squared projection on the frozen discovery axis=0.524157, p=0.00009999;
- absolute discovery/reserve leading-axis alignment about 0.986.

Strict 20% gate:

- discovery N=75 and reserve N=65;
- discovery/reserve concentration and frozen-axis transport all passed the isotropic reference at p=0.00009999.

Conclusion at this stage: recurrent geometry exists, but its source was not yet identified.

### B. Coarse-state-preserving construction null

Within each cohort and each frozen coarse morph state, normalized nine-colour palette rows were permuted across the already selected H2 species while every species x coarse-morph row count was preserved. Label-free two-means and Delta vectors were then rebuilt.

Primary 10% gate:

- discovery lambda1=0.541412 versus structured-null median 0.476228, p=0.001;
- reserve frozen-axis projection=0.524157 versus structured-null median 0.482719, p=0.001.

Strict 20% gate:

- discovery concentration p=0.001;
- reserve frozen-axis transport p=0.001.

Conclusion: the recurrent axis cannot be reduced to four-state composition plus the global state-to-palette mapping alone.

### C. White-axis removal and non-white-only audit

A canonical zero-sum white contrast was defined as

`normalize([1, -1/8, -1/8, -1/8, -1/8, -1/8, -1/8, -1/8, -1/8])`.

After projecting this contrast out of every species Delta axis, apparent residual concentration still looked strong against an isotropic reference, but it did **not** exceed the coarse-state-preserving structured null.

Primary 10%:

- residual discovery lambda1=0.388546; structured-null p=1.0;
- reserve projection on frozen residual axis=0.343953; structured-null transport p=1.0.

Strict 20%:

- residual discovery lambda1=0.400773; structured-null p=1.0;
- reserve residual-axis projection=0.344670; structured-null transport p=0.962.

Restricting to species whose top two H1 coarse states both excluded white also failed the discovery structured-null concentration test:

- primary: N=28/19 discovery/reserve, discovery p=1.0;
- strict: N=10/11, discovery p=0.671.

Conclusion: there is no supported recurrent hue direction beyond the white/non-white component under the construction control.

### D. Fixed achromatic–chromatic target

The fixed white-versus-equal-nonwhite statistic was

`W = mean_i (u_i^T q_white)^2`.

Primary 10%:

- discovery N=152, W=0.514625; structured-null median=0.430808; p=0.001;
- reserve N=129, W=0.514586; structured-null median=0.466546; p=0.001.

Strict 20%:

- discovery N=75, W=0.542355; structured-null median=0.443626; p=0.001;
- reserve N=65, W=0.510517; structured-null median=0.469943; p=0.008.

Verdict: `WHITE_AXIS_TARGETED_SUPPORT_PRIMARY_AND_STRICT`.

## Allowed current H2 claim

Across the existing high-depth validation cohorts, within-species flower-colour polymorphism is disproportionately aligned with an achromatic–chromatic (white versus non-white) axis, beyond that expected from the frozen coarse-state composition and global state-to-palette mapping.

This is a descriptive cross-species pattern in the existing cohorts. Because the specific axis was isolated after opening the initial H2 geometry, the existing-cohort result is not an untouched confirmatory test.

## Hard non-claims

Do not claim from H2 alone that:

- a general recurrent non-white hue axis exists;
- the axis represents a specific anthocyanin, carotenoid, or other pigment pathway;
- white morphs are derived, ancestral, losses, or gains;
- pollinators, climate, life form, or phylogeny cause the axis;
- the result estimates global polymorphism prevalence across 42,111 species;
- the 369/363 species are an unbiased global sample.

## Prospective handoff

The achromatic–chromatic `q_white` vector, statistic `W`, H1/H2 admission rules, and coarse-state-preserving structured null are now frozen before any U100 high-depth expansion result is opened.

The exact high-depth opportunity ceiling remains U100=4,730 species. Failure to reach H1/H2 measurement gates is missing/underidentified, not monomorphism.

H3 may proceed on the existing deep cohorts, but ecological or phylogenetic covariates must not be used to redefine H1, H2, `q_white`, or the H2 selected species sets.
