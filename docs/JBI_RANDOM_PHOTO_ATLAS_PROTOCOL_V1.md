# JBI random photo atlas v1 — prospective ecological validation protocol

Status: **frozen before opening any image pixel from the new source window** (2026-09-03).

This analysis is a new exploratory/validation experiment. It is **not** a rescue, continuation or reinterpretation of the frozen terminal 200-species / 60,000-photo atlas. No terminal selected photograph, terminal colour measurement, terminal coordinate join or terminal downstream result may enter this experiment.

## Ecological questions

The experiment asks three ordered questions.

1. **Recurrent transition zones:** when flowering photographs are sampled without fixing species in advance, do the same broad geographic flower-colour transition zones recur across independent random replicates?
2. **Climate concordance:** are recurrent flower-colour transitions stronger across adjacent map cells that differ more strongly in broad climate?
3. **Species-conditioned persistence:** do recurrent transition zones remain stronger than expected when morph labels are randomized within species while preserving the exact sampled geography and species composition?

A result may be `supported`, `unsupported` or `not_evaluable`. `not_evaluable` never becomes `unsupported`.

## Independent source window

Primary source: iNaturalist API observations satisfying all of the following:

- created from **2026-08-28 00:00:00 UTC through 2026-09-02 23:59:59 UTC**;
- Angiospermae (`taxon_id = 47125`);
- research grade;
- photographs present;
- public georeference present;
- not captive/cultivated;
- Plant Phenology controlled term (`term_id = 12`) equal to Flowers (`term_value_id = 13`);
- non-obscured public coordinates;
- observation taxon resolved at species or infraspecific rank;
- one deterministic first photograph per observation.

The frozen terminal atlas used the iNaturalist Open Data snapshot dated 2026-08-27. The new source window begins after that snapshot and is therefore source-time independent of the terminal 60,000-photo selection. A direct observation/photo-ID overlap guard is additionally applied whenever terminal IDs are available to the workflow.

Metadata are retrieved completely for the closed source window before image pixels are opened. A hard metadata safety ceiling of 100,000 observations is fixed; reaching that ceiling makes the source universe `not_evaluable` rather than permitting truncation.

## Preimage reservoir

Before any new image pixel is opened, exactly **12,000 observations** are selected from the eligible metadata universe. If fewer than 12,000 eligible observations exist, the experiment is `not_evaluable` and no images are measured.

Selection is independent of colour and uses only latitude, longitude and deterministic hashes:

- primary grid: 24 longitude bins × 12 equal-area `sin(latitude)` bins;
- cells are traversed in deterministic hash order;
- observations within cells are traversed in deterministic hash order;
- round-robin selection across occupied cells continues until exactly 12,000 records are frozen;
- seed/salt: `20260903` / `jbi-random-photo-atlas-v1`.

Taxon identity is retained only for dominance control and species-conditioned nulls. No species list is fixed in advance and no species is required for inclusion.

## Frozen image measurement

The 12,000 frozen reservoir photographs are measured with the already-qualified **ROI-v4** flower estimator, without changing its detector weight, EfficientSAM revision, ROI contract or colour-summary definitions.

Frozen estimator identity:

- model: `jbi-atlas-roi-estimator-v4`;
- trained detector SHA-256: `f1aaeec4664fe2c178e5cf2bc1f508977bef3e4aa7b40613026cb8ae3de789d5`;
- EfficientSAM revision: `d525f622e6f640acf5a0fc37c7ca1f243da5bde0`;
- flower summaries: CIELAB L/a/b mean, SD and q10/q50/q90 as already emitted by ROI-v4.

Images are downloaded only ephemerally by measurement workers and are not uploaded as workflow artifacts.

A photograph is morph-evaluable only when ROI-v4 returns `automated_colour_state_admitted` and finite flower L/a/b q50 values.

## Coarse morph classification

Coarse morphs are assigned deterministically from the ROI-v4 flower-median CIELAB values, converted to sRGB/HSV with the frozen transform used by `skimage.color.lab2rgb`. No map or ecological result is inspected before these rules are fixed.

Let HSV hue be in degrees, saturation `S` and value `V` in [0,1].

- `white`: `S < 0.20` and `V >= 0.72`;
- `yellow`: `40 <= hue < 80` and `S >= 0.20`;
- `orange_red`: `hue < 40` and `S >= 0.35`;
- `pink_purple`: (`hue >= 300` or (`hue < 40` and `0.20 <= S < 0.35`)) and `V >= 0.45`;
- `blue`: `180 <= hue < 300` and `S >= 0.20`;
- `other_uncertain`: all remaining admitted measurements.

`other_uncertain` is retained for QC but excluded from the five-morph composition denominator. No photograph is reassigned after map inspection.

Primary morph composition uses exactly five states: `white`, `yellow`, `orange_red`, `pink_purple`, `blue`.

## Random replicate sampling

Number of replicates: **200**.

Each replicate attempts to select exactly **4,000 five-morph-evaluable photographs** from the frozen 12,000-photo reservoir using deterministic pseudorandom priorities derived from seed `20260903 + replicate_index`.

Hard contribution controls applied in every replicate:

- maximum 25 selected photographs per spatial cell;
- maximum 2 selected photographs from the same species key within a spatial cell;
- no image appears more than once in a replicate.

If a replicate cannot reach exactly 4,000 photographs under those frozen caps, the whole experiment is `not_evaluable`; the denominator or caps are not relaxed.

Replicate-level species dominance is audited. The maximum global contribution of any one species is reported, but species are not otherwise balanced globally.

## Spatial support and transition intensity

The same frozen 24 × 12 equal-area grid is used for inference.

A cell is evaluable in a replicate when it contains at least **8 selected five-morph photographs**.

For each evaluable cell, the five-state morph-composition vector is calculated with no pseudocount.

Adjacency is rook adjacency:

- north/south adjacent `sin(latitude)` rows;
- east/west adjacent longitude columns;
- longitude wraps at ±180°.

An edge is evaluable in a replicate only when both cells are evaluable.

Primary transition intensity is base-2 **Jensen–Shannon divergence** between the two five-morph composition vectors. Dominant-morph switching is a secondary descriptive indicator.

Within each replicate, strong transitions are the top **10%** of JSD values among evaluable edges, using the empirical 90th percentile and `JSD >= threshold`.

## Boundary persistence

For each grid edge:

- `evaluable_count` = number of replicates in which both cells were evaluable;
- `strong_count` = number of those replicates in which the edge was a strong transition;
- `persistence = strong_count / evaluable_count`.

An edge is confirmably evaluated only when `evaluable_count >= 100` of 200 replicates.

A recurrent edge is a confirmably evaluated edge with `persistence >= 0.60`.

Connected recurrent edges are grouped by shared grid-cell endpoints. The primary recurrent-zone statistic is the size of the largest connected recurrent-edge component.

## Species-conditioned null

The species-unfixed atlas is not allowed to treat species turnover itself as flower-colour convergence.

For each of **199** null ensembles, using seed `20260903 + 100000 + null_index`:

- coordinates, observation identities, species identities and the 200 replicate membership sets are fixed;
- five-state morph labels are permuted **within species** across the frozen morph-evaluable reservoir;
- the full cell-composition, JSD, replicate-specific 90th-percentile strong-transition, persistence and connected-component calculation is repeated.

No morph is moved between species.

The null statistic is the largest recurrent connected-component size at the same `evaluable_count >= 100` and `persistence >= 0.60` thresholds.

### Decision: recurrent common transition zones

- `supported`: observed largest recurrent component contains at least 3 edges **and** exceeds the 97.5th percentile of the 199 species-conditioned null component sizes;
- `unsupported`: all evaluability gates pass but the conjunction is not met;
- `not_evaluable`: source, measurement, replicate or edge-evaluability gates fail.

### Decision: persistence after species conditioning

- `supported`: the same observed-vs-within-species-null conjunction above passes;
- `unsupported`: full null is evaluable but it does not pass;
- `not_evaluable`: the null cannot be completed under the frozen design.

## Climate concordance

Climate is attached only after the morph-transition surface is frozen.

Primary static environmental variables are WorldClim 2.1 at 10 arc-minutes:

- BIO1 annual mean temperature;
- BIO4 temperature seasonality;
- BIO12 annual precipitation;
- BIO15 precipitation seasonality.

For each grid cell, climate is sampled at the equal-area cell centre. For each confirmably evaluated edge, the absolute difference between the two cell-centre values is calculated for each variable.

For each climate variable separately, a binomial GLM is fit to `strong_count / evaluable_count` using:

- standardized absolute climate difference as the focal predictor;
- standardized great-circle distance between cell centres as the fixed covariate.

The four focal climate p-values are Holm-adjusted. Report the odds ratio per 1-SD increase in absolute climate difference and its 95% confidence interval.

### Decision: climate-concordant recurrent boundaries

- `supported`: at least one of the four frozen climate predictors has OR > 1, Holm-adjusted p <= 0.05 and 95% CI lower bound > 1;
- `unsupported`: climate data and all four models are evaluable but none meets the full conjunction;
- `not_evaluable`: transition surface, climate extraction or any required model family is incomplete.

No variable may be substituted or promoted outside this four-variable family after results are observed.

## Geographic sensitivities

Global inference is primary. Continent-within and biome-within summaries may be reported as descriptive sensitivities only if they can be computed without changing the global sample, grid, morph rules or transition thresholds. They cannot replace a global decision.

## Claim boundary

Even a supported result establishes recurrent photographic flower-colour turnover under this sampling design and concordance with broad environmental gradients. It does not establish a universal evolutionary boundary, morph-specific fitness, causal climate selection, pollinator mechanism or taxonomic convergence.

If recurrent transitions disappear under the within-species label null, the appropriate conclusion is that apparent global boundaries are substantially explained by species composition rather than a cross-species flower-colour transition field.

## No post-result substitutions

The following are frozen before new pixels are opened and cannot be changed in response to the result:

- source window and iNaturalist filters;
- 12,000-photo reservoir denominator;
- 24 × 12 equal-area grid;
- seed/salt;
- ROI-v4 estimator identity;
- morph definitions;
- 200 replicates and 4,000-photo denominator;
- cell and species contribution caps;
- minimum cell support = 8;
- JSD transition metric;
- top-10% strong-transition rule;
- edge confirmability = 100 replicates;
- persistence threshold = 0.60;
- recurrent-zone minimum = 3 connected edges;
- 199 within-species null ensembles;
- 97.5th-percentile null decision threshold;
- four WorldClim predictors;
- climate GLM and Holm family;
- supported / unsupported / not_evaluable rules.
