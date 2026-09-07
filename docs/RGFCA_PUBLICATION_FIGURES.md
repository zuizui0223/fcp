# RGFCA publication figures: discovery evidence

Status: exploratory discovery figure set, not a submission-ready package.
These figures read only the completed discovery result at immutable commit
`29584f3ad7ae0cd99a1d8f43459252af38f615da`. They never read reserve measurements,
legacy six-species outcomes or the literature-derived 34-species data.

## Chart contracts

| Figure | Question and supported takeaway | Form, grain and full denominator |
|---|---|---|
| 1: Discovery photo-derived colour map | Where are the evaluated observations, and what colour mixtures were measured? Coverage is geographically uneven; the display is descriptive, not a pooled-species test. | Static geographic scatter; one point per eligible photo; all 21,424 photos from 369 species; no taxon legend, fitted field or boundary overlay. |
| 2: Within-species distance-colour association | How large and heterogeneous is the exploratory association, and where does its species-equal mean fall relative to random labelling? The mean is small and positive; independent replication remains pending. | Two histograms: all 369 species-level Spearman coefficients and all 999 synchronized global null means. Observed mean and zero marked explicitly; full species-effect range retained. |

Figure 1 uses longitude/latitude in degrees, not an equal-area projection or an
inference grid. No spatial smoothing or geographic exclusion is allowed.
Natural Earth low-resolution land geometry bundled with GeoPandas 0.14.4 is
display context only; its source files are hashed. Points use a deterministic
photo-ID hash order. Point overlap is not an abundance or occupancy estimate.

Figure 1's explicit palette is an intentional four-category semantic exception:
white `#F4F1E8`, yellow/orange `#DBAD34`, red/pink `#B84A78`, blue/purple `#655AB0`.
Display RGB is the convex mixture of these anchors weighted by each frozen
four-component vector. This is not the original photographic RGB, calibrated
reflectance or pollinator-perceived colour. A four-anchor key and direct labels
explain the encoding. A map need not be fully interpretable in grayscale; all
inferential distinctions in Figure 2 also use line style and direct labels.
Figure 2 uses one blue root `#477A9F` plus charcoal/gray neutrals, with white
backgrounds, dark keylines, zero-based counts and no significance colouring of
individual species. These are third-party FCP publication figures: no OpenAI
branding is used.

The final surface is the scientific manuscript's static PNG/PDF figure pair,
not a dashboard. Figure 1 is 11 by 5.6 inches; Figure 2 is 11 by 4.8 inches.
Exports are 240-dpi PNG and PDF in `docs/figures/`; the source/output SHA-256
manifest is `docs/supporting/rgfca_publication_figure_manifest_v1.json`.
QA includes independent reconstruction of the pooled statistic and p-value,
full cohort/ID/vector validation, two-render byte equality in one pinned runtime,
and visual inspection of the final PNGs. Cross-platform PDF/PNG byte equality is
not assumed; input bytes and numerical results must agree exactly/tightly.

## Figure legends

**Figure 1. Geographic distribution of discovery colour measurements.**
Each point represents one of 21,424 classifiable photographs in the 369-species
inferential frame. Point colour is a display mixture of four palette anchors,
not calibrated reflectance. Species names are suppressed only in the display;
all statistical comparisons retain species identity. No inference can be read
from a region's visual point density or apparent pooled-colour turnover.
Natural Earth outlines are geographic reference, not a tested boundary layer.

**Figure 2. Exploratory within-species distance-colour association.**
(a) Distribution of species-specific rank correlations between all within-species
geographic pairwise distances and Jensen-Shannon colour dissimilarities. Each
species has equal inferential weight; pairs are not independent observations.
(b) Distribution of 999 species-equal mean statistics after reassigning complete
colour vectors among the fixed coordinates within each species. The observed
mean is 0.0270213; the upper-tail Monte Carlo permutation p-value, including the
observed configuration, is 0.001. The shaded null interval is not a confidence
interval for the observed mean. This postoutcome exploratory candidate is not
independent replication, evidence of a shared boundary, or a flower-specific
ecological mechanism; observer, season and matched-background checks are required.

## Explicit omissions

A real flower-crop photo bar is not manufactured from palette swatches. The
discovery result stores measurements rather than a provenance-verified collection
of redistributable ROI crops. Until licensed crops and their measurement IDs are
available, Figure 1 is a colour-measurement map only. Independent reserve outcomes
and any matched-background results are intentionally absent from this version.
Their future figures must follow the frozen, complete-cohort analysis gates.
