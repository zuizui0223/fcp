# RGFCA: independent reserve-species replication

The discovery atlas now has a positive but small exploratory within-species
distance–colour association: **369 species, 21,424 photos, mean rho = 0.0270213,
999-permutation p = 0.001**. This motivates a new prospective replication; it
does not rescue the non-supported common-boundary results.

## Cohort fixed before images

Use the entire complementary 500 taxa from the original 1,000-taxon metadata
pool, all 100 frozen candidate photographs per taxon. The original hash-selected
500 taxa reproduce the discovery measurement IDs exactly. No favourable species
selection, remeasurement of discovery, replacement of failed images or extra
acquisition after outcomes is allowed. This is a **new 50,000-photo measurement**,
not 50,000 already analysed photos relabelled as validation.

The [metadata audit](supporting/rgfca_reserve_metadata_audit_v1.json) checks exact
Git-object inputs and verifies zero photo/observation overlap with discovery,
the older photo-first measurements, H9 fresh metadata and its exclusion ledger.
All stated positional accuracies are at most 4,998 m; observer caps are two photos
per species (median 95 distinct observers per species). Calendar-quarter coverage
has median four quarters, but one-quarter species also exist.

At 0.25 degrees, median occupied cells = 96, minimum = 7; at 0.5 degrees, median
= 78, minimum = 4. Some species have no pairs within 100 or 250 km. The full
[500-species geometry table](../data/frozen/rgfca_reserve_geometry_audit_v1.csv)
reports 100/250/500-km support. These are metadata diagnostics, not evidence of
uniform global sampling or a colour-selected resolution. The present replication
uses the same all-pairs estimand as discovery, not a new boundary-scale search.

## Fixed measurement and inference

1. Reuse exactly the qualified ROI-v4 measurement infrastructure and frozen
   four-component palette probabilities, with separate reserve-only blind IDs.
   Workers receive only a blind ID, filename and licence. Images remain ephemeral.
2. Require all 256 terminal partitions before joining colour to coordinates.
   Retain acquisition failures and unclassifiable records; admit all species
   with at least 40 classifiable photos, requiring at least 250 eligible species.
3. Primary: equal-species mean Spearman correlation of all-pairs geographic
   distance and colour JSD, with 999 within-species complete-vector permutations.
   A positive mean and upper-tail p < 0.05 indicate directional replication.
4. Required robustness: different-observer pairs, and within-calendar-quarter
   colour permutations. Robust replication requires all three tests to pass;
   neither sensitivity can replace a failed primary. Full details and fixed
   seeds are in the [contract](supporting/rgfca_reserve_replication_contract_v1.json).
5. Before any flower-specific interpretation, require a matched background
   control. Count the same twelve palette anchors in both flower and the frozen
   surrounding annulus during each image's first decode. Test geographic
   distance against flower JSD minus background JSD using joint same-photo
   permutations. This is specified before reserve pixels and cannot be replaced
   by an easier control after results. A background correlation may reflect real
   ecology as well as photographic conditions; neither comparison is causal.

The first CI gate is metadata and synthetic implementation validation only.
Pixel measurement requires a later hash-bound authorization referencing the
successful preflight. No previous authorization is edited or reused. Preserve
failed gates and every negative or non-evaluable result; no selective reruns.

## Interpretation and publication boundary

This validates a **photo-derived within-species spatial association across an
admitted iNaturalist frame**, not all plants, genetic colour polymorphism,
species prevalence, a shared global boundary, pollinator causation or an
environmental mechanism. New species and photo IDs share the platform and
measurement model, so systematic colour/illumination, identification, observer,
habitat and phenology biases can recur. Quarter stratification is only an
observation-season proxy. A replicated small effect warrants ecological
interpretation with these limits, not a claim of a universal cause.

The 6-species and 34-species studies remain immutable legacy. Publication still
requires completed replication, measurement/bias assessment, species-free atlas
figures, species-conditioned inference, an explicit full exploration record,
manuscript and reproducibility bundle. Statistical significance alone does not
complete the research goal.
