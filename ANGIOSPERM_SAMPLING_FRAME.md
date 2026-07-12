# Angiosperm-wide sampling frame

## Goal

Test why floral-colour polymorphism is observed in some angiosperm species but not others without conditioning the dataset only on already famous polymorphic taxa.

The empirical design therefore has two linked samples.

1. **Positive-evidence sample**: species with literature evidence for floral-colour variation, coded for spatial organization of morphs.
2. **Background sample**: species sampled from the angiosperm phylogeny independently of whether polymorphism has been reported.

The second sample is essential. A literature-only collection of polymorphic species cannot estimate the probability of polymorphism across pollination systems.

## Primary outcome hierarchy

`outcome_class` is coded as:

- `monomorphic_confirmed`: adequate evidence for one floral-colour state across the sampled scope.
- `geographic_mosaic`: two or more colour morphs occur across the species, but sampled populations are predominantly locally monomorphic.
- `mixed`: both locally mixed and locally differentiated/near-monomorphic populations occur.
- `local_coexistence`: two or more colour morphs coexist within the same population and flowering period.
- `unknown`: insufficient evidence. This is not treated as monomorphic.

A separate binary variable `species_level_colour_variation` records whether more than one colour state is documented anywhere in the species.

## Pollination hierarchy

Do not collapse the first analysis directly to bee/bird/wind. Code three levels.

### Level A: pollen vector

- `abiotic`
- `biotic`
- `mixed`
- `unknown`

### Level B: dominant mode

- `wind`
- `water`
- `insect`
- `bird`
- `bat`
- `other_vertebrate`
- `ambophilous`
- `mixed_biotic`
- `unknown`

### Level C: functional pollinator guilds

Multi-label field, for example:

- bee
- butterfly
- moth
- fly
- beetle
- wasp
- bird
- bat
- other

`pollinator_guild_count` and `pollination_specialization` are derived separately.

## Floral display opportunity

Abiotic and biotic species are not directly comparable unless the opportunity for visible floral-colour variation is coded.

Required variables:

- `showy_perianth`: yes/no/unknown
- `display_structure`: petal_or_tepal / bract / stamen / inflorescence / reduced / other
- `human_visible_colour_display`: yes/no/unknown

Primary tests of pollination mode should be repeated in:

1. all angiosperms;
2. only species with a human-visible floral display;
3. only species with homologous display structures (for example petal/tepal display).

This separates a true difference in colour-polymorphism dynamics from the simpler fact that many wind-pollinated flowers have reduced, inconspicuous floral displays.

## Range size is an opportunity axis

`range_size` must not be interpreted as a balancing mechanism. It can affect:

- environmental heterogeneity;
- number of populations;
- mutation/origin opportunity;
- persistence opportunity;
- literature and observation effort.

Therefore the main models should decompose range size into measurable mediators where possible.

## Sampling design

### Tier 1: broad phylogenetic background

Sample species across major angiosperm orders/families, stratified by pollen vector. Target a minimum of:

- 300 biotically pollinated species;
- 100 wind-pollinated species;
- all feasible ambophilous and water-pollinated species in the sampled clades.

Do not hand-pick monomorphic controls. Draw them from a predefined taxonomic sampling frame and then code the evidence state.

### Tier 2: matched contrasts

For every verified polymorphic species, select 1-3 phylogenetically close species with similar growth form and geography where possible. These matched contrasts improve power but do not replace the broad background sample.

## Core comparative models

### Model 1: any documented species-level colour variation

`variation_present ~ pollen_vector + display_opportunity + range_size + sampling_effort + phylogeny`

### Model 2: spatial organization conditional on variation

Ordinal or multinomial outcome:

`geographic_mosaic / mixed / local_coexistence`

Predicted by:

`environmental_heterogeneity + connectivity + pollinator_diversity + specialization + range_size + sampling_effort + phylogeny`

### Model 3: biotic versus abiotic selection architecture

Test whether biotic pollination predicts colour variation after restricting the dataset to species with comparable visible display opportunity.

The theoretical expectation is not simply `animal pollination -> polymorphism`. Biotic pollination creates a route for colour-mediated selection, but the sign of that selection can be balancing or directional. Strong consistent pollinator selection can therefore favour monomorphism, whereas heterogeneous or frequency-dependent selection can favour polymorphism.

## Bias controls

At minimum code:

- number of population-level studies;
- number of occurrence records or equivalent distribution-data proxy;
- number of floras/monographs searched;
- whether colour variation was a focal topic of the source;
- sample size where reported;
- cultivated/ornamental status.

Cultivated variation must not be counted as natural polymorphism unless natural populations are documented separately.

## Immediate implementation sequence

1. Expand the positive-evidence queue across major angiosperm families.
2. Build a taxonomically stratified background sampling frame.
3. Code pollination mode and display opportunity before coding colour outcome.
4. Treat `unknown` explicitly rather than converting missing evidence to monomorphism.
5. Run the first analysis on broad `variation_present`, then the spatial-organization analysis on the verified polymorphic subset.
