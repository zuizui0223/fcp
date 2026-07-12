# Theory-guided comparative synthesis protocol

## Objective

Test why some flowering-plant species are recorded as floral-colour polymorphic whereas others are recorded as monomorphic, while separating ecological predictors from origin opportunity and observation bias.

## Unit of analysis

Primary unit: species.

Secondary nested unit when data permit: population within species.

The species-level response is not treated as a perfect observation of the underlying evolutionary state.

## Response variables

### Primary

`observed_polymorphism`:

- 1: at least two discrete, naturally occurring floral-colour morphs reported within a taxonomic species;
- 0: explicitly surveyed and reported as monomorphic;
- NA: absence of evidence only.

A crucial rule is that an unreported polymorphism is not automatically coded 0.

### Secondary

- number of discrete colour morphs;
- within-species colour variance when quantitative spectral data exist;
- number/proportion of polymorphic populations;
- morph-frequency evenness.

## Four-filter covariate architecture

### A. Origin opportunity

Proxies:

- geographic range size;
- census/occupancy proxies;
- life span and generation time;
- hybridization/introgression evidence;
- ploidy;
- known pigment-pathway lability when available.

Prediction: more opportunity for origin can increase recorded polymorphism even without stronger balancing selection.

### B. Invasion / balancing environment

Proxies:

- spatial environmental heterogeneity across the range;
- pollinator-community heterogeneity;
- habitat heterogeneity;
- climatic heterogeneity;
- strength/proxies of negative frequency dependence where directly measured.

Prediction: heterogeneity is more likely to associate with polymorphism when not erased by high effective connectivity.

### C. Persistence / homogenization

Proxies:

- dispersal syndrome;
- pollen and seed dispersal distances;
- habitat connectivity;
- population fragmentation;
- selfing rate / mating system;
- effective population size proxies.

Primary theory-guided interaction:

`heterogeneity × connectivity`.

The sign cannot be reduced to a universal raw-connectivity main effect because gene flow can both spread variants and homogenize local adaptation.

### D. Observation process

Mandatory where possible:

- number of populations sampled;
- number of individuals sampled;
- number of independent publications/records;
- herbarium or citizen-science record count;
- whether colour was an explicit study target;
- spectral measurement versus human colour category.

Prediction: recorded polymorphism rises with sampling effort even when the true evolutionary state is held constant.

## Primary model family

A phylogenetically informed latent-state model is preferred.

Observation model:

\[
Y_i^{obs}\sim Bernoulli(z_i q_i),
\]

where \(z_i\) is latent true polymorphism and \(q_i\) is detection probability.

Ecological state model:

\[
logit\,P(z_i=1)=
\alpha+\beta_H H_i+\beta_C C_i+\beta_{HC}H_iC_i+
\beta_T T_i+\beta_{TB}T_iB_i+
\beta_O O_i+u_{phylogeny}.
\]

Detection model:

\[
logit(q_i)=\eta_0+\eta_1\log(1+n_{pop,i})+\eta_2\log(1+n_{ind,i})+\eta_3 I(\text{explicit colour study}).
\]

If the data cannot identify a latent detection model, use sensitivity analyses rather than pretending observed zeros are true zeros.

## Confirmatory hypotheses

H1. Spatial environmental heterogeneity is positively associated with floral-colour polymorphism after observation effort is controlled.

H2. The heterogeneity effect depends on connectivity/dispersal structure.

H3. Temporal environmental variance alone has no universal positive effect; positive associations should be strongest in species with buffering/storage-like life histories or spatial structure.

H4. Sampling effort strongly predicts recorded polymorphism and must be modeled.

H5. Range size has two separable pathways: ecological heterogeneity and opportunity/detection. Analyses should not interpret a positive range-size coefficient as balancing selection without decomposition.

## Negative controls / falsification

- Replace environmental heterogeneity with randomly permuted within-biome heterogeneity.
- Test whether sampling effort predicts traits that should not depend on ecological balancing but are similarly under-recorded.
- Refit after excluding species whose polymorphism is known only from horticultural cultivars.
- Refit under strict and broad definitions of polymorphism.
- Test whether the main result survives equal-effort subsampling.

## Inclusion boundary

Include naturally occurring, heritable or repeatedly stable floral-colour variants within species.

Exclude by default:

- ontogenetic colour change within a flower;
- post-pollination colour change;
- purely environmentally induced transient colour shifts;
- horticultural cultivars without natural populations;
- taxonomic mixtures or unresolved species complexes unless explicitly modeled.

## Immediate data-building strategy

1. Build a positive set of well-documented polymorphic species.
2. For each positive species, select phylogenetically matched candidate controls with explicit floral descriptions and comparable study effort.
3. Record evidence for both presence and genuine absence rather than coding database silence as monomorphism.
4. Expand to a broader sample only after the coding protocol yields acceptable inter-rater agreement.

This matched-first strategy is preferable to a naive all-angiosperm scrape because the dominant initial threat is ascertainment bias, not sample size.
