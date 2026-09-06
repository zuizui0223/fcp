# RGFCA methodological novelty positioning

Status: working positioning note; **not** a claim that generic Monte Carlo, bootstrap aggregation, spatial thinning or spatial bootstrap is new.

## Claim to avoid

> We introduce the first ecological method to repeatedly resample spatial biodiversity records.

That claim is not defensible. Ecological modelling already uses bootstrap aggregation, range bagging, block bootstrap, spatial thinning/balancing and spatially aware ensemble methods.

Relevant precedents include:

- Drake JM. 2015. *Range bagging: a new method for ecological niche modelling from presence-only data*. Journal of the Royal Society Interface 12:20150086. doi:10.1098/rsif.2015.0086. Range bagging repeatedly subsamples occurrence records/environmental dimensions and aggregates an ensemble of marginal niche estimates.
- Fithian W et al. 2015. *Bias correction in species distribution models: pooling survey and collection data for multiple species*. Methods in Ecology and Evolution 6. doi:10.1111/2041-210X.12242. Spatial block bootstrap is used to accommodate spatial autocorrelation in inference.
- Steen VA et al. 2021. *Spatial thinning and class balancing: Key choices lead to variation in the performance of species distribution models with citizen science data*. Methods in Ecology and Evolution. doi:10.1111/2041-210X.13525. Spatial bias and spatial thinning/balancing of opportunistic citizen-science data are explicit methodological concerns.
- Recent spatial-bagging work also incorporates spatial dependence or spatial weighting into ensemble learning; therefore 'bagging a map' alone is not a novelty claim.

## What is methodologically distinctive here

The RGFCA contribution is the **estimand and workflow combination**, not resampling in isolation.

### 1. The base learner is a bounded multispecies world-map realization

Each observed realization contains a balanced subset of eligible species and a balanced fixed number of photographs per species, retained at true coordinates. The inferential object is not a species-distribution prediction and not a bootstrap confidence interval around one fitted model; it is a realization-specific global field of flower-colour discontinuity conditional on geographic sampling opportunity.

### 2. Equal-species contribution is enforced before global mapping

Data-rich species cannot dominate merely through record count. Species and photographs are scheduled so long-run inclusion counts are balanced, while one computational unit remains bounded.

### 3. The field is opportunity corrected

The numerator is colour-discontinuity support and the denominator is the geographic opportunity for sampled species to contribute an edge at that location. Thus recurrent signal is distinguished from regions that simply contain many sampled edges.

### 4. The primary output is recurrence of geographic zones

The atlas estimates `P_hot(x)`: the probability that a cell reappears as a top-decile colour-discontinuity hotspot across repeated world-map realizations. Persistent connected components are defined prospectively and receive neutral IDs before named geography or mechanisms are inspected.

### 5. The null preserves the whole sampled geometry

For null realizations, species membership, selected photographs, coordinates, graph geometry and the Monte Carlo schedule stay fixed. Only complete colour vectors are permuted within species. The target comparison is therefore recurrent observed colour geography versus recurrent colour geography expected from the same opportunistic sampling geometry with within-species colour-location association removed.

### 6. Global discovery is itself audited as a repeated sampling problem

The species universe is not assumed from one API draw. Repeated metadata-only discovery is explicitly audited for accumulation, page reuse and independent-round overlap before flower-colour pixels are opened. The V1 cache diagnostic is part of this methodological point: repeated calls are not automatically independent realizations.

## Defensible novelty sentence

> We treat global flower-colour biogeography as a balanced repeated-atlas problem: bounded multispecies realizations of public flower photographs are repeatedly mapped at their true coordinates, converted to opportunity-corrected colour-discontinuity fields, and summarized by the recurrence of geographic hotspots against a species-conditioned matched null.

A stronger 'first' claim should be made only after a formal literature search specifically targeting multispecies trait-boundary atlases, repeated biodiversity maps and recurrence-based biogeographic zoning.
