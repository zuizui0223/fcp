# RGFCA methodological novelty positioning

Status: working positioning note; **not** a claim that generic Monte Carlo, bootstrap aggregation, spatial thinning or spatial bootstrap is new.

Updated 7 September 2026 after a bounded primary-source audit. Read/access
ceilings and full references are in the [image/ecology audit](RGFCA_IMAGE_ECOLOGY_LITERATURE_AUDIT.md)
and [statistical audit](RGFCA_STATISTICAL_LITERATURE_AUDIT.md). This is not a
systematic review. The active [manuscript](RGFCA_MANUSCRIPT.md) now acknowledges
the direct image-derived floral precedents as well as the statistical ones.

## Claim to avoid

> We introduce the first ecological method to repeatedly resample spatial biodiversity records.

That claim is not defensible. Ecological modelling already uses bootstrap aggregation, range bagging, block bootstrap, spatial thinning/balancing and spatially aware ensemble methods.

Relevant precedents include:

- [Drake (2015)](https://doi.org/10.1098/rsif.2015.0086): range bagging targets environmental niche support.
- [Fithian et al. (2015)](https://doi.org/10.1111/2041-210X.12242): joint survey/collection modelling and spatial block bootstrap.
- [Steen et al. (2021)](https://doi.org/10.1111/2041-210X.13525): thinning and class-balancing performance depends on the modelling/evaluation setting.

Image-first flower-colour work also has direct precedents:

- [Luong et al. (2023)](https://doi.org/10.1002/aps3.11546): landscape colour analysis in *Erysimum*.
- [McKenzie, Church and Hopkins (2026)](https://doi.org/10.1086/739413): high-throughput geographic colour phenotyping in *Monarda fistulosa*.
- [McKenzie, Berardi and Hopkins (2025)](https://doi.org/10.1016/j.cub.2025.03.035): multispecies colour/phenology and seasonal pollinator context in North America.

Accordingly, neither photograph volume, within-species geographic colour
variation nor adding pollinator distributions is by itself a priority claim.
Published versions are cited; unread final methodological details are not
inferred from earlier preprints.

## What is methodologically distinctive here

The RGFCA contribution is the **estimand and workflow combination**, not resampling in isolation.
The points below describe the proposed design, not demonstrated superiority to
those methods. Discovery currently supports only a small exploratory
distance-colour association. Sharedness, flower specificity, independent
replication and ecological mechanism do not become established by this wording.

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

This does not preserve the colour field's spatial autocorrelation or remove all
selection/measurement bias. Testing one field's spatial structure and testing
association between two autocorrelated fields need different nulls.
[Guillot and Rousset (2013)](https://doi.org/10.1111/2041-210X.12018).

### 6. Global discovery is itself audited as a repeated sampling problem

The species universe is not assumed from one API draw. Repeated metadata-only discovery is explicitly audited for accumulation, page reuse and independent-round overlap before flower-colour pixels are opened. The V1 cache diagnostic is part of this methodological point: repeated calls are not automatically independent realizations.

## Defensible novelty sentence

> We treat global flower-colour biogeography as a balanced repeated-atlas problem: bounded multispecies realizations of public flower photographs are repeatedly mapped at their true coordinates, converted to opportunity-corrected colour-discontinuity fields, and summarized by the recurrence of geographic hotspots against a species-conditioned matched null.

A stronger 'first' claim should be made only after a formal literature search specifically targeting multispecies trait-boundary atlases, repeated biodiversity maps and recurrence-based biogeographic zoning.
