# FCP image-first atlas: global boundary data sources and prospective inference design

## Scope and decision boundary

This note evaluates primary data sources for the next FCP atlas gates. It does not inspect candidate image pixels, revise the frozen six-species Stage A/B result, revise the 34-species comparison, or reopen the completed three-species negative validation. The already frozen 50-species, 20,200-observation metadata cohort remains closed pending estimator qualification.

The purpose is to predeclare a defensible route from a larger image atlas to one of three ordered conclusions:

1. flower-colour transitions show cross-species geographic concentration;
2. if geographic concentration is not supported, the same frozen colour fields show or do not show concordance with independently frozen environmental boundaries;
3. independently of those outcomes, the same fields show or do not show concordance with a separately frozen pollinator-biogeographic overlay, provided that overlay passes its own coverage gate.

“Direction change” must therefore mean moving through this fixed decision tree, not changing cohorts, scales, overlays, or null models until a small p-value appears. A branch can end as `supported`, `not_supported`, or `not_evaluable`. A negative result is a conclusion.

## Recommendations at a glance

| Resource or method | Decision | Atlas role | Data availability and claim ceiling |
|---|---|---|---|
| iNaturalist Open Data monthly snapshot | **Adopt** | Primary bulk image and metadata acquisition after estimator qualification | Dated metadata snapshots and resized licensed images are openly addressable. This is a community-observation sample, not a census of all flowers or plants. |
| iNaturalist API | **Hold** | Metadata feasibility, small refreshes, and ID resolution only | Official guidance says the API is not for bulk download. It cannot be the primary path for a multi-cohort image atlas. |
| GBIF occurrence download with DOI | **Adopt** | Immutable occurrence provenance for independent plant or pollinator overlays | A download DOI fixes the mediated records and contributing datasets, not their accuracy, completeness, or absence process. |
| WorldClim 2.1 | **Adopt** | Primary macroclimate boundary overlay | Global approximately 1-km 1970–2000 climatologies. Supports concordance with long-term climate gradients, not event-date weather or climatic causation. |
| CHELSA 2.1 | **Hold** | One frozen orographic-climate sensitivity | Kilometer-scale global climate predictors are available, but choosing between CHELSA and WorldClim after seeing colour results would create an analysis multiverse. |
| Copernicus DEM GLO-30 | **Adopt** | Primary terrain overlay after frozen aggregation | Global 30-m DEM with documented access/licensing. Supports elevation/slope/ruggedness concordance, not a historical or mechanistic explanation. |
| ESA WorldCover 2021 | **Adopt** | Primary categorical vegetation/land-cover overlay | Global 10-m, 11-class product. Its 2021 snapshot is not observation-date habitat and cannot identify natural vegetation or pollinator habitat by itself. |
| RESOLVE Ecoregions 2017 | **Adopt** | Primary independent categorical ecoregion boundary | Downloadable CC BY 4.0 polygons for 846 terrestrial ecoregions nested in 14 biomes and 8 realms. Concordance is descriptive biogeography, not evidence that ecoregions caused flower colour. |
| EarthEnv 1-km consensus land cover | **Hold** | Prespecified legacy sensitivity only | Continuous class prevalence is useful at atlas scale, but inputs are older and the product is CC BY-NC 4.0. |
| Fresh GBIF *Bombus* occurrence regionalization | **Hold pending coverage gate** | Candidate pollinator-biogeographic overlay | Can be made independent of flower colour by freezing a DOI, taxonomy, filters, grid, and clustering before colour joins. Presence-only bias and uneven global coverage can make the overlay `not_evaluable`. |
| Orr et al. global bee checklist/occurrence synthesis | **Hold** | External coverage benchmark and coarse bee-region sensitivity | The synthesis documents severe spatial and taxonomic bias in public occurrence databases. It is not a ready-made high-resolution *Bombus* boundary layer. |
| NHM *Bombus* maps | **Do not adopt for inference** | Taxonomic/biogeographic background | Authoritative context, but not a versioned, machine-readable, globally complete inferential surface. |
| GloBI interaction data | **Hold for exploratory guild annotation** | Possible post-confirmatory annotation of known plant–pollinator interactions | Stable versioned data exist, but heterogeneous interaction records do not supply standardized survey effort or absences. They cannot be a global pollinator-boundary primary test. |
| Repeated random cohorts with one nested null | **Adopt** | Scale-out and repeatability assessment | All cohorts, seeds, exclusions, stopping rules, and the single across-cohort statistic must be frozen. Repeating tests until one rejects is prohibited. |
| Unrestricted permutation for environment/pollinator concordance | **Do not adopt** | None | It destroys spatial autocorrelation and can make smooth spatial fields appear too surprising. Use a validated spatially constrained null. |

## 1. iNaturalist and GBIF acquisition

### 1.1 Primary image source: iNaturalist Open Data

The official iNaturalist Open Data repository describes a monthly snapshot containing observation, observer, photo, taxon, and project tables, with image objects addressable by photo ID. It provides `large` images with maximum dimension 1,024 px and `original` images with maximum dimension 2,048 px. It also warns that live deletions and identification changes can leave temporary mismatches between the monthly tables and object store. At the time of the official documentation, the bucket contains more than 400 million licensed organism photographs ([iNaturalist Open Data documentation](https://github.com/inaturalist/inaturalist-open-data)).

Each photograph retains its own Creative Commons licence. The official documentation states that photographers retain copyright unless a photo is CC0, and that use must follow the recorded licence and attribution requirements ([iNaturalist Open Data licensing](https://github.com/inaturalist/inaturalist-open-data#licenses)). Therefore the frozen photo manifest must retain at least snapshot date, observation ID, photo ID, object key or URL, photo licence, observer attribution fields, file size, extension, and acquired-file SHA-256.

Recommended implementation:

- use a dated monthly metadata archive rather than a moving `latest` pointer;
- freeze all admitted observation and photo IDs before fetching pixels;
- fetch one declared image size for measurement and retain original URLs only as provenance;
- verify licence, object existence, content length, decode status, and SHA-256 without silently replacing failed photos;
- record live-snapshot mismatches as explicit acquisition failures.

Claim ceiling: iNaturalist records are user-submitted encounters with organisms and one observation may have multiple evidence photos. The open bucket is therefore not “all flower photographs,” and neither a Plantae/Angiospermae identification nor a flowering annotation guarantees that an evaluable flower-tissue ROI exists. Atlas inclusion remains a measurement result, not a metadata fact.

### 1.2 API is for bounded feasibility, not atlas-scale scraping

iNaturalist’s official API practices allow pagination to about 10,000 results for many endpoints and direct users needing more than 10,000 records to exports, the weekly GBIF dataset, or other bulk products. The same page asks clients to stay near one request per second and about 10,000 API requests per day, and warns that downloading more than 5 GB of media per hour or 24 GB per day may result in a permanent block ([iNaturalist API Recommended Practices](https://www.inaturalist.org/pages/api+recommended+practices)). The developer page likewise says the API supports applications rather than scraping and recommends datasets for large-scale use ([iNaturalist Developers](https://www.inaturalist.org/pages/developers)).

Decision: keep the API adapter for metadata-only feasibility and small identifier refreshes. Do not shard requests across IP addresses or use GitHub Actions matrices to evade service limits. Bulk cohort pixels must come from the open-data object store after the cohort manifest is frozen.

### 1.3 GBIF as the immutable occurrence layer

GBIF’s official download API creates asynchronous downloads and assigns a DOI once a download succeeds ([GBIF occurrence download API](https://techdocs.gbif.org/en/data-use/api-downloads)). Darwin Core Archive downloads can include a `multimedia.txt` table, while Simple and Parquet formats expose interpreted occurrence fields ([GBIF download formats](https://techdocs.gbif.org/en/data-use/download-formats)).

GBIF standardizes dataset licences to CC0, CC BY, or CC BY-NC, but explicitly does not guarantee the accuracy of mediated biodiversity data ([GBIF terms](https://www.gbif.org/terms)). A GBIF DOI is consequently a provenance and reproducibility object, not a quality certificate.

Recommended implementation for an independent *Bombus* or other pollinator overlay:

- submit a frozen predicate and retain the returned DOI and download key;
- retain `datasetKey`, `occurrenceID`, taxon fields, basis of record, event date, coordinates, coordinate uncertainty, geospatial issue flags, establishment means where available, and licence;
- exclude fossils, absences, records without coordinates, records above the frozen uncertainty ceiling, and records with disqualifying geospatial flags;
- deduplicate stable occurrence identifiers and declared spatial/date duplicates;
- never replace the DOI after seeing flower-colour results.

## 2. Independent environmental boundary layers

All environmental rasters must be downloaded, versioned, reprojected, and aggregated before joining to colour fields. The join manifest must hash the source files and the derived 100/250/500-km products. Only one primary layer per environmental family is permitted; named alternatives are sensitivities, not a menu.

### 2.1 Macroclimate: WorldClim 2.1 primary

WorldClim 2.1 provides global land climatologies for 1970–2000 at resolutions down to 30 arc seconds, including 19 bioclimatic variables in GeoTIFF format ([WorldClim historical climate data](https://www.worldclim.org/data/worldclim21.html)). Its primary paper describes approximately 1-km global climate surfaces ([Fick & Hijmans 2017](https://doi.org/10.1002/joc.5086)).

Adopt a small, biologically interpretable, predeclared vector rather than all 19 correlated variables. A defensible primary set is annual mean temperature, temperature seasonality, annual precipitation, and precipitation seasonality. Standardize globally using land cells, aggregate to each frozen atlas grid, and define boundary intensity from the norm of the standardized local gradient. The exact variables, standardization population, neighbor operator, and gradient norm must be frozen before colour is joined.

Claim ceiling: this tests whether flower-colour transition intensity is spatially concordant with long-term macroclimate transitions. The 1970–2000 normal is not observation-date weather, cannot separate correlated environmental drivers, and cannot establish climate adaptation.

### 2.2 CHELSA 2.1 as a single frozen climate sensitivity

CHELSA 2.1 offers global kilometer-scale climate data and derived BIOCLIM+ variables; the official source identifies version, status, licence, and model citations ([CHELSA 2.1 model and datasets](https://www.chelsa-climate.org/models/chelsa)). Its foundational climatology paper documents topography-informed downscaling ([Karger et al. 2017](https://doi.org/10.1038/sdata.2017.122)).

Decision: hold one CHELSA-derived replication of the exact WorldClim variable set as an orographic sensitivity. Do not search CHELSA’s much larger derived-variable catalogue after viewing atlas results. Agreement can show robustness to climate product; disagreement lowers the claim ceiling.

### 2.3 Terrain: Copernicus DEM GLO-30

The Copernicus Data Space documents worldwide GLO-30 and GLO-90 DEM products, free-licence availability, source acknowledgement, dataset identifiers, and DOI `10.5270/ESA-c5d3d65` ([Copernicus DEM collection](https://dataspace.copernicus.eu/explore-data/data-collections/copernicus-contributing-missions/collections-description/COP-DEM)).

Adopt GLO-30 to derive elevation, slope, and terrain ruggedness under a frozen algorithm. Aggregate those derived surfaces to the atlas grids and calculate one multivariate terrain-boundary intensity. Because GLO-30 is a digital surface model assembled from multiple elevation sources, do not interpret a local gradient as pure bare-earth topography without validation.

Claim ceiling: concordance with terrain discontinuity; no inference about dispersal barriers, historical isolation, or selection without separate data.

### 2.4 Vegetation/land cover: ESA WorldCover 2021

ESA WorldCover provides global 2020 and 2021 land-cover maps at 10-m resolution based on Sentinel-1 and Sentinel-2 data ([ESA WorldCover](https://esa-worldcover.org/en)). The 2021 product contains 11 land-cover classes and is described as open and free ([ESA 2021 release](https://esa-worldcover.org/en/release-worldcover-map-2021)).

Adopt the 2021 v200 product and aggregate each class to within-cell proportions. Define boundary intensity from compositional dissimilarity between neighboring cells; freeze the distance metric and water/ice handling. Do not treat class boundaries as direct vegetation-species turnover or as conditions at the date of every iNaturalist photograph.

EarthEnv’s global 1-km consensus land-cover product supplies continuous prevalence for 12 classes and may be useful as a legacy sensitivity, but it integrates older remote-sensing inputs, covers land only to 56°S, and is distributed under CC BY-NC 4.0 ([EarthEnv consensus land cover](https://www.earthenv.org/landcover)). It should not displace WorldCover after outcomes are known.

### 2.5 Categorical biogeography: RESOLVE Ecoregions 2017

The official Ecoregions 2017 site provides a downloadable shapefile under CC BY 4.0 ([Ecoregions 2017 data](https://ecoregions.appspot.com/)). The dataset delineates 846 terrestrial ecoregions grouped into 14 biomes and 8 realms ([RESOLVE Ecoregions 2017 catalogue](https://developers.google.com/earth-engine/datasets/catalog/RESOLVE_ECOREGIONS_2017)); the accompanying primary paper explains the updated global framework ([Dinerstein et al. 2017](https://doi.org/10.1093/biosci/bix014)).

Adopt the polygon boundaries as an independent categorical overlay. Use a frozen distance-to-boundary or neighboring-cell class-change statistic, and report realm, biome, and ecoregion resolutions as one prespecified nested family rather than three unrelated opportunities to reject.

Claim ceiling: ecoregions synthesize climate, habitat, and biotic assemblages. Alignment is an omnibus biogeographic association and cannot identify which component generated flower-colour transitions.

## 3. Pollinator biogeography

### 3.1 No ready-to-use, high-resolution global *Bombus* region layer was verified

The Natural History Museum’s *Bombus* resource is authoritative taxonomic and biogeographic context and describes roughly 249 species with strong Northern Hemisphere and montane concentration ([NHM worldwide *Bombus* patterns](https://www.nhm.ac.uk/research-curation/research/projects/bombus/decline.html); [NHM Paul Williams profile](https://www.nhm.ac.uk/our-science/people/paul-williams.html)). It is not, however, a current versioned polygon/raster product with a frozen global sampling model. It should be cited for background, not converted by hand into an inferential boundary.

The most defensible candidate is therefore a prospectively derived *Bombus* community-turnover regionalization from a frozen GBIF download. It is independent of the flower-colour outcome only if the DOI, taxonomic backbone, filters, grid, effort thresholds, clustering method, number-of-region rule, and uncertainty analysis are frozen before any colour join.

### 3.2 Mandatory *Bombus* coverage gate

Orr et al. compiled 5,857,811 bee occurrence records from several public sources, but fewer than 16% passed all filters. Their public-data map strongly followed sampling effort, with major gaps in Asia, the Middle East, and Africa; the checklist-based data were more representative ([Orr et al. 2021, article and supplements](https://www.sciencedirect.com/science/article/pii/S0960982220315967)). That primary study is direct evidence that a visually detailed global bee map can still be dominated by digitization effort.

Before a *Bombus* overlay is admitted, freeze and pass at least:

- minimum occupied cells and minimum records per retained species;
- maximum contribution by any dataset, country, observer/source, and decade;
- coverage by realm and by atlas opportunity cells;
- taxonomic reconciliation against a declared *Bombus* checklist;
- sensitivity to specimens versus observations and to predeclared date windows;
- stability of boundaries under spatial block bootstrap and dataset-leave-one-out analysis;
- independence from iNaturalist flower-photo IDs and from all measured colour values.

If a realm fails coverage, it is `not_evaluable`; the rule must not interpolate a global boundary merely to fill the map. A passed layer supports only concordance with an occurrence-derived *Bombus* assemblage transition. It does not establish that bumblebees visited the photographed plants, saw the measured colours, selected the floral phenotype, or caused the transition.

### 3.3 Broader bees and interaction databases

Orr et al.’s country/checklist synthesis can be frozen as a coarse broad-bee sensitivity or as an external benchmark for whether a new *Bombus* occurrence grid reproduces known coverage gaps. It should not be represented as a high-resolution pollinator boundary.

GloBI integrates open species-interaction datasets, retains source citations, and publishes stable versioned snapshots with DOI `10.5281/zenodo.3950589` for data-intensive work ([GloBI data access](https://www.globalbioticinteractions.org/data); [Poelen et al. 2014](https://doi.org/10.1016/j.ecoinf.2014.08.005)). Its records are heterogeneous evidence of reported interactions, not standardized global surveys with known detection effort. GloBI may annotate which atlas plant taxa have recorded pollinator guilds after the primary tests, but it cannot define the primary global pollinator boundary or turn missing interactions into absences.

A 2026 primary compilation explicitly documents geographic biases and gaps in sampled plant–pollinator networks and releases its metadata and reproducible workflow ([Brito et al. 2026 data](https://doi.org/10.5061/dryad.ffbg79d8z)). This supports a conservative decision: do not use Web-of-Life-style network availability as if it were uniform global interaction coverage.

## 4. Repeated random cohorts and valid nulls

### 4.1 Cohort repetition is one designed experiment, not repeated chances to reject

Prospective registration separates confirmatory tests from hypotheses generated by exploration ([Nosek et al. 2018](https://doi.org/10.1073/pnas.1708274114)). A multiverse analysis can expose sensitivity to reasonable processing choices, but only when the complete set of choices is shown rather than selectively reported ([Steegen et al. 2016](https://doi.org/10.1177/1745691616658637)).

Recommended prospective cohort design:

1. Build one metadata-only eligible species frame from a dated iNaturalist Open Data snapshot under the existing positional-accuracy, licence, flowering-annotation, observer/cell, spatial-cell, month, and geometry rules.
2. Freeze the eligible-frame hash, exclusions, genus cap, sample tiers, seed-generation rule, number of cohorts, overlap rule, and a hard compute budget before pixels are opened.
3. Prefer disjoint species cohorts sampled without replacement. If metadata feasibility shows that disjoint cohorts are impossible, freeze one explicit overlapping-cohort design and include the overlap in the null; do not switch after colour results.
4. Run every frozen cohort. No stopping after the first supportive or non-supportive cohort.
5. Give each species equal weight and each cohort a predeclared weight. Report the full distribution of cohort effects and one global summary statistic, not a count of nominally significant cohorts.
6. Generate a single Monte Carlo null by rerunning the complete hierarchy—within-species colour randomization or surrogate generation, field reconstruction, boundary scoring, cohort aggregation, and cross-cohort summary—for every null replicate.
7. Treat cohort-to-cohort heterogeneity as an estimand. A result concentrated in one geography, clade, data-rich region, or cohort cannot be called global repeatability.

The final numeric number of cohorts must follow metadata-only feasibility and a prospective precision/compute calculation. It must not follow observed colour effect sizes. A defensible starting target is at least five frozen cohorts, but this is not yet a contract value.

### 4.2 Monte Carlo p-values

Phipson and Smyth show that randomly sampled permutation p-values should not be zero and derive valid exact calculations for Monte Carlo tests ([Phipson & Smyth 2010](https://doi.org/10.2202/1544-6115.1585)). Use

`p = (b + 1) / (B + 1)`,

where `b` is the number of null statistics at least as extreme as the observed statistic and `B` is the number of random null replicates. Freeze `B`, the seed stream, tie rule, and one-/two-sided direction. For a 0.01-scale inferential threshold, 9,999 replicates give a minimum attainable p-value of 0.0001 and should remain the minimum unless a precision calculation justifies more.

### 4.3 Spatially constrained null for overlay concordance

The existing within-species random-mark null is appropriate for asking whether measured colour is spatially organized on the frozen observation graph. It is not automatically sufficient for asking whether two already spatially smooth fields align: unrestricted shuffling can remove spatial autocorrelation and yield a null that is too narrow.

Spatial randomization research shows that the permutation mechanism must encode the relevant autocorrelation. The Floating Grid Permutation Technique was developed to control spatial autocorrelation in ecological randomization tests ([Radersma & Sheldon 2015](https://doi.org/10.1111/2041-210X.12390)). Moran spectral randomization can preserve multiscale autocorrelation for irregularly spaced data and had correct or conservative type-I error under stationary simulated fields, but was sensitive to linear trends ([Wagner & Dray 2015](https://doi.org/10.1111/2041-210X.12407)).

Required estimator-qualification simulation before adopting an environment or pollinator concordance test:

- simulate null colour fields on the exact species graphs with the observed opportunity masks;
- include stationary fields, broad linear/latitudinal trends, range-edge truncation, uneven sampling, season strata, and observer/cell clustering;
- compare unrestricted within-species mark permutation, a frozen Moran-spectral surrogate, and any proposed shift/rotate or graph-constrained alternative;
- require calibrated false-positive rates across 100/250/500-km scales and adequate recovery of injected boundary alignments;
- freeze one primary null before real atlas colour–overlay joins; alternatives remain diagnostics.

If no null controls false positives under the atlas geometry, environmental and pollinator concordance are `not_evaluable`, not rescued with a simpler permutation.

### 4.4 Multiplicity across pivots, overlays, and scales

The three inference domains and every overlay are known in advance, so their multiplicity exists even if only one result is emphasized. Use one frozen familywise procedure. Recommended implementation:

- construct one maximum statistic across geographic concentration, the declared environmental families, and the admitted pollinator family inside each complete null replicate;
- calculate adjusted global and branch-level p-values from that same joint null;
- retain 100 km as primary and 250/500 km as mandatory sensitivities, not three independent claims;
- retain WorldClim, terrain, WorldCover, and ecoregions as a declared environmental family; do not add variables after unblinding;
- admit the pollinator family only if its colour-blind coverage gate passes, and record failure as `not_evaluable` before colour joins.

This joint-null approach permits the narrative to move from spatial to environmental to pollinator results without turning the move into a new uncounted chance to reject.

## 5. Proposed gate sequence

1. **Estimator qualification:** pass an independent flower-tissue ROI benchmark and downstream signal-recovery simulations already required by the atlas protocol.
2. **Snapshot freeze:** choose and hash a dated iNaturalist Open Data snapshot; never use `latest` in the final bundle.
3. **Random-cohort feasibility:** construct the eligible species frame and freeze the number, size, seed schedule, overlap rule, and stop conditions without image pixels.
4. **Environmental freeze:** download and hash WorldClim 2.1, Copernicus DEM GLO-30, WorldCover 2021 v200, and RESOLVE Ecoregions 2017; generate all scale products without colour.
5. **Pollinator feasibility:** request a citable GBIF *Bombus* download, freeze filtering and regionalization, and make a colour-blind `pass`/`not_evaluable` coverage decision.
6. **Null qualification:** choose the spatially constrained concordance null using simulations only.
7. **Pixel opening and measurement:** run all frozen image cohorts with explicit acquisition and ROI failures; no replacements informed by colour.
8. **One complete inference run:** compute geographic, environmental, and admitted pollinator statistics plus one joint multiplicity-corrected null.
9. **Publication outcome:** publish support, non-support, heterogeneity, and `not_evaluable` branches together with the species-free display and species-conditioned inference.

## Final claim ceiling

At maximum, the expanded atlas can establish that continuous, automated flower-candidate colour transitions are repeatably spatially organized within sampled species and are more geographically concentrated, or more concordant with independently defined environmental or pollinator-biogeographic transitions, than under the frozen species-conditioned spatial null.

It cannot establish that the atlas includes all global flowering plants, that a computer-generated candidate ROI is flower tissue without benchmark evidence, that an environmental layer is the causal driver, that a recorded pollinator visited the photographed flower, that pollinators selected the observed colours, or that any detected boundary is universal. Those causal and universal claims require independent trait validation, interaction data, and a design that identifies mechanism rather than spatial concordance.
