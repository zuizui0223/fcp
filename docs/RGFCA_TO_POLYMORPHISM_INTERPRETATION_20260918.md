# RGFCA to species-level polymorphism — conceptual interpretation

Date: 2026-09-18 JST

Status: reporting-only interpretation note. This document does not rerun RGFCA or the current polymorphism analyses.

## 1. What RGFCA is

RGFCA stands for **Repeated Global Flower-Colour Atlas**.

It was designed as a global, image-first biogeographic framework for asking whether within-species flower-colour discontinuities recur in the same broad geographic regions across independent plant species.

The core methodological unit was not one pooled global dataset. It was a sequence of bounded, balanced world-map realizations:

- globally eligible species were repeatedly sampled under equalized inclusion schedules;
- a fixed number of classifiable photographs was sampled within each included species;
- photographs were kept at their observed geographic coordinates;
- within-species flower-colour discontinuity was estimated without pooling species identities;
- each species contributed equally to an opportunity-corrected global discontinuity field;
- hotspot recurrence across repeated realizations was compared with a species-conditioned null that preserved species identities, coordinates, graph geometry, sampling schedules and measurement missingness while breaking colour-location association within species.

The primary RGFCA shared-geography estimand therefore asked:

> Do independent species repeatedly place strong within-species flower-colour discontinuities in the same broad geographic regions more often than expected from the same geographic sampling opportunity?

RGFCA also retained a separate species-specific spatial estimand:

`rho_i = Spearman(pairwise geographic distance, pairwise flower-colour dissimilarity)`

which asks how strongly colour is geographically organized within species, regardless of whether different species share the same boundary.

## 2. What RGFCA achieved

RGFCA produced the sampling and measurement infrastructure that now underlies the active polymorphism paper:

- repeated equal-area metadata discovery and the 42,111-species opportunity frame;
- metadata-only capacity scans;
- fixed-photo high-depth species sampling;
- discovery/reserve species-disjoint cohorts;
- observer caps and geographic maximin selection;
- location-blind image measurement;
- flower/background paired measurement in reserve;
- species-specific within-species spatial organization statistics;
- species-conditioned spatial nulls;
- explicit observation-bias and measurement-missingness diagnostics.

The high-depth RGFCA resource contains 1,000 species × 100 raw photographs, divided into 500-species discovery and reserve source cohorts. After the frozen >=40-classifiable gate, 369 discovery and 363 reserve species enter the current D-based analyses.

## 3. What did not become the positive core of the current paper

The original RGFCA goal was a **shared global geography** of strong within-species colour transitions.

That estimand was not promoted into the current paper's main positive claim. The upstream record includes:

- primary recurrent-field G1 concentration: p = 0.070;
- species-disjoint commonness: p = 0.856;
- reserve whole-frame flower-minus-background differential: mean rho = 0.0044773, p = 0.087;
- no environmental process block passing its fixed five-block Holm gate;
- later sharedness qualification failing its prospective power/identifiability gate.

These results do not establish absence of flower-colour biogeography. They show that a single shared cross-species boundary field was not supported/identifiable strongly enough under the frozen RGFCA designs to serve as the current paper's biological spine.

## 4. The key pivot

RGFCA revealed an important asymmetry:

- **shared geographic location across species** was weak or unsupported;
- **species-specific spatial organization** was measurable and heterogeneous.

The active polymorphism paper changes the estimand from:

`Where in the world do species share colour boundaries?`

to:

1. `How much within-species flower-colour diversity does each species carry?` — D;
2. `Along which direction in colour space is that diversity organized?` — H2;
3. `Do species with greater diversity also show stronger internal geographic organization?` — D–spatial association;
4. `Can the between-species differences be reduced to broad ancestry or sampled geographic opportunity?` — H3.

This is not a new data-generation campaign for the legacy discovery/reserve analyses. It is a reparameterization of the same high-depth RGFCA measurement resource into species-level phenotypes and species-level comparative estimands.

## 5. What the current paper inherits from RGFCA

### Directly inherited data

- discovery measured-photo table;
- reserve measured-photo table;
- species coordinates and observer IDs;
- four-state coarse colour calls;
- normalized nine-colour palette measurements;
- technical-failure and ambiguity states;
- reserve paired flower/background measurements;
- species-specific spatial organization statistics.

### Directly inherited design principles

- outcome-blind species/photo acquisition;
- fixed raw-photo denominators;
- observer contribution caps;
- geographic maximin photo selection;
- location-blind measurement;
- species-disjoint reserve validation;
- species-conditioned nulls;
- no-rescue/no-replacement principles.

### Not inherited as a positive claim

- one universal shared global flower-colour boundary;
- globally recurrent named zones;
- a confirmed climate/pollinator/terrain mechanism;
- unbiased global prevalence.

## 6. How the fresh third cohort differs

The prospective third cohort is not part of the original 1,000-species RGFCA high-depth resource.

It was selected from the later 42,111-species RGFCA-derived opportunity frame after excluding:

- all 1,000 legacy discovery/reserve species;
- all 500 P500-selected species.

The resulting candidate universe contained 3,230 species. Exactly 500 were selected outcome-blind, and fresh metadata yielded 499 species × 100 authorized rows. Previously used observation/photo IDs were excluded.

Thus the current paper combines:

- **legacy RGFCA data** for H1, legacy H2 discovery/audit, D–spatial organization and H3;
- **new species and fresh photo IDs** for prospective H2 confirmation.

The third cohort is still within the same iNaturalist source/opportunity universe and measurement system, so it is species/photo-disjoint prospective confirmation rather than independent-source replication.

## 7. Biological interpretation of the pivot

The current evidence changes the level at which generality is found.

The strongest generality is not a common geographic boundary shared by many species. It is instead:

- a reproducible species-level amount of flower-colour diversity;
- a recurrent white-versus-nonwhite achromatic–chromatic geometry;
- a positive relationship between the amount of diversity and the strength of species-specific geographic organization.

A concise conceptual interpretation is therefore:

> **The common rule is more evident in phenotype space than in geographic space. Species repeatedly use a similar achromatic–chromatic direction of flower-colour variation, but the locations and processes that maintain or sort that variation are species-specific.**

This interpretation does not identify the maintenance mechanism. Climate, pollinator turnover, gene flow, dispersal, drift, demographic history and mating system remain candidate processes rather than demonstrated causes.

## 8. Relation to the current paper's two kinds of contribution

### Methodological contribution

RGFCA supplies the scalable sampling/measurement architecture:

- repeated metadata discovery;
- equal-area geographic opportunity;
- fixed high-depth photo budgets;
- observer and sampling controls;
- blind image measurement;
- matched nulls and background controls;
- explicit measurement missingness.

The current paper adds the species-level phenotype and confirmation layer:

- observer-disjoint reliability of D;
- separation of target discovery from prospective confirmation;
- construction-preserving H2 null;
- frozen third-cohort q_white/W test.

### Ecological contribution

The current paper extracts biological regularities from that infrastructure:

- within-species flower-colour diversity differs reproducibly among species;
- the recurrent colour-space geometry is white versus non-white;
- more polymorphic species show stronger species-specific geographic colour organization;
- these differences are not simply explained by sampled geographic span or broad tree-wide phylogenetic conservation.

## 9. Claim boundary

Do not describe the current paper as if RGFCA had demonstrated a universal global flower-colour boundary and the present study merely explains it.

The supported relation is the reverse:

> RGFCA created the global sampling/measurement framework and showed that species-specific spatial information was recoverable even when cross-species shared geography was not. The current paper uses that heterogeneity as the biological object of study.
