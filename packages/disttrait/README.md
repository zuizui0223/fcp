# disttrait

`disttrait` is a small Python package for **validated species-level distributional trait inference** from repeated individual observations.

It extracts the general inferential core used by the FCP flower-colour polymorphism study without hard-coding flower colours, iNaturalist, or the RGFCA world-map application.

## v0.1 scope

The package provides reusable components for:

- species-level categorical diversity using Gini-Simpson diversity;
- observer-disjoint split-half reliability;
- Hellinger-transformed two-mode trait geometry;
- fixed contrast / axis-alignment statistics;
- construction-preserving row permutation within fixed strata;
- pairwise great-circle distances;
- pairwise Jensen-Shannon trait dissimilarity;
- species-specific spatial organization `rho_i`;
- matched trait-minus-background spatial organization;
- species-level distribution-versus-spatial association;
- partial-rank adjustment with matched spatial nulls.

The package deliberately does **not** claim a new standalone statistic. Most components are established statistics. The reusable contribution is an inference architecture that separates measurement validity, species-level trait distributions, trait geometry, spatial organization, matched nulls, and prospective confirmation.

## Relationship to FCP

The active New Phytologist manuscript remains an ecological application. `disttrait` is the reusable methods layer.

Flower-colour-specific objects such as the four frozen colour states, ROI-v4/EfficientSAM measurement, and the white-versus-nonwhite `q_white` target remain outside the generic package.
