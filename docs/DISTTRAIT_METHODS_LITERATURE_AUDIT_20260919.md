# disttrait methods-paper literature audit — 2026-09-19

Status: source audit for `docs/DISTTRAIT_METHODS_MANUSCRIPT_V0_1_20260919.md`.

Purpose: map manuscript-level background and methodological claims to
representative literature without using citations to overstate what the
`disttrait` benchmark itself establishes.

## 1. Intraspecific trait variation and distributional trait thinking

### Bolnick et al. 2011

Daniel I. Bolnick et al. "Why intraspecific trait variation matters in
community ecology." Trends in Ecology & Evolution 26:183–192.

DOI: `10.1016/j.tree.2011.01.009`

Use in manuscript:

- supports the ecological importance of phenotypic variation among individuals
  within populations/species;
- motivates retaining within-species distributions rather than only species
  means.

Does not support:

- the specific `disttrait` spatial statistic;
- reliability of opportunistic observations.

### Siefert et al. 2015

Andrew Siefert et al. "A global meta-analysis of the relative extent of
intraspecific trait variation in plant communities." Ecology Letters
18:1406–1419.

DOI: `10.1111/ele.12508`

Use in manuscript:

- demonstrates that intraspecific trait variation can constitute a substantial
  share of trait variation in empirical communities;
- supports the motivation for explicitly representing ITV in comparative
  ecology.

Does not support:

- global universality of any particular ITV fraction outside the study's plant
  community scope.

### Carmona et al. 2016

Carlos P. Carmona, Francesco de Bello, Norman W. H. Mason & Jan Lepš.
"Traits Without Borders: Integrating Functional Diversity Across Scales."
Trends in Ecology & Evolution 31:382–394.

DOI: `10.1016/j.tree.2016.02.003`

Use in manuscript:

- supports distributional/probabilistic representations of traits across
  organismal-to-macroecological scales;
- provides conceptual precedent for treating trait distributions rather than a
  single species value.

### Carmona 2019

Carlos P. Carmona. "Trait probability density (TPD): measuring functional
diversity across scales based on TPD with R." Ecology 100:e02876.

DOI: `10.1002/ecy.2876`

Use in manuscript:

- reinforces explicit inclusion of intraspecific trait variability in a
  probabilistic trait-distribution framework.

The manuscript should not imply that `disttrait` is an implementation of TPD;
the frameworks share a distributional motivation but target different
estimands.

## 2. Opportunistic / citizen-science observation processes

### Bird et al. 2014

Tomas J. Bird et al. "Statistical solutions for error and bias in global
citizen science datasets." Biological Conservation 173:144–154.

DOI: `10.1016/j.biocon.2013.07.037`

Use in manuscript:

- supports the general claim that large citizen-science data contain sampling
  error/bias and that metadata about the sampling process matter;
- motivates separating observation-process validity from downstream inference.

### Di Cecco et al. 2021

Grace J. Di Cecco, Vijay Barve, Michael W. Belitz, Brian J. Stucky,
Robert P. Guralnick & Allen H. Hurlbert. "Observing the Observers: How
Participants Contribute Data to iNaturalist and Implications for Biodiversity
Science." BioScience 71:1179–1188.

DOI: `10.1093/biosci/biab093`

Use in manuscript:

- supports treating observer behaviour as part of the iNaturalist observation
  process;
- motivates observer restrictions and observer-disjoint validation when such
  metadata exist.

### Laitly et al. 2021

Alexandra Laitly, Corey T. Callaghan, Kaspar Delhey & William K. Cornwell.
"Is color data from citizen science photographs reliable for biodiversity
research?" Ecology and Evolution 11:4071–4083.

DOI: `10.1002/ece3.7307`

Use in manuscript:

- application-specific precedent that colour extracted from citizen-science
  photographs requires validation and can contain useful biological signal;
- relevant to the FCP origin of `disttrait`.

Do not use it to claim that image-derived colour is unbiased at the individual
level; the paper reports both utility and additional uncontrolled variation.

## 3. Permutation / Monte Carlo calibration

### Phipson & Smyth 2010

Belinda Phipson & Gordon K. Smyth. "Permutation P-values Should Never Be
Zero: Calculating Exact P-values When Permutations Are Randomly Drawn."
Statistical Applications in Genetics and Molecular Biology 9:Article 39.

DOI: `10.2202/1544-6115.1585`

Use in manuscript:

- supports the plus-one Monte Carlo probability used for randomly drawn
  permutation/null worlds.

Does not validate the exchangeability assumption of a particular permutation
scheme; that remains an application-specific scientific assumption.

## 4. Distance-matrix inference and its limits

### Legendre & Fortin 2010

Pierre Legendre & Marie-Josée Fortin. "Comparison of the Mantel test and
alternative approaches for detecting complex multivariate relationships in the
spatial analysis of genetic data." Molecular Ecology Resources 10:831–844.

DOI: `10.1111/j.1755-0998.2010.02866.x`

### Guillot & Rousset 2013

Gilles Guillot & François Rousset. "Dismantling the Mantel tests."
Methods in Ecology and Evolution 4:336–344.

DOI: `10.1111/2041-210x.12018`

Use in manuscript:

- establishes that distance-matrix correlations/permutations should not be
  treated as a generic substitute for regression or as automatically valid in
  the presence of spatial structure;
- motivates explicit statement of the exact estimand and permutation unit.

Critical boundary:

- `disttrait` is **not** presented as a generic Mantel test;
- these references do **not** prove the validity of the `disttrait` matched
  null;
- `disttrait` conditions within species, fixes observed coordinates and
  complete trait rows, and validates its specific null by simulation under
  explicit DGPs;
- the MNAR benchmark demonstrates a setting in which conditional randomization
  cannot recover the latent biological null.

## 5. Random-effects summaries and finite-sample calibration

### DerSimonian & Laird 1986

Rebecca DerSimonian & Nan Laird. "Meta-analysis in clinical trials."
Controlled Clinical Trials 7:177–188.

DOI: `10.1016/0197-2456(86)90046-2`

Use in manuscript:

- source for the simple random-effects heterogeneity estimator used in the
  species-slope summary.

### Viechtbauer 2010

Wolfgang Viechtbauer. "Conducting Meta-Analyses in R with the metafor
Package." Journal of Statistical Software 36(3):1–48.

DOI: `10.18637/jss.v036.i03`

Use in manuscript:

- general reference for fixed/random-effects meta-analytic models and
  heterogeneity analysis.

### Röver, Knapp & Friede 2015

Christian Röver, Guido Knapp & Tim Friede. "Hartung-Knapp-Sidik-Jonkman
approach and its modification for random-effects meta-analysis with few
studies." BMC Medical Research Methodology 15:99.

DOI: `10.1186/s12874-015-0091-1`

Use in manuscript:

- supports the general caution that random-effects inference can be sensitive
  to few units and heterogeneous precisions;
- provides literature context for checking rather than assuming large-sample
  calibration.

Critical boundary:

The concrete `disttrait` claim that its analytic DL/Q calibration reached
0.175/0.30 null rejection is established by the frozen v0.8 benchmark, not by
these external papers.

## 6. Multivariate trait-space metrics

### Gower 1971

J. C. Gower. "A General Coefficient of Similarity and Some of Its
Properties." Biometrics 27:857–871.

DOI: `10.2307/2528823`

Use in manuscript:

- historical example that multivariate resemblance/distance requires an
  explicit definition and can be constructed differently for different data
  types.

Critical boundary:

`disttrait 0.12` uses Euclidean distance for continuous multivariate traits;
it does not implement Gower distance as its default. The citation is therefore
used to support the broader point that metric choice is substantive, not to
describe the implemented algorithm.

## 7. External empirical sources

### San Francisco street trees

The frozen `disttrait` fixture is a historical snapshot of the San Francisco
Public Works street-tree inventory distributed via TidyTuesday. The municipal
dataset is public and uses the Open Data Commons Public Domain Dedication and
License (PDDL) 1.0.

The current municipal inventory has changed since the frozen 2020 snapshot;
the methods paper must cite the frozen source/provenance and must not present it
as a current tree census.

### ShareTraitDatabase

Irene Martorelli, Brett G. Olivier, Raimon Cuxart-Erruz, Félix P. Leiva,
Matty P. Berg, Jacintha Ellers & Wilco C.E.P. Verberk. ShareTraitDatabase
v1.0.0. Zenodo.

DOI: `10.5281/zenodo.16537297`

Licence: CC BY 4.0.

### Shokri et al. 2022

Milad Shokri, Francesco Cozzoli, Fabio Vignes, Marco Bertoli,
Elisabetta Pizzul & Alberto Basset. "Metabolic rate and climate change across
latitudes: evidence of mass-dependent responses in aquatic amphipods."
Journal of Experimental Biology 225:jeb244842.

DOI: `10.1242/jeb.244842`

Underlying data DOI: `10.17605/OSF.IO/56TNH`.

Use in manuscript:

- source biological study for the external `Gammarus insensibilis` fixture.

The generic `disttrait` null result does not test or overturn the source
study's mass/temperature hypotheses.

## 8. Manuscript claim policy after this audit

The literature should support background, precedent and known methodological
cautions. The package's performance claims must remain anchored to its frozen
benchmarks and empirical receipts.

In particular:

- do not cite Mantel critiques as proof that the package solves spatial
  autocorrelation;
- do not cite random-effects literature as proof that the package-specific
  permutation calibration is valid;
- do not cite citizen-science literature as proof that FCP observation bias is
  absent;
- do not cite ITV reviews as evidence for the paper's simulated FPR/power
  values;
- distinguish source-study biological conclusions from the deliberately generic
  transport estimands.
