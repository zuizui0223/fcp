# H3 phylogenetic method reference audit

Checked 2026-09-15. Scope: bibliographic identity and bounded support for Pagel's
lambda and phylogenetic linear regression. No measurements, analyses, frozen
outcomes, software installations or opening permissions were changed. This note
does not certify the historical runtime or validate the FCP implementation.

## Pagel: covariance pattern, not an identified mechanism

Pagel, M. (1999). Inferring the historical patterns of biological evolution.
*Nature*, **401**, 877–884. https://doi.org/10.1038/44766

Identity and access: the [publisher record](https://www.nature.com/articles/44766)
confirms title, author, publication date (28 October 1999), pages and DOI. Only
the abstract and bibliographic information were accessible. The publisher PDF
redirect was unavailable; a university-hosted copy also failed to open. The
original mathematical derivation was therefore not full-text audited here.

For the operational covariance definition, the official
[ape `corPagel` documentation](https://stat.ethz.ch/CRAN/web/packages/ape/refman/ape.html#corPagel)
states that lambda multiplies Brownian-motion off-diagonal covariances while
leaving variances unchanged. It cites Pagel (1999). This supports describing a
fitted covariance transformation, not scaling every branch without qualification.
The ape function is a definition reference here, not a claim that FCP used ape
to fit H3 lambda.

Interpretation boundary (inference from this model definition): a fitted lambda
describes compatibility of the measured trait with a tree-derived covariance
family. It is not by itself evidence identifying inheritance, adaptation,
pollinator causation or another generating mechanism. A nonsignificant test is
not equivalence to no phylogenetic structure. In regression, the covariance
concerns variation left after the specified predictors; it must not be conflated
with an unadjusted trait-signal estimate.

## Ho and Ane: likelihood computation and phylogenetic regression

Ho, L. S. T., & Ané, C. (2014). A linear-time algorithm for Gaussian and
non-Gaussian trait evolution models. *Systematic Biology*, **63**(3), 397–408.
https://doi.org/10.1093/sysbio/syu005

The [author's publication list](https://pages.stat.wisc.edu/~ane/publis.html)
confirms title, journal, pages and associated phylolm software. The
[deposited primary abstract and bibliography](https://pubmed.ncbi.nlm.nih.gov/24500037/)
confirm authors and DOI and describe efficient likelihood/parameter calculations
for tree-structured covariance models, including Pagel's lambda. This supports
citing the computational method and software, not a universal calibration
guarantee for rank-based FCP inference. The publisher full-text redirect failed;
proofs and simulations were not full-text audited.

The [official phylolm manual](https://cran.r-project.org/web/packages/phylolm/phylolm.pdf),
version 2.6.5, function `phylolm`, pp. 14–16, documents lambda as an error-covariance
model, ML as the default, and coefficient covariance conditional on the fitted
phylogenetic parameter. It lists defaults `boot=0`, `measurement_error=FALSE`,
and lambda bounds `[1e-7, 1]`. It also warns that lambda and an added measurement
error component are not jointly identifiable in this implementation. These are
documentation facts, not verified FCP run settings or grounds for changing a
frozen analysis. A separate image-measurement validation problem is not repaired
by fitting lambda.

## H3a and H3b are different statistical operations

The official [phytools `phylosig` manual](https://cran.r-project.org/web/packages/phytools/phytools.pdf),
pp. 184–185, distinguishes lambda likelihood-ratio testing against lambda zero
from K randomization. For lambda it returns `lambda`, `logL`, `logL0` and `P`;
`nsim` belongs to the K randomization. It describes optimization between zero
and a tree-specific upper limit. Do not transplant phylolm's documented bounds
to phytools or label the lambda likelihood-ratio P value a permutation P value.

Direct local-source verification: H3a uses
`phytools::phylosig(method="lambda", test=TRUE)`, whereas H3b fits
`phylolm::phylolm(y_rank ~ span_rank + nclass_rank + nobs_rank, model="lambda")`
on tree-matched, centered average ranks. H3b
explicitly computes `2 * pnorm(-abs(beta / se))`. That is a two-sided normal
approximation using the fitted coefficient standard error, not an exact
permutation test or automatically the package summary test. This distinction
comes from the FCP implementation, not from attributing that custom calculation
to Ho and Ané. The inspected H3a runner's `lambda_test` is at lines 111–119,
Git blob `39b4d653ef97106474c9a106316ffffda9f5b54c`; the H3b runner's
`rank_pgls` is at lines 100–131, Git blob
`23aecf021bcde80d97ca6ded7926373adc5464b0`. Both belong to publication
commit `6cfe52cc373ee65e0309d482945b353e998364b0` and were read, not rerun.

Retained local original-result copies report R 4.6.1, ape 5.8.1 and jsonlite
2.0.0 for both runs; H3a reports phytools 2.5.2 and H3b phylolm 2.6.5.
The inspected result-file SHA256 values are
`73edde717c094cbe64295ef989ed893086ad3e7679a720bab5b88949d56a6c95`
(H3a) and `d7fadae50a3e12794978f5ff6eea3e08bdbe7a0e485a49f48b44bdd5f9ea8f79`
(H3b). These are recorded runtime versions, not a newly reproduced environment
or proof of all optimizer defaults. Original artifact IDs and archive digests
are listed in Supplement S3.

Current online manuals do not establish which package versions were installed
in the historical Actions runs. Archived session information or exact package
artifacts are required before asserting historical defaults or dependency pins.
Neither these references nor a successful optimizer removes uncertainty about
tree placement, measurement validity, sampling, model adequacy or missing taxa.
