# H1–H3 statistical reference audit

Checked 2026-09-15. Scope: four requested primary papers, bibliographic identity and bounded methodological support. No FCP data, models, images, frozen outcomes, or opening permissions were changed. This is not a complete statistical bibliography or an independent validation of the FCP implementation.

## Lin: concordance, not rank agreement

Lin, L. I.-K. (1989). A concordance correlation coefficient to evaluate reproducibility. *Biometrics*, **45**(1), 255–268. https://doi.org/10.2307/2532051

Identity: [journal issue record](https://www.jstor.org/stable/i343420) gives Lawrence I-Kuei Lin, title, pages and DOI; [primary abstract deposited in PubMed](https://pubmed.ncbi.nlm.nih.gov/2720055/) confirms year, journal and pages. The abstract defines agreement relative to the identity line.

Application: cite for H1 CCC as identity-line concordance, distinct from rank-order repeatability. It does not establish that repeated image measurements are biologically accurate, nor justify FCP thresholds, observer independence, rarefaction or its bootstrap procedure. Those require design/code evidence.

Access limit: bibliographic record and abstract inspected; the [JSTOR article endpoint](https://www.jstor.org/stable/2532051) did not expose the full article. Original formula derivation and interval methods were not full-text audited here. Do not describe this as a full-text verification of the estimator implementation.

## Phipson and Smyth: nonzero randomization probabilities

Phipson, B., & Smyth, G. K. (2010). Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn. *Statistical Applications in Genetics and Molecular Biology*, **9**(1), Article 39. https://doi.org/10.2202/1544-6115.1585

Access: [author-hosted corrected paper](https://gksmyth.github.io/pubs/PermPValuesPreprint.pdf), 12 pages, marked published 31 October 2010 and corrected 9 February 2011; sections 4–6 inspected. [Bibliographic record](https://pubmed.ncbi.nlm.nih.gov/21044043/) corroborates authors, title and DOI.

Application: supports reporting `(b + 1)/(B + 1)` rather than zero for randomization exceedances. Section 4 derives the Monte Carlo expression; section 5 treats sampling without replacement. Section 6 distinguishes with-replacement sampling, where this expression is valid but can be conservative, not universally the exact exhaustive tail probability.

Limit: this numerical correction presupposes a valid null generation mechanism. It does not repair selection-conditioned inference, missing exchangeability, observer/season dependence, or post-hoc hypothesis choice. Cite it for probability computation, not as validation of FCP's null model. No original FCP draws were audited in this reference check.

## Blomberg, Garland and Ives: phylogenetic pattern, not cause

Blomberg, S. P., Garland, T., Jr., & Ives, A. R. (2003). Testing for phylogenetic signal in comparative data: behavioral traits are more labile. *Evolution*, **57**(4), 717–745. https://doi.org/10.1111/j.0014-3820.2003.tb00285.x

Access: [original article PDF hosted at the University of Wisconsin](https://pages.cs.wisc.edu/~larget/botany940/Blomberg_etal2003.pdf), 29 pages; introduction (pp. 717–718) and K definition (pp. 722–723, equation 4) inspected. [PubMed bibliographic record](https://pubmed.ncbi.nlm.nih.gov/12778543/) verifies the DOI and author suffix. Publisher access attempts failed; the PDF is the original paper, not a secondary summary.

Application: K compares observed covariance pattern with Brownian-motion expectations for a specified topology and branch lengths. The paper explicitly separates phylogenetic resemblance from its generating mechanism and notes sensitivity to trait/tree error.

Limit: neither K nor its significance establishes inheritance, adaptation, ecological causation, or absence of phylogenetic structure after a nonsignificant result. This paper is not the original reference for Pagel's lambda, rank-PGLS, or FCP-specific multiple-testing corrections; those remain separate citation/implementation checks.

## Jin and Qian: assembled plant phylogenies

Jin, Y., & Qian, H. (2022). V.PhyloMaker2: An updated and enlarged R package that can generate very large phylogenies for vascular plants. *Plant Diversity*, **44**(4), 335–339. https://doi.org/10.1016/j.pld.2022.05.005

Access: [PubMed record](https://pubmed.ncbi.nlm.nih.gov/35967255/) confirms authors Yi Jin and Hong Qian, title, volume, pages and DOI. Indexed text from the [deposited original article](https://pmc.ncbi.nlm.nih.gov/articles/PMC9363651/) exposed its abstract, Box 1 and species-list preparation passage; direct opening encountered a browser challenge. The [author's supplementary repository](https://github.com/jinyizju/Supplementary-files-of-the-V.PhyloMaker2-paper) independently links the original article, but its supplementary files were not read.

Application: supports identifying the named plant-tree construction package and its use of a supplied species list, backbone and placement scenario. The exposed species-list passage describes optional relative-based attachments.

Limit (inference from that construction): assembled placements are not new molecular estimates for every added tip. Running S1–S3 does not itself supply uncertainty-calibrated posterior phylogenies. The exact scenario algorithms, FCP backbone version, taxonomy joins, unmatched taxa and downstream covariance require separate code/artifact evidence. Do not claim the full article or supplementary scenario definitions were inspected here.
