# Hesperis: bounded observation-subsampling audit

2026-10-12. Retrospective, already exposed to source summaries and their entropy decomposition. This is neither preregistration nor a new biological sample. Fixed before evaluating the new probability-bound outputs.

Goal: decide whether the published equal-stand Shannon gap can be attributed entirely to within-stand random subsampling. Keep the parent input, full-year cohorts, colour labels, and primary entropy output unchanged. Do not infer a time trend or causes of differentiation.

Target: equality of the **finite census compositions** of all stands within each year. This is NOT a test of a shared superpopulation with demographic drift, selection, gene flow, or an urbanization slope. A census has zero observation-subsampling error even if it contains one individual; demographic randomness is not set to zero.

Reported design: stands smaller than 100 are censused; larger stands receive a random sample of 100, with some larger samples. Treat every stand as a census OR a simple random sample without replacement of at least n0=100. This is conditional on the published protocol; individual compliance cannot be audited without raw counts.

For each year, use the smallest Jensen gap b compatible with published rounding (frequency +/-0.0005 and mean Shannon entropy +/-0.005; simplex preserved). Subtract 1e-12 numerically to make the input conservative.

Under identical finite census compositions, D_g=KL(q_g||p) is zero for census groups, otherwise its MGF is at most (1-t/n0)^-(k-1), 0<=t<n0. Gross-Nesme/Hoeffding convex comparison extends Agrawal's multinomial bound to sampling without replacement. The sample Jensen gap B <= sum(w_g D_g). Convexity of exp yields E exp(t B) <= sum w_g E exp(t D_g), **without assuming between-stand independence**. Chernoff minimization gives U(b)=1 for n0*b<=k-1, otherwise exp[-n0*b+(k-1)*(1+ln(n0*b/(k-1)))]. Use k=3, natural logarithms. Uniform in p; do not plug the observed mean in as a known null probability.

Primary n0=100. Prespecified diagnostic scenarios n0 in [1,25,50,75,100]. These are alternative simple-random-sampling minimums, NOT plug-in effective sample sizes or empirically measured design effects. Bonferroni over the three yearwise upper bounds; not independence of years. Invert the same bound to report the minimum n0 sufficient for all three yearwise Bonferroni bounds <=0.05. All computations are source-exposed and supplementary.

Verification: fail tests first; formula boundary/invalid input checks; exact small multinomial and multivariate-hypergeometric enumerations; perfectly dependent groups; census-vs-subsample guards; direct KL identity. Recompute rounding minimum and probabilities with an independent 60-digit Decimal implementation that does not import the analysis. Test all declared conditions and retain both supportive and non-supportive scenarios. Prior parent tests also rerun. Full repository suite, if available, reported separately.

Caveats: no corrected effect size or interval for true diversity; no phenotype misclassification correction; no claim of mechanistic selection, genotypes, individual plant independence beyond actual probability sampling, population independence, or multi-species generality. Unequal/proximity/colour-dependent sampling, missed plants, incomplete census, or protocol exceptions can invalidate this observation-only bound. Global research multiplicity is not resolved by adjusting the three new yearwise comparisons.

References: Maunder et al. 2026, DOI 10.1093/aob/mcag035, Methods geographical survey / Table 2 / Results Neutral predictions 1 & 2. Agrawal 2020, DOI 10.1109/TIT.2020.2996134, Theorems I.2-I.3. Gross & Nesme 2010, arXiv:1001.2738, Section D, Lemma 2 and following convex comparison.
