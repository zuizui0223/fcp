# A conditional bound on observation-subsampling-only heterogeneity

## Scope and notation

This is an application/combination of established probability inequalities, not a claim to a new theorem or a new biological discovery. For a fixed collection of stands, P_g is the true finite census composition of reproductive plants, Q_g its recorded sample composition, and fixed weights w_g sum to one. The target null is P_g=P for every g. This null concerns the actual finite stands. It is **not** a demographic model where small stands are random realizations of a shared gene pool.

A correctly classified full census gives Q_g=P_g exactly; no observation-subsampling uncertainty is added to it. For incompletely censused stands we assume a simple random sample without replacement with a fixed, or colour-independent, sample size n_g>=n0. Sampling designs may be dependent between stands. No random selection of stands from a wider flora is assumed or claimed.

## 1. One-stand exponential bound

Agrawal (2020), Theorem I.3, gives for k-category independent sampling with replacement:

E exp[t D(Q_g||P)] <= (1-t/n_g)^[-(k-1)], 0<=t<n_g.

D is KL divergence with natural logarithms. By the convex comparison in Gross & Nesme (2010), Section D, sampling without replacement has no larger expected value of a convex function of the sample sum. The function x -> exp[t D(x/n_g||P)] is convex for t>=0. Thus the same inequality holds for a simple random sample of finite census labels without replacement. Treat equal-coloured individuals as distinct indexed items, with their labels given by a linear map to category-indicator vectors. Categories with P_k=0 cannot occur in a genuine sample and are omitted; using the original k only weakens the bound. A full census has D=0 and satisfies the bound for every t.

Consequently each stand satisfies E exp[t D_g] <= (1-t/n0)^[-(k-1)] for 0<=t<n0, whether it is a census or an eligible subsample.

## 2. Fixed-weight aggregation without independence

For the observed Jensen-Shannon gap B=H(sum_g w_g Q_g)-sum_g w_g H(Q_g), direct expansion yields

sum_g w_g D(Q_g||P) = B + D(sum_g w_g Q_g||P).

The final term is nonnegative, so B<=sum_g w_g D_g. Applying convexity of exp **before** taking expectation gives

E exp(t B) <= E exp(t sum_g w_g D_g)
             <= sum_g w_g E exp(t D_g)
             <= (1-t/n0)^[-(k-1)].

This uses no between-stand independence. In particular, there is no multiplication by the number of stands in the exponential rate. Shared sites and temporal autocorrelation are not converted into independent replicates. It is invalid to obtain a stronger rate by multiplying marginal generating functions unless independence is additionally justified; we do not do so.

## 3. Tail probability and unknown common composition

Exponential Markov gives P(B>=b) <= inf_{0<=t<n0} exp(-t b)(1-t/n0)^[-(k-1)].

Let a=k-1. If n0*b<=a, return 1. Otherwise the minimizing t is n0-a/b and

U(b;n0,k) = exp[-n0*b+a*(1+ln(n0*b/a))].

This is a conservative analytic p-value upper bound under the fixed null/design, not the exact tail distribution. It holds uniformly in P, so the estimated overall mean is not plugged in as though it were the known null composition. No colour counts or population histories are reconstructed. The expectation counterpart is E B <= ln[1+(k-1)/n0]; an expectation bound alone is not used as a significance test.

## 4. Rounded published input

We reuse the parent three all-stand annual summaries. For each printed frequency, allow +/-0.0005, maintaining nonnegativity and sum one; allow +/-0.005 for the printed average entropy. Concave entropy reaches its minimum at a vertex of this box/simplex polytope. Set b to this minimum total entropy minus the maximum allowed average within entropy. Because the actual unrounded B is at least b, evaluating the decreasing U at b is conservative.

We subtract a further 1e-12 nats from b and multiply computed U by 1+1e-12 (capped at 1) as floating-point margins. The independent Decimal verifier reconstructs the vertices and unpadded bound and checks this numerical conservatism.

Three yearwise bounds receive Bonferroni factor three. This does not require independent years. It only addresses these three declared comparisons, not the entire project-wide history of exploratory analyses. The present source-exposed calculation is supplementary, not a preregistered confirmatory test.

## 5. Application and boundaries

The paper reports whole-stand counting below 100 plants and random samples of 100 (sometimes more) for larger stands. The design therefore implies census-or-n>=100, **not** n>=100 for every stand, and not n=96 for every stand. The global range 1-423 is not a reason to model a one-plant census as a one-plant subsample of a larger stand.

The result can rule against pure observation-subsampling under the reported design. It cannot exclude drift/founder effects, spatially varying recruitment, finite demographic formation, gene flow, plasticity, or selection. Biological randomness can create real differences in finite stand composition. Rejecting equal finite compositions is not rejecting demographic stochasticity or neutral evolution.

Nor does it provide a bias-corrected effect size, a confidence interval on true between-stand diversity, a verified temporal trend, or a test of association with geographic distance/urbanization. The old 21.6-26.3% shares remain uncorrected descriptive values. Unknown individual design compliance, missed plants, systematic colour misclassification, convenience sampling, correlated plant selection or colour-dependent stopping can invalidate the observation model. The n0=50/25 scenarios remain simple-random-sampling assumptions, not estimated effective sizes that correct arbitrary dependence.

## Sources inspected

1. Maunder et al. (2026). DOI 10.1093/aob/mcag035. Geographical survey of flower colour variation; Table 2 All stands; Results Neutral predictions 1 & 2. https://academic.oup.com/aob/article/137/5/1290/8488381
2. Agrawal (2020). Finite-Sample Concentration of the Multinomial in Relative Entropy. Theorems I.2-I.3, PDF p.2, equations visually checked. DOI 10.1109/TIT.2020.2996134. https://arxiv.org/pdf/1904.02291
3. Gross & Nesme (2010). Note on sampling without replacing from a finite collection of matrices. Section D, Lemma 2 and convex-comparison derivation, pp.2-3, extracted text inspected; screenshot service failed for these pages. https://arxiv.org/pdf/1001.2738
