# Hesperis: one joint observation-error bound (sampling and flower-colour coding)

**Status:** retrospective mathematical sensitivity extension, 2026-10-12. Source-observed error rates are unavailable. This is an application of existing information-theoretic / concentration inequalities, **not a new probability theorem, new field measurement or evolutionary confirmation**. The analytical derivation below is the operative proof; no test statistic has been fitted to unrecovered stand-by-year records.

## 1. What was actually observed

Maunder et al. (2026, DOI 10.1093/aob/mcag035), Table 2 and the reported yearly mean Shannon `H′`, provide the equal-stand mean of recorded morph frequencies (white, pink, purple) and the equal-stand mean of their entropies. Previous PR173 recovered, under ordinary published rounding, a **lower bound** on the observed stand–morph Shannon Jensen gap for each year:

| Year | Number of stands | Observed reported `B_Z` | Lower rounded bound `b` (nats) |
|---|---:|---:|---:|
| 2021 | 132 | 0.2538601243 | 0.24814489869438183 |
| 2022 | 112 | 0.2305232935 | 0.22496448574712613 |
| 2023 | 116 | 0.2232261736 | 0.21776296903477543 |

All `b` are read directly from an unmodified, SHA-256-locked PR173 JSON (`980b50bfe99f69b31b57d273a2995c814e280cf552a52a1eefaf17e828b4faf1`); this extension does not refit the source or construct hidden stand counts. The 2022 inconsistency noted in PR175 is **for conditional polymorphic-stand means**, not the all-stand summaries used here.

## 2. Three statistical layers are not interchangeable

- `P_g` is the **true finite census composition** of all flowering/reproductive plants in stand `g` (unobserved in the public source for larger stands).
- `Q_g` is the empirical composition among the **randomly sampled plants when all their morphs are classified correctly**; below, sampled true category `X`.
- `R_g` is the observed empirical composition after possible morph classification errors; below, coded category `Z`.

The stand index `G` is drawn uniformly across the year's **fixed observed stands**, and one of its sampled plants is then drawn uniformly. Hence, exactly, `B_X=I(G;X)` and `B_Z=I(G;Z)`, where entropy uses natural logarithms. These are sample-level Jensen gaps; neither is genetic Fst. Stand set selection is conditioned on and not modeled as random geographical population sampling.

## 3. A deterministic coupling for arbitrary stand-varying category errors

Let `e = Pr(X != Z)` under the two-stage stand-uniform experiment above. This is the **realized equal-stand-weighted sampled-plant misclassification fraction**, not a source-measured value. The category error process may depend on morph, stand, observer, season and sampling composition. The chain rule and binary error indicator yield

```
 B_Z = I(G;Z) <= I(G;X)+H(Z|X)
               <= B_X+h2(e)+e ln(k-1),    k=3.
```

Define the increasing ceiling `Phi_k(eps)=h2(eps)+eps ln(k-1)` for `eps<=1-1/k`, and `ln k` otherwise. Thus on **any realized observation** with `e<=eps`, one has `B_X >= B_Z-Phi_k(eps)`. There is no independence requirement between stand membership and misclassification.

This is an information-theoretic error accounting identity (conditional entropy/Fano-type bound), *not* an assumption of a known common observer confusion matrix.

## 4. The sampling-only null bound, reused without multiplying p-values

The deliberately **narrow** null states that every stand's **actual finite population** has exactly the same colour composition, `P_g=P`, where `P` may be any categorical vector. Under the source-reported design, stands with fewer than 100 reproductive plants are completely counted, while in larger stands a colour-independent simple-random sample *without replacement* of at least `n0` plants is assumed. A complete census adds no observation-subsampling noise even if it has just one plant. No particular number of independent stands is assumed.

Agrawal's gamma-MGF bound for a `k`-category empirical KL divergence, extended to finite sampling without replacement by Hoeffding's convex-order comparison, implies for each incompletely sampled stand:

```
 E[ exp(t KL(Q_g || P)) ] <= (1-t/n0)^(-(k-1)),  0 <= t < n0.
```

For a complete census, `KL(Q_g||P)=0` under this null. The identity

```
 sum_g w_g KL(Q_g||P) = B_X + KL(sum_g w_g Q_g || P) >= B_X
```

and convexity of the exponential **before** expectation give

```
 E exp(t B_X) <= sum_g w_g E exp(t KL(Q_g||P))
                <= (1-t/n0)^(-(k-1)).
```

This deliberately does **not** multiply moment-generating functions across stands, so no between-stand statistical independence is needed. Chernoff minimization leads to the sampling-null one-year tail ceiling

```
 U(b; n0,k) = 1                                  if n0*b <= k-1 or b <= 0,
            = exp[-n0*b+(k-1)*(1+ln(n0*b/(k-1)))] otherwise.
```

This repeats the earlier PR174 approach on the *same all-stand source* and is an **upper bound on a tail probability**, not the exact null p-value or a corrected diversity estimator. It assumes the actual field subsampling followed the stated design. If stands were preferentially sampled by colour or actual sample sizes were below `n0` without being censuses, the bound may not apply.

## 5. A valid combined statement (without independence multiplication)

On any outcome meeting `B_Z >= b` and `e<=eps`, Step 3 implies `B_X>= b-Phi_3(eps)`. Therefore the two-event inclusion gives, under the sampling null,

```
 Pr_null( B_Z >= b  AND  e <= eps )
 <= Pr_null( B_X >= b - Phi_3(eps) )
 <= U(b-Phi_3(eps); n0, 3).
```

If the classification error budget holds **pathwise for every possible observed sample**, then `Pr(e>eps)=0` and this is also an *unconditional* upper bound on the observed-statistic tail. More generally, if independently justified error-model validation can certify `Pr(e>eps)<=eta` for the observation process, a union bound gives

```
 Pr_null(B_Z >= b) <= eta + U(b-Phi_3(eps);n0,3).
```

**Without either type of error budget, only the joint-event bound exists**, not an unconditional p-value. The reported source does not certify any such budget. If a future audit determines `e<=eps` in the *realized field sample*, its observed event can be tested against the *joint-event probability* above; one must **not** call this `Pr(B_Z>=b | e<=eps)`.

For the three annual checks, use `min(1,3[eta+U])` by Bonferroni. Independence across years is not needed. The adjustment protects only this specific three-year family, not the many exploratory analyses performed throughout the project. No combined 95% effect-size confidence interval is produced.

## 6. Bounded scenario and what its sign means

`n0=100` means complete census for small stands, and a colour-independent random sample of **at least 100** from larger ones; it does **not** say every stand was sampled 100 times. The main hypothetical errors are equal-stand-weighted `eps=0%, 0.25%, 0.5%, 1%, 1.5%, 2%, 2.5%, 3%, 4%, 5%`. The effective minimum `n0=75,50,25` are **alternative designs**, not estimates of effective sample size or a cure for unknown colour-dependent observation bias. No cut-off has been selected to maximize publication significance: the complete grid, including failures, is archived.

Under `n0=100, eps=2%, eta=0`, the conservative three-year-Bonferroni upper tail is:

- 2021: 0.001245449651,
- 2022: 0.008710597630,
- 2023: 0.015690633971.

These values are **conditional on a guaranteed <=2% realized label-error budget**. No actual misclassification estimate is reported. The critical maximum error rates for the 0.05 Bonferroni threshold are 3.031%, 2.487%, and 2.323%, respectively, under n0=100. Under n0=25 even zero error does not cross 0.05 with this coarse inequality; that is inconclusive, not evidence of equality.

For an alternative probabilistic error budget `Pr(e>eps)<=0.01`, adding three times 0.01 produces yearwise adjusted upper bounds `0.031245, 0.038711, 0.045691` for the 2%/100 scenario. This is a **mathematical example** only: no source observation validates eta.

## 7. Biological nonclaims

Rejecting the **common actual finite census composition plus observation-error-only** model does *not* exclude biological demographic stochasticity, genetic drift, multiple horticultural introductions, founder events, migration/gene flow, plasticity, assortative mating, or selection. It cannot prove any of those mechanisms. It does not infer true stand-level morph frequencies, provide a sampling-bias-corrected percentage of total diversity, or test the outcome-selected subset of polymorphic stands from PR175.

We do not combine the sibling PR174 sampling p-upper with the PR176 misclassification lower bound by multiplying p-values, nor claim either result was a separate independent biological replication. The derivation is a **single event-inclusion bound** under jointly stated assumptions. It is retrospective and source-exposed, not registered prospectively.

## Source references

- Maunder et al. (2026), Annals of Botany, DOI 10.1093/aob/mcag035: https://academic.oup.com/aob/article/137/5/1290/8488381
- Agrawal (2020), *Finite-Sample Concentration of the Multinomial in Relative Entropy*: https://arxiv.org/abs/1904.02291
- Gross & Nesme (2010), *Note on sampling without replacing from a finite collection of matrices*, Hoeffding convex-comparison argument: https://arxiv.org/abs/1001.2738
- Earlier executed source-replay proofs: https://github.com/zuizui0223/fcp/pull/174 and https://github.com/zuizui0223/fcp/pull/176
