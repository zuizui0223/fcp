# A bounded classification-error sensitivity for stand-level flower-colour differentiation

## Fixed observational object

Maunder et al. (2026), DOI `10.1093/aob/mcag035`, supply all-stand annual mean white/pink/purple shares and mean stand-level Shannon entropy. The existing 2026-10-12 PR #173 analysis already computed the following **uncorrected plug-in observational** Shannon partition with **equal stand weights**:

\[B_Z=H\left(\frac1m\sum_{g=1}^m q_g\right)-\frac1m\sum_{g=1}^m H(q_g),\]

where `q_g` is the distribution of the three *recorded* colours among the sampled plants in stand `g`; the individual `q_g` were **not** reconstructed or recovered. From the known printed summary, `B_Z` has a rounding-conservative lower bound in the parent frozen `result.json` (SHA256 `980b50bfe99f69b31b57d273a2995c814e280cf552a52a1eefaf17e828b4faf1`). This analysis keeps those bounds unchanged.

Define the random experiment: choose a stand `G` with probability `1/m`; choose one of the sampled plants in that stand uniformly; denote its **correct** recorded-at-sampling-time flower colour as `X`, and the observed human-coded colour as `Z`. Thus `B_Z=I(G;Z)`. This is an exact identity for these **sample frequencies** and group weights; it does not identify the composition of all reproductive plants in a stand.

## 1. Arbitrary stand-varying errors, with an upper bound

Let `e=Pr(X != Z)` under the random experiment above: the **equal-stand-weighted fraction of misclassified sampled individuals**. This may be different in each stand and may depend on the actual morph or its environment. It is not a study-reported measured error rate.

By the chain rule and positivity of conditional mutual information,

\[ I(G;Z)\le I(G;X,Z)=I(G;X)+I(G;Z\mid X)\le I(G;X)+H(Z\mid X). \]

Let `E=1{X != Z}`. Conditional entropy obeys

\[H(Z\mid X)\le H(E)+H(Z\mid E,X)\le h_2(e)+e\ln(k-1)\]

where `k=3`, `h_2` is binary entropy in natural logarithms, and the second term is zero for correct labels and at most `ln(k-1)` for an incorrect label. No between-stand independence is needed; classification channels can vary arbitrarily with stand. For a stated **upper bound** `epsilon>=e`, use the monotonically increasing envelope

\[\Phi_k(\epsilon)=\begin{cases}h_2(\epsilon)+\epsilon\ln(k-1), & \epsilon\le 1-1/k,\\ \ln k,&\epsilon>1-1/k.\end{cases}\]

Therefore

\[\boxed{I(G;X) \ge \max\{0,B_Z-\Phi_3(\epsilon)\}.}\]

Since the printed source has rounding error, the implementation conservatively substitutes the source-frozen **lower** `B_Z` bound. The error threshold is the solution to `Phi_3(epsilon*) = B_lower` (approximately 5.42%, 4.78%, 4.58% for 2021–2023). Below that hypothetical error rate, a positive Shannon differentiation component **among correctly classified sampled plant labels** is guaranteed by this deterministic inequality. This is a **necessary** amount of potential mislabelling if equal *sample* true-colour compositions are to account for all observed differentiation, not an exact minimum or an estimated actual misclassification rate.

For the true pooled entropy `H(X)`, conditional entropy similarly gives `H(X)<=H(Z)+Phi_3(epsilon)`. Hence a conservative lower bound on the latent sampled-label between fraction is `(B_lower - Phi)_+ / min(ln 3, Hgamma_original_upper+Phi)`. This is a **deterministic sensitivity bound**, not a frequentist confidence interval. If one cannot certify `epsilon`, it is not a lower-confidence bound on biological differentiation.

## 2. A single common, stand-invariant classification channel

If for each true colour `x`, the recording channel `Pr(Z=z|X=x,G=g)` equals `Pr(Z=z|X=x)` in every stand, then `G -> X -> Z` is a Markov chain and the data-processing inequality guarantees `I(G;Z)<=I(G;X)` **for the channel-induced joint distributions**, regardless of the channel's error rate. But this does **not** mean an empirical plug-in estimate from finite random observations can never be upward biased: finite observation subsampling has to be handled separately.

## 3. What is deliberately not combined

Previous sibling PR #174 bounded **observation-subsampling-only** explanations for the all-stand plug-in Shannon gap under a common *finite census* colour composition and reported census/random-sampling design. It was a probability upper bound under assumptions about stand sampling. This new label-error bound instead relates recorded sampled-plant categories to their unobserved correct categories. **Do not add the separate p-values, multiply their supports, or announce a joint misclassification-plus-subsampling confidence level.** An actual joint probability model would require more information about conditional misclassification, sample sizes, stand-specific field compliance and observer error.

The 2022 inconsistency recorded in #175 involves the **polymorphic-stand conditional summaries**. It does not invalidate the all-stand values used in the frozen #173 Shannon partition; neither is it repaired in this analysis.

## 4. Boundaries

- Colour classification in the original paper used reference photographs but did **not provide a validated misclassification rate or stand-independent confusion matrix** in the recovered main article.
- Source article is one study of one species with 132/112/116 annual stands; the same 105 matched stands are not the basis of the equal-stand-year summary.
- `X` describes **sampled plants**. The underlying true *census* composition `P_g` may differ from the sample even with perfect classification. No census-error-corrected effect size is produced.
- Do not equate pooled between-stand entropy with neutral genetic `F_ST`, heritable morph segregation, selection, evolutionary adaptation or fitness effects.
- The error thresholds are hypothetical **posthoc sensitivity** diagnostics, not measured effect sizes, p-values, confidence limits or new biological replication.

## Sources

1. Maunder et al. (2026), *Annals of Botany*, DOI `10.1093/aob/mcag035`, Methods and Table 2: https://academic.oup.com/aob/article/137/5/1290/8488381
2. PR #173 parent source and Shannon partition: https://github.com/zuizui0223/fcp/pull/173
3. Jensen–Shannon = stand-colour mutual information is a standard identity; see https://www.stat.ethz.ch/CRAN/web/packages/philentropy/vignettes/Information_Theory.html
4. Conditional-entropy/Fano inequality and data processing are standard information-theory results. The self-contained proof above is the operative source; it is **not presented as a new theorem**.
