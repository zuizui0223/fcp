# Pre-result protocol: adversarial geographically biased flower-colour photo-label error (2026-10-10)

## Motivation and explicit scope

The frozen New Phytologist FCP manuscript reports statistically repeated geographic allocation of **photographed, visible** four-state within-species flower-colour diversity across three disjoint species cohorts. A separate technical photo-success negative-control audit (PR #149) did not find the same technical classification clustering in all three. That test **cannot** rule out geographically correlated *colour-specific* misclassification.

This strictly post-outcome audit asks a different counterfactual: **How many strategically chosen swaps of two existing same-species photo colour labels would suffice to abolish the positive mean local colour depletion, if errors could be adversarially aligned with geography?** It does NOT assert any such errors exist in real photos.

## Source and frozen baseline

Recover exact existing high-depth colour-measurement rows and enforce SHA256:

- Discovery and validation commit: 5142f7951af0dde5364bb047a566d67e8c479e51
- Third commit: 7e538e5c51c05a7cc47b2fcf53eea92634c8a863
- Four biological photo label categories only: white, yellow-orange, red-pink, blue-purple.
- Retain species with >=40 classifiable photo records; evaluate those with >=30 original same-species within-50km photograph pairs.
- Baseline MUST match source manuscript exploratory receipt exactly: n_species discovery/validation/third = 166/181/204; mean absolute local colour discordance depletion = 0.020529254583812922 / 0.018672971642749295 / 0.01468491968437185 (tolerance 1e-10).

## Attack and algorithm

At each step identify the **single best** within-species pair of different coarse-colour classified photo labels for exchange that increases pairwise local colour discordance by the greatest number of local edges per eligible species weight. All source photo sites and between-photo geographic pair graphs remain frozen. Swapping the labels of two records in the **same species** holds its exact sampled global four-state colour proportions and Gini-Simpson pair diversity unchanged; however it changes the *geographic allocation* of the existing colours.

The local-edge change of swapping photo i with label a and photo j with b is exactly

    C_i(a) - C_i(b) + C_j(b) - C_j(a),

where C_i(k) is the number of <=50km same-species neighbours of photo i bearing label k. The i-j edge, if present, cancels from this formula. Verify exact recomputation after every accepted swap. Only strictly positive gains are allowed. Do not reuse any altered source photo in a later swap. Across all species choose the largest reduction in the **equal-species mean** depletion; tie-break deterministically. Stop at the first mean <=0 or 10% of originally classifiable photographs belonging to species with sufficient local pairs, whichever occurs first.

Report the baseline, all realized constructive pair swaps, full monotone attack trajectory, budgets 0/0.1/0.25/0.5/1/2/5/10%, and first realized nonpositive mean if any. Edited photo fraction is 2*number_of_accepted_swaps / number_of_classifiable_photos_in_50km_evaluable_species; the same photo cannot be counted twice. All labels are hypothetical; primary source data remain unchanged.

## Strong inference firewall

- This is a **constructive adversarial perturbation**, not a statistical null, a distribution-free lower bound, a mathematically optimal minimum error rate, a measured misclassification rate, or a calibrated observational robustness guarantee.
- A small edit fraction for the worst-case geography-aware adversary does NOT imply the real classifier makes that many errors. Such an adversary is granted the observed spatial colour outcome and can target unusually influential photos; photo-realistic error mechanisms generally cannot target the spatial statistic in this way.
- A large fraction would only exclude the specifically constructed pair-swap attack up to its permitted search budget; it would not rule out other errors, including morph-dependent omission, correlated illumination, whole-site relabeling, species mix-ups or within-plant colour plasticity.
- Fractions cannot be compared directly to PR #149 binary-classification geographical depletion or treated as the fraction of signal attributable to nonclassifiable photos.
- The resulting post-attack label configurations have been chosen by optimizing on the outcome. Do NOT attach the original permutation p-values to them and do NOT claim a valid new post-attack unadjusted significance test.
- A fresh photo-based botanical adjudication with blinded annotators and ideally calibrated flower-colour standards remains required to estimate real white/nonwhite and hue-specific confusion, and individual-linked genotype/fitness evidence is necessary to claim inherited polymorphism or selection.
- Source H1/H2, existing full-paper claims, 42,111 global metadata opportunity frame, and unopened future 2,000+730 taxon cohorts do not change.

**Status before CI completes:** algorithm and synthetic controls, no quantitative biological outcome yet.
