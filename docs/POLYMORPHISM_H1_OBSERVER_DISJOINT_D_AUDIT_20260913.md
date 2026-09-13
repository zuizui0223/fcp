# H1 observer-disjoint D audit — execution note

Date: 2026-09-13 JST
Branch: `analysis/polymorphism-42111-h1-h2-gates-20260912`

## Purpose

Determine whether the previously reported observer-disjoint reproducibility value (~0.971) is a direct replication of the current four-state polymorphism score

`D = 1 - sum_k p_k^2`,

computed independently from observer-disjoint subsets within species, or whether it came from a different colour/measurement representation.

## Decision rule

1. If a repository result already computes four-state D independently in observer-disjoint subsets, recover the exact implementation/result and use it as the H1 measurement-validity evidence.
2. If not, do not reinterpret another reliability statistic as D reliability. Freeze and run a dedicated observer-disjoint D analysis instead.
3. The dedicated analysis must keep the four biological states fixed (`white`, `yellow_orange`, `red_pink`, `blue_purple`), exclude `mixed_uncertain` from the D numerator/denominator, and require sufficient classifiable observations in both observer-disjoint halves.
4. H2/H3 outcomes must not be used to tune the H1 split, threshold, or statistic.

This note is only an audit/execution marker; it does not itself establish H1 support.
