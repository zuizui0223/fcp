#!/usr/bin/env python3
"""Robustness analyses for the mixed-preserving C/S climatic analysis.

Analyses:
1. add an outcome-independent v2.2 exact-name literature-attention covariate;
2. restrict to species with exactly one resolved positive-evidence source, equalising the
   positive-source count without pretending that this count is independent observation effort;
3. fit a non-ordinal three-state multinomial model (C-only / S-only / C+S).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

METRICS = [
    'temperature_breadth',
    'moisture_breadth',
    'climatic_heterogeneity',
    'pca_dispersion',
    'pca_hull_area',
]
OUTCOMES = ['C_local_coexistence_documented', 'S_spatial_segregation_documented']
STATE_ORDER = ['local_coexistence_only', 'spatial_segregation_only', 'coexistence_and_segregation']


def zscore(x: pd.Series) -> pd.Series:
    x = pd.to_numeric(x, errors='coerce')
    sd = x.std(ddof=0)
    if not np.isfinite(sd) or sd <= 0:
        return pd.Series(np.nan, index=x.index, dtype=float)
    return (x - x.mean()) / sd


def prep(d: pd.DataFrame, metric: str) -> pd.DataFrame:
    x = d.sort_values('canonical_name', kind='stable').reset_index(drop=True).copy()
    x['metric_z'] = zscore(x[metric])
    x['effort_z'] = zscore(np.log1p(pd.to_numeric(x['n_climate_cells'], errors='coerce')))
    if 'n_v22_exact_name_records' in x:
        x['literature_attention_z'] = zscore(np.log1p(pd.to_numeric(x['n_v22_exact_name_records'], errors='coerce')))
    return x


def fit_binary(x: pd.DataFrame, outcome: str, predictors: list[str]):
    d = x.dropna(subset=[outcome, 'family', *predictors]).copy()
    if len(d) < 20 or d[outcome].nunique() < 2 or d.family.nunique() < 2:
        return None, d
    X = sm.add_constant(d[predictors], has_constant='add')
    fit = sm.GLM(pd.to_numeric(d[outcome]), X, family=sm.families.Binomial()).fit(
        cov_type='cluster', cov_kwds={'groups': d.family}
    )
    return fit, d


def binary_rows(data: pd.DataFrame, analysis: str, literature_adjusted: bool) -> list[dict]:
    rows = []
    for outcome in OUTCOMES:
        for metric in METRICS:
            x = prep(data, metric)
            predictors = ['metric_z', 'effort_z'] + (['literature_attention_z'] if literature_adjusted else [])
            fit, d = fit_binary(x, outcome, predictors)
            if fit is None:
                rows.append({'analysis': analysis, 'outcome': outcome, 'metric': metric, 'status': 'not_estimable', 'n_species': len(d)})
                continue
            b = float(fit.params['metric_z'])
            se = float(fit.bse['metric_z'])
            lo, hi = b - 1.96 * se, b + 1.96 * se
            row = {
                'analysis': analysis,
                'outcome': outcome,
                'metric': metric,
                'status': 'complete',
                'n_species': int(len(d)),
                'n_families': int(d.family.nunique()),
                'n_positive': int(pd.to_numeric(d[outcome]).sum()),
                'n_negative': int(len(d) - pd.to_numeric(d[outcome]).sum()),
                'estimate': b,
                'odds_ratio': float(np.exp(b)),
                'odds_ratio_ci_low': float(np.exp(lo)),
                'odds_ratio_ci_high': float(np.exp(hi)),
                'wald_p_clustered': float(fit.pvalues['metric_z']),
                'formula': f"{outcome} ~ metric_z + effort_z" + (' + literature_attention_z' if literature_adjusted else ''),
            }
            if literature_adjusted:
                row['literature_attention_estimate'] = float(fit.params['literature_attention_z'])
                row['literature_attention_p'] = float(fit.pvalues['literature_attention_z'])
            rows.append(row)
    return rows


def multinomial_rows(data: pd.DataFrame) -> list[dict]:
    rows = []
    cats = pd.Categorical(data['organization_state'], categories=STATE_ORDER, ordered=False)
    if (cats.codes < 0).any():
        raise ValueError('Unexpected organization_state')
    for metric in METRICS:
        d = prep(data, metric)
        d['state_code'] = pd.Categorical(d['organization_state'], categories=STATE_ORDER).codes
        d = d.dropna(subset=['metric_z', 'effort_z', 'family']).copy()
        X = sm.add_constant(d[['metric_z', 'effort_z']], has_constant='add')
        fit = sm.MNLogit(d.state_code, X).fit(
            method='newton', maxiter=200, disp=False,
            cov_type='cluster', cov_kwds={'groups': d.family}
        )
        # Statsmodels columns correspond non-baseline categories 1 and 2 vs category 0.
        comparisons = {
            0: ('spatial_segregation_only', 'local_coexistence_only'),
            1: ('coexistence_and_segregation', 'local_coexistence_only'),
        }
        for column, (numerator, reference) in comparisons.items():
            b = float(fit.params.loc['metric_z', column])
            se = float(fit.bse.loc['metric_z', column])
            lo, hi = b - 1.96 * se, b + 1.96 * se
            rows.append({
                'analysis': 'three_state_multinomial',
                'metric': metric,
                'numerator_state': numerator,
                'reference_state': reference,
                'n_species': int(len(d)),
                'n_families': int(d.family.nunique()),
                'estimate': b,
                'relative_odds_ratio': float(np.exp(b)),
                'relative_odds_ratio_ci_low': float(np.exp(lo)),
                'relative_odds_ratio_ci_high': float(np.exp(hi)),
                'wald_p_clustered': float(fit.pvalues.loc['metric_z', column]),
                'converged': bool(fit.mle_retvals.get('converged', False)),
                'formula': 'state ~ metric_z + effort_z; non-ordinal multinomial',
            })
    return rows


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument('--dataset', required=True)
    p.add_argument('--attention', required=True)
    p.add_argument('--outdir', required=True)
    args = p.parse_args()

    data = pd.read_csv(args.dataset)
    attention = pd.read_csv(args.attention)
    data = data.merge(attention[['canonical_name', 'n_v22_exact_name_records']], on='canonical_name', how='left', validate='one_to_one')
    if data.n_v22_exact_name_records.isna().any():
        raise SystemExit('Missing literature-attention count for at least one species')

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    adjusted = pd.DataFrame(binary_rows(data, 'literature_attention_adjusted', literature_adjusted=True))
    adjusted.to_csv(outdir / 'cs_literature_attention_adjusted.csv', index=False)

    one_source = data.loc[pd.to_numeric(data['n_resolved_sources'], errors='coerce') == 1].copy()
    one_source_results = pd.DataFrame(binary_rows(one_source, 'one_resolved_positive_source_stratum', literature_adjusted=False))
    one_source_results.to_csv(outdir / 'cs_one_source_stratum.csv', index=False)

    multinomial = pd.DataFrame(multinomial_rows(data))
    multinomial.to_csv(outdir / 'cs_three_state_multinomial.csv', index=False)

    manifest = {
        'n_species': int(len(data)),
        'n_families': int(data.family.nunique()),
        'state_counts': data.organization_state.value_counts().to_dict(),
        'one_source_stratum_n': int(len(one_source)),
        'one_source_C_positive': int(one_source.C_local_coexistence_documented.sum()),
        'one_source_S_positive': int(one_source.S_spatial_segregation_documented.sum()),
        'literature_attention_definition': 'exact accepted canonical binomial mentions in all 12,064 v2.2 title/abstract records; independent of C/S labels but conservative for synonyms',
        'n_resolved_sources_guard': 'n_resolved_sources is outcome-path-derived and is NOT treated as an independent documentation-effort covariate; it is used only to define an equal-positive-source-count stratum.',
        'multinomial_reference': 'local_coexistence_only',
        'semantic_guard': 'All analyses concern documented evidence states, not demonstrated biological absence.',
    }
    (outdir / 'cs_robustness_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(manifest, indent=2))

    if len(adjusted) != 10 or len(one_source_results) != 10 or len(multinomial) != 10:
        raise SystemExit('Unexpected robustness output row count')
    if not (adjusted.status == 'complete').all():
        raise SystemExit('Literature-attention adjusted model not estimable')
    if not (one_source_results.status == 'complete').all():
        raise SystemExit('One-source stratum model not estimable')
    if not multinomial.converged.all():
        raise SystemExit('At least one multinomial fit failed to converge')


if __name__ == '__main__':
    main()
