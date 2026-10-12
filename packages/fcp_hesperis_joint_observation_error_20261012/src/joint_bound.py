#!/usr/bin/env python3
"""Conditional joint classification + sampling observation-null upper bound.

Source-exposed POSTHOC sensitivity; only published summary and source-locked
parent results. No raw field counts or measured label-error rates were used.
Run: python src/joint_bound.py --outdir results
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
PARENT_SHA = '980b50bfe99f69b31b57d273a2995c814e280cf552a52a1eefaf17e828b4faf1'
ERROR_RATES = (0.0, 0.0025, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05)
SAMPLE_MINIMA = (25, 50, 75, 100)
ALPHA_FAMILY = 0.05
FAMILY_COMPARISONS = 3
K = 3


def verify_parent(path: Path) -> dict:
    path = Path(path)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != PARENT_SHA:
        raise ValueError(f'PR173 parent summary hash mismatch: {digest}')
    source = json.loads(path.read_text(encoding='utf-8'))
    if (source['schema'] != 'fcp_hesperis_summary_shannon_partition_v1'
            or source['weighting'] != 'equal stand within each year'
            or source['logarithm'] != 'natural'):
        raise ValueError('Unanticipated population, weights or entropy base')
    annual = source['hesperis']
    if [x['year'] for x in annual] != [2021, 2022, 2023]:
        raise ValueError('Expected exactly the three published annual summaries')
    if [x['n_stands'] for x in annual] != [132, 112, 116]:
        raise ValueError('Unexpected stand opportunities')
    for row in annual:
        lo, hi = row['rounding_bounds']['H_between']
        if not (0 < lo <= row['H_between'] <= hi <= math.log(K)):
            raise ValueError('Bad parent Shannon bounds')
    return {'sha256':digest,'annual':annual,'n_annual':len(annual)}


def binary_entropy(p: float) -> float:
    if not 0 <= p <= 1:
        raise ValueError('binary event mass must be in [0,1]')
    if p in (0, 1):
        return 0.0
    return -p * math.log(p) - (1-p) * math.log1p(-p)


def fano_ceiling(epsilon: float, k: int = K) -> float:
    """Monotonic worst-case H(Z|X) for an upper error-budget on sampled labels."""
    if not 0 <= epsilon <= 1 or k < 2:
        raise ValueError('Error-budget in [0,1] and at least two colour classes')
    if epsilon >= (k - 1) / k:
        return math.log(k)
    return binary_entropy(epsilon) + epsilon * math.log(k-1)


def sample_null_upper(gap: float, min_sample: int, k: int = K) -> float:
    """No-between-stand-independence PR174 Chernoff upper U(gap;n0,k)."""
    if not math.isfinite(gap) or min_sample < 0 or k < 2:
        raise ValueError('Invalid null input')
    if min_sample == 0 or gap <= 0 or min_sample * gap <= k-1:
        return 1.0
    a = k-1
    ln_u = -min_sample * gap + a * (1+math.log(min_sample*gap/a))
    return min(1.0, math.exp(ln_u))


def joint_year_upper(gap_lower: float, min_sample: int, epsilon: float, eta: float=0.0,
                     k: int=K) -> float:
    """The null event inclusion: P(B_recorded>=b) <= eta + U(b-Phi(eps)).

    An unconditional statement needs Pr(actual equal-stand-weighted
    misclassification fraction > eps) <= eta for every resampling world.
    Pathwise guarantee takes eta=0. Neither is supplied by the source.
    """
    if not 0<=eta<=1:
        raise ValueError('failure probability eta must lie in [0,1]')
    corrected_gap = gap_lower - fano_ceiling(epsilon,k)
    return min(1.0, eta + sample_null_upper(corrected_gap,min_sample,k))


def joint_family_upper(gap_lower: float, min_sample: int, epsilon: float,
                       eta: float=0.0, k: int=K, comparisons: int=FAMILY_COMPARISONS) -> float:
    if comparisons<1:raise ValueError('Invalid multiplicity count')
    return min(1.0, comparisons*joint_year_upper(gap_lower,min_sample,epsilon,eta,k))


def max_error_for_family_alpha(gap_lower:float,min_sample:int,alpha:float=ALPHA_FAMILY,
                              eta:float=0.0,k:int=K,comparisons:int=FAMILY_COMPARISONS):
    """Largest admissible epsilon in [0,(k-1)/k] with conditional upper <= alpha.

    None means even epsilon=0 gives no rejection certificate.
    """
    if not 0<alpha<1:raise ValueError('alpha must be strictly between zero and one')
    if joint_family_upper(gap_lower,min_sample,0.0,eta,k,comparisons)>alpha:
        return None
    lo,hi=0.0,1-1/k
    for _ in range(84):
        mid=(lo+hi)/2
        if joint_family_upper(gap_lower,min_sample,mid,eta,k,comparisons)<=alpha:
            lo=mid
        else:hi=mid
    return lo


def empirical_stand_info(groups: list[list[int]],k: int=K)->float:
    """Equal-stand weighted I(G;colour), with one uniform sampled record/stand."""
    if not groups or any(not group for group in groups):
        raise ValueError('Nonempty sampled labels needed for every stand')
    m=len(groups)
    joint=[0.0]*k
    within=0.0
    for group in groups:
        n=len(group)
        p=[group.count(j)/n for j in range(k)]
        if any(x not in range(k) for x in group):raise ValueError('Unknown morph state')
        for j in range(k):joint[j]+=p[j]/m
        within-=sum(x*math.log(x) for x in p if x>0)/m
    global_h=-sum(x*math.log(x) for x in joint if x>0)
    return max(0.0,global_h-within)


def analyze(parent:Path):
    record=verify_parent(parent)
    out={
        'schema':'fcp_hesperis_joint_observation_error_20261012_v1',
        'inference_status':'CONDITIONAL_JOINT_SAMPLING_AND_CLASSIFICATION_BUDGET',
        'source_doi':'10.1093/aob/mcag035',
        'parent_PR173_sha256':record['sha256'],
        'prior_sampling_only_PR174':'https://github.com/zuizui0223/fcp/pull/174',
        'prior_classification_only_PR176':'https://github.com/zuizui0223/fcp/pull/176',
        'null_hypothesis':'same_true_finite_census_colour_composition_all_stands',
        'sample_size_semantics':'census_or_simple_random_without_replacement_n_ge_n0',
        'classification_error_semantics':'pathwise realized equal-stand-weighted sampled-label mismatch <= epsilon; optional probability eta of violation',
        'k_colours':K,'n_years':record['n_annual'],'familywise_adjustment':'Bonferroni across three annual comparisons, no independence needed',
        'alpha_family':ALPHA_FAMILY,'error_budget_failure_eta':0.0,
        'sample_minimum_scenarios':list(SAMPLE_MINIMA),
        'error_upper_scenarios':list(ERROR_RATES),
        'annual':{},
        'measured_misclassification_rate':None,'raw_stand_counts_recovered':False,
        'new_field_measurements':False,'new_hypothesis_test':'retrospective conditional null upper-bound extension, not preregistered confirmation',
        'limits':[
            'No empirical observation of classification error, and no proof that any epsilon budget is satisfied.',
            'A p-upper under joint assumptions is not the product of sibling PR p-values.',
            'Null concerns identical finite census compositions, not a common source gene pool or population demographic process.',
            'No validation of actual field sampling, colour-dependent sampling, cryptic genets, error structure, or field observers.',
            'Cannot infer corrected entropy effect size, temporal trend, fitness selection, drift exclusion, genetic Fst or adaptation.',
            'The sample-null model may be invalid if field counting is not census or colour-independent random sampling.',
            'The three-year correction does not account for all historical exploratory questions or model checks.',
        ]
    }
    for a in record['annual']:
        year=a['year'];b=a['rounding_bounds']['H_between'][0]
        row={
            'n_stands':a['n_stands'], 'published_between_nats':a['H_between'],
            'published_between_fraction':a['between_fraction'],
            'rounding_conservative_between_lower_nats':b,
            'critical_error_n100':max_error_for_family_alpha(b,100),
            'critical_error_n75':max_error_for_family_alpha(b,75),
            'critical_error_n50':max_error_for_family_alpha(b,50),
            'critical_error_n25':max_error_for_family_alpha(b,25),
            'scenario':{}
        }
        for n in SAMPLE_MINIMA:
            row['scenario'][n]={}
            for eps in ERROR_RATES:
                delta=fano_ceiling(eps)
                adj=joint_family_upper(b,n,eps)
                row['scenario'][n][f'{eps:.4f}']={
                    'epsilon_upper':eps,'fano_bound_nats':delta,
                    'minimum_correctly_classified_sample_gap_nats':max(0,b-delta),
                    'test_gap_after_classification_budget_nats':b-delta,
                    'sampling_only_tail_upper':sample_null_upper(b-delta,n),
                    'unadjusted_joint_upper':joint_year_upper(b,n,eps),
                    'adjusted_bound':adj,
                    'reject_joint_null_under_pathwise_budget_at_0_05':adj<=ALPHA_FAMILY,
                }
        out['annual'][year]=row
    return out


def write_csv(result,path:Path):
    with Path(path).open('w',newline='',encoding='utf-8') as f:
        wr=csv.DictWriter(f,fieldnames=['year','stand_n','min_sample_n0','epsilon_upper','fano_ceiling','B_lower','B_after_error','p_upper_one_year','p_upper_three_year_Bonferroni','rejection_certificate'])
        wr.writeheader()
        for year,r in result['annual'].items():
            for n,row in r['scenario'].items():
                for eps,v in row.items():
                    wr.writerow({'year':year,'stand_n':r['n_stands'],'min_sample_n0':n,
                      'epsilon_upper':eps,'fano_ceiling':v['fano_bound_nats'],
                      'B_lower':r['rounding_conservative_between_lower_nats'],
                      'B_after_error':v['test_gap_after_classification_budget_nats'],
                      'p_upper_one_year':v['unadjusted_joint_upper'],
                      'p_upper_three_year_Bonferroni':v['adjusted_bound'],
                      'rejection_certificate':v['reject_joint_null_under_pathwise_budget_at_0_05']})


def main():
    p=argparse.ArgumentParser();p.add_argument('--parent',type=Path,default=PROJECT/'data'/'parent_entropy_result.json')
    p.add_argument('--outdir',type=Path,default=PROJECT/'results')
    a=p.parse_args();a.outdir.mkdir(exist_ok=True,parents=True)
    result=analyze(a.parent)
    (a.outdir/'result.json').write_text(json.dumps(result,indent=2,ensure_ascii=False,sort_keys=True)+'\n',encoding='utf-8')
    write_csv(result,a.outdir/'joint_sensitivity.csv')
    print(json.dumps({'status':'COMPLETE_CONDITIONAL_SCENARIO',
        'n_years':result['n_years'],'critical_at_n100':{y:r['critical_error_n100'] for y,r in result['annual'].items()},
        'measured_error_rate':None}))

if __name__=='__main__':main()
