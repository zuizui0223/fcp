#!/usr/bin/env python3
"""Exact summary-identifiable Shannon partition; no reconstructed stand records.

Equal-stand mean proportions and equal-stand mean entropy MUST refer to the same
cohort/year. Rounding enclosures are not sampling confidence intervals.
"""
import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path


def entropy(p):
    p = tuple(float(x) for x in p)
    if not p or any(not math.isfinite(x) or x < 0 for x in p):
        raise ValueError('invalid probability')
    if not math.isclose(sum(p), 1., abs_tol=1e-10, rel_tol=0):
        raise ValueError('probabilities must sum to one')
    return -sum(x * math.log(x) for x in p if x > 0)


def decompose(p, alpha):
    gamma = entropy(p)
    alpha = float(alpha)
    if not math.isfinite(alpha) or alpha < 0 or alpha > gamma + 1e-10:
        raise ValueError('mean within entropy incompatible with pooled entropy')
    beta = gamma - alpha
    return {'H_total': gamma, 'H_within': alpha, 'H_between': beta,
            'between_fraction': beta / gamma if gamma > 1e-14 else None,
            'effective_beta': math.exp(beta)}


def entropy_bounds(center, radius):
    """Concave entropy over a box intersected with the simplex, any dimension."""
    c = [float(x) for x in center]
    if not c or radius < 0 or not math.isfinite(radius) or any(not math.isfinite(x) for x in c):
        raise ValueError('invalid rounding box')
    lo = [max(0., x-radius) for x in c]
    hi = [min(1., x+radius) for x in c]
    if any(a > b for a,b in zip(lo,hi)) or sum(lo) > 1+1e-12 or sum(hi) < 1-1e-12:
        raise ValueError('infeasible rounding box')
    vertices = []
    for free in range(len(c)):
        fixed = [i for i in range(len(c)) if i != free]
        for bits in itertools.product((0,1), repeat=len(fixed)):
            p = [0.] * len(c)
            for i,b in zip(fixed,bits): p[i] = (lo,hi)[b][i]
            p[free] = 1-sum(p)
            if lo[free]-1e-12 <= p[free] <= hi[free]+1e-12:
                p[free] = max(0., min(1., p[free]))
                vertices.append(entropy(p))
    if not vertices: raise ValueError('no feasible rounding vertices')
    # KKT solution: clip a shared water level to each probability interval.
    left, right = 0., 1.
    for _ in range(100):
        level = (left+right)/2
        if sum(min(b,max(a,level)) for a,b in zip(lo,hi)) > 1: right = level
        else: left = level
    pmax = [min(b,max(a,(left+right)/2)) for a,b in zip(lo,hi)]
    return min(vertices), entropy(pmax)


def summary_result(row):
    if int(row['n_stands']) != int(row['alpha_n_stands']):
        raise ValueError('mean frequency and mean entropy cohort mismatch')
    p = [float(row[k]) for k in ('white','pink','purple')]
    alpha = float(row['mean_shannon'])
    out = decompose(p,alpha)
    gmin,gmax = entropy_bounds(p,.0005)
    amin,amax = max(0.,alpha-.005),alpha+.005
    if amax > gmin: raise ValueError('rounding envelope overlaps incompatible summary values')
    out.update(year=int(row['year']), n_stands=int(row['n_stands']),
        inference='descriptive_summary_identity_no_p_value',
        rounding_bounds={'H_total':[gmin,gmax], 'H_within':[amin,amax],
            'H_between':[gmin-amax,gmax-amin],
            'between_fraction':[1-amax/gmin,1-amin/gmax],
            'effective_beta':[math.exp(gmin-amax), math.exp(gmax-amin)]})
    return out


def census_results(path):
    """Same entropy metric on old Linanthus counts, descriptive only."""
    with Path(path).open() as f: rows = list(csv.DictReader(f))
    result = []
    for year in sorted({int(r['year']) for r in rows}):
        r = [x for x in rows if int(x['year']) == year]
        if len(r) != 4 or len({(x['transect'],x['side']) for x in r}) != 4:
            raise ValueError('Linanthus must retain four unique fixed sections')
        counts = [[int(x['blue']),int(x['white'])] for x in r]
        if any(x < 0 for pair in counts for x in pair): raise ValueError('negative count')
        ns = [sum(x) for x in counts]
        if min(ns) < 20:
            result.append({'year':year,'N':sum(ns),'status':'HOLD_NO_PLANTS' if sum(ns)==0 else 'HOLD_SPARSE_SECTIONS'})
            continue
        ps = [[v/n for v in x] for x,n in zip(counts,ns)]
        for method,w in [('equal_section',[.25]*4),('count_weighted',[n/sum(ns) for n in ns])]:
            mean = [sum(a*p[k] for a,p in zip(w,ps)) for k in range(2)]
            out = decompose(mean,sum(a*entropy(p) for a,p in zip(w,ps)))
            out.update(year=year,N=sum(ns),weighting=method,status='OK',
                sections_both_colours=sum(all(v>0 for v in pair) for pair in counts))
            result.append(out)
    return result


def run(summary, counts, outdir):
    with Path(summary).open() as f: rows = list(csv.DictReader(f))
    if len(rows)!=3 or {int(r['year']) for r in rows}!={2021,2022,2023}:
        raise ValueError('exact three published years required')
    results = [summary_result(r) for r in rows]
    res = {'schema':'fcp_hesperis_summary_shannon_partition_v1',
        'status':'SUMMARY_IDENTITY_COMPLETE_RAW_CENSUS_NOT_RECOVERED',
        'source_doi':'10.1093/aob/mcag035',
        'source_locator':'Table 2 all stands mean morph frequencies; Results neutral predictions 1 and 2 mean H-prime',
        'weighting':'equal stand within each year', 'logarithm':'natural',
        'source_p_half_rounding_unit':.0005,'source_alpha_half_rounding_unit':.005,
        'hesperis':results,'linanthus_same_metric':census_results(counts),
        'source_hashes':{Path(p).name:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in (summary,counts)},
        'limits':['Rounding intervals only, not sampling confidence intervals.',
            'No per-stand census records or sample denominators were recovered or manufactured.',
            'Empirical plug-in entropy, not bias-corrected population diversity.',
            'Annual sets differ and overlap; three years are not three independent replicates.',
            'The 105-stand longitudinal H-prime means cannot be mixed with all-stand colour shares.',
            'Linanthus four transect-sections and Hesperis stands have incomparable sampling scales.',
            'No multi-species general law, causal mechanism, or new p-value.',
            'The original Linanthus Gini-Simpson results are unchanged.']}
    outdir = Path(outdir); outdir.mkdir(parents=True,exist_ok=True)
    (outdir/'result.json').write_text(json.dumps(res,indent=2,ensure_ascii=False)+'\n')
    for x in results: print(json.dumps(x))
    return res


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--summary',type=Path,required=True)
    p.add_argument('--counts',type=Path,required=True)
    p.add_argument('--outdir',type=Path,required=True)
    a=p.parse_args(); run(a.summary,a.counts,a.outdir)
