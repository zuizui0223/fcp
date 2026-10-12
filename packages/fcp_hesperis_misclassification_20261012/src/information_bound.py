#!/usr/bin/env python3
"""Bound classification-error effects on published Hesperis stand-colour entropy.

Source-exposed *hypothetical* sensitivity, not a measured misclassification rate,
not original census recovery, not error-corrected biological variation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

PARENT_SHA = '980b50bfe99f69b31b57d273a2995c814e280cf552a52a1eefaf17e828b4faf1'
SOURCE_DOI = '10.1093/aob/mcag035'
SOURCE_SCHEMA = 'fcp_hesperis_summary_shannon_partition_v1'


def fano_entropy(e: float, k: int = 3) -> float:
    """Upper bound on H(recorded colour | true colour), in nats.

    e is an *upper bound* on the weighted fraction misclassified. For e >=
    1-1/k the entropy ceiling log(k) is the relevant monotone envelope.
    """
    if not isinstance(k, int) or k < 2:
        raise ValueError('at least two categories are needed')
    if not math.isfinite(e) or not 0 <= e <= 1:
        raise ValueError('classification error upper bound must be finite in [0,1]')
    if e == 0:
        return 0.0
    if e >= 1-1/k:
        return math.log(k)
    return -e*math.log(e)-(1-e)*math.log1p(-e)+e*math.log(k-1)


def safe_threshold(B_min: float, k: int = 3) -> float:
    """Largest weighted error epsilon with Fano(epsilon)<B_min (at equality)."""
    if not math.isfinite(B_min) or B_min < 0 or B_min > math.log(k)+1e-12:
        raise ValueError('observed Shannon between-component outside entropy range')
    if B_min <= 0:
        return 0.0
    if B_min >= math.log(k):
        return 1-1/k
    lo, hi = 0.0, 1-1/k
    for _ in range(65):
        mid = (lo+hi)*0.5
        if fano_entropy(mid,k) < B_min:
            lo = mid
        else:
            hi = mid
    return (lo+hi)*0.5


def mutual_info(w: list[float], group_morph_probs: list[list[float]]) -> float:
    """I(stand;colour) in nats, keeping externally supplied stand weights."""
    if len(w) != len(group_morph_probs) or not w or abs(sum(w)-1)>1e-9:
        raise ValueError('stand weights must match and sum to one')
    k=len(group_morph_probs[0]); allrows=[]
    for wt,vec in zip(w,group_morph_probs):
        if wt < 0 or len(vec)!=k or not all(math.isfinite(x) and x>=0 for x in vec) or abs(sum(vec)-1)>1e-9:
            raise ValueError('invalid stand category distribution')
        allrows.append([wt*v for v in vec])
    pooled=[sum(row[j] for row in allrows) for j in range(k)]
    def entropy(vec):
        return -sum(v*math.log(v) for v in vec if v>0)
    return max(0., entropy(pooled)-sum(wt*entropy(vec) for wt,vec in zip(w,group_morph_probs)))


def apply_stand_channels(true_probs: list[list[float]], channels: list[list[list[float]]]) -> list[list[float]]:
    """Apply each stand's stochastic true-colour -> recorded-colour channel."""
    if len(true_probs)!=len(channels):
        raise ValueError('number of channels does not match stands')
    out=[]
    for p,C in zip(true_probs,channels):
        k=len(p)
        if len(C)!=k or any(len(r)!=k or not all(0<=t<=1 and math.isfinite(t) for t in r) or abs(sum(r)-1)>1e-9 for r in C):
            raise ValueError('each channel must have stochastic rows')
        if abs(sum(p)-1)>1e-9:
            raise ValueError('invalid true colour composition')
        out.append([sum(p[a]*C[a][b] for a in range(k)) for b in range(k)])
    return out


def verify_parent(path: Path) -> dict:
    b=Path(path).read_bytes()
    sha=hashlib.sha256(b).hexdigest()
    if sha!=PARENT_SHA:
        raise ValueError('parent result provenance SHA256 mismatch')
    d=json.loads(b)
    if d.get('schema')!=SOURCE_SCHEMA:
        raise ValueError('parent result schema mismatch')
    if d.get('source_doi')!=SOURCE_DOI or d.get('weighting')!='equal stand within each year':
        raise ValueError('parent source or weighting mismatched')
    if len(d.get('hesperis',[]))!=3 or {x['year'] for x in d['hesperis']}!={2021,2022,2023}:
        raise ValueError('parent time coverage drift')
    if {x['year']: x['n_stands'] for x in d['hesperis']} != {2021:132,2022:112,2023:116}:
        raise ValueError('parent census stand denominators mismatch')
    for a in d['hesperis']:
        if not 0 < a['rounding_bounds']['H_between'][0] <= a['H_between'] <= a['rounding_bounds']['H_between'][1] < math.log(3):
            raise ValueError('invalid parent between-group Shannon component')
    return {'status':'PASS','parent_sha256':sha,'n_years':3,'source_not_raw_stand_counts':True}


def analyze_parent(path: Path, errors: list[float]) -> dict:
    provenance=verify_parent(path)
    d=json.loads(Path(path).read_text())
    values=sorted(set(float(e) for e in errors))
    if not values or values[0]<0 or values[-1]>1:
        raise ValueError('error scenario grid must be in [0,1]')
    a=[]
    for row in sorted(d['hesperis'], key=lambda x:x['year']):
        bmin=float(row['rounding_bounds']['H_between'][0]); bmax=float(row['rounding_bounds']['H_between'][1])
        hgamma_max=float(row['rounding_bounds']['H_total'][1]); threshold=safe_threshold(bmin,3)
        scenarios=[]
        for e in values:
            phi=fano_entropy(e,3)
            B_true_min=max(0,bmin-phi)
            Hgamma_true_max=min(math.log(3),hgamma_max+phi)
            ratio_min=B_true_min/Hgamma_true_max if Hgamma_true_max>0 else 0
            scenarios.append({'error_rate_upper':e,'fano_information_ceiling_nats':phi,
                              'guaranteed_true_sample_B_min':B_true_min,
                              'guaranteed_true_sample_fraction_min':ratio_min,
                              'assumption':'sampled correct labels differ from recorded colours on at most this equal-stand-weighted fraction'})
        a.append({'year':row['year'],'n_stands':row['n_stands'],
                  'B_observed_nats':row['H_between'],'B_original_rounding_min':bmin,
                  'B_original_rounding_max':bmax,'Hgamma_original_rounding_max':hgamma_max,
                  'threshold_error_rate':threshold,
                  'assumptions':'three colour classes; stand weights equal; unknown possibly stand-varying misclassification; parent publication rounding envelope',
                  'scenarios':scenarios})
    return {
        'schema':'fcp_hesperis_error_information_sensitivity_v1',
        'scope':'DETERMINISTIC_LABEL_ERROR_SENSITIVITY_NOT_FIELD_DATA',
        'source_doi':SOURCE_DOI,'parent':provenance,
        'verified_empirical_error_rate':None,'n_color_categories':3,
        'annual':a,
        'theorem':{'identity':'B_shannon=I(stand;recorded_colour)',
                   'general_bound':'I(stand;true_sample_colour) >= max(0,B_obs - h2(e)-e*ln(k-1))',
                   'constant_channel':'if recorded_colour is independent of stand conditional on true_colour, then I(stand;recorded)<=I(stand;true)',
                   'proof':'I(G;Z)<=I(G;X,Z)=I(G;X)+I(G;Z|X)<=I(G;X)+H(Z|X)<=I(G;X)+h2(e)+e*ln(k-1)',
                   'error_weight':'equal stand weight, then uniform among the sampled records of the stand',
                   'validity':'latent correctly identified morphs of the SAME SAMPLED PLANTS, not true finite-stand census compositions'},
        'nonclaims':['The original study did not measure a true misclassification rate; epsilon scenarios are not observed.',
                     'These lower bounds do not remove finite sampling noise or ascertainment bias.',
                     'The threshold is a necessary lower bound on misclassification if all properly classified sample compositions were equal; it is not a confidence interval or p-value.',
                     'The data processing result requires stand-invariant classification channel; stand-specific errors need a rate bound.',
                     'Do not combine this deterministic bound with PR174 sampling-only p-upper as if it were a calibrated joint model.',
                     'No causal selection, drift, population genetics, cross-species generality, or longitudinal effect is established.']
    }


def main():
    pa=argparse.ArgumentParser()
    pa.add_argument('--parent',type=Path,required=True)
    pa.add_argument('--outdir',type=Path,required=True)
    args=pa.parse_args()
    result=analyze_parent(args.parent,[0,.005,.01,.02,.03,.04,.05,.075,.10])
    args.outdir.mkdir(parents=True,exist_ok=True)
    (args.outdir/'result.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    import csv
    with (args.outdir/'annual_sensitivity.csv').open('w',newline='') as f:
        writer=csv.writer(f)
        writer.writerow(['year','n_stands','B_observed','B_rounding_min','required_misclassification_to_fully_explain','scenario_error_upper','fano_ceiling','correct_sample_B_lower','correct_sample_between_fraction_lower'])
        for a in result['annual']:
            for s in a['scenarios']:
                writer.writerow([a['year'],a['n_stands'],a['B_observed_nats'],a['B_original_rounding_min'],a['threshold_error_rate'],s['error_rate_upper'],s['fano_information_ceiling_nats'],s['guaranteed_true_sample_B_min'],s['guaranteed_true_sample_fraction_min']])
    print(json.dumps({'status':'PASS','three_year_error_thresholds':{x['year']:x['threshold_error_rate'] for x in result['annual']},'source_sha256':result['parent']['parent_sha256'],'new_p_values':False},ensure_ascii=False))

if __name__=='__main__':main()
