#!/usr/bin/env python3
"""Source-authenticated, opportunity-preserving comparison of two field flower-colour datasets.

The outcome is within-site quadrat allocation (not neutral Fst, selection, or environmental R2).
The source inventory and quadrat sizes are held fixed throughout every exact-label null.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import zipfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

SIS_MD5 = '008ec196aa5acb769dbfd36e835e32ba'
BOE_INDIV_MD5 = '85990c8a30a8b9ef3a2c045bb6ffb72a'
BOE_PLOT_MD5 = '1c284e1543450aaacf161cb7400d0a5a'
SOURCE_DOIS = {'sisyrinchium': '10.5061/dryad.j8q34', 'boechera': '10.5061/dryad.q0032'}
R = 9999
SEEDS = {'sisyrinchium':2026101307, 'boechera':2026101308}


def original_md5(path: Path, expected: str) -> str:
    d = hashlib.md5(Path(path).read_bytes()).hexdigest()
    if d != expected:
        raise ValueError(f'ORIGINAL_SOURCE_MD5_MISMATCH: {path.name} {d} != {expected}')
    return d


def load_sisyrinchium(zip_path: Path) -> pd.DataFrame:
    original_md5(zip_path, SIS_MD5)
    with zipfile.ZipFile(zip_path) as z:
        bad = z.testzip()
        if bad is not None:
            raise ValueError('CORRUPT_SOURCE_ZIP '+str(bad))
        if z.namelist().count('data/frequency.txt') != 1:
            raise ValueError('MISSING_OR_DUPLICATE_FREQUENCY_SOURCE')
        source = z.read('data/frequency.txt')
    d = pd.read_csv(io.BytesIO(source), sep='\t', encoding='utf-8-sig', engine='python')
    keys=['year','population','line','quad']
    required = keys + ['distance','humidity','v','w','l','s']
    if list(d.columns) != required:
        raise ValueError('SOURCE_SCHEMA_DRIFT '+str(d.columns.tolist()))
    if len(d) != 512 or d.duplicated(keys).any() or d[keys].isna().any().any():
        raise ValueError('SOURCE_DUPLICATE_OR_NAN_KEYS')
    if not set(d.year).issubset({'Y11','Y13'}):
        raise ValueError('SOURCE_YEAR_SCHEMA_CHANGED')
    for c in ['v','w','l','s']:
        y=pd.to_numeric(d[c],errors='raise')
        if not y.notna().all() or not np.all(np.equal(np.floor(y),y)) or not (y>=0).all():
            raise ValueError('INVALID_OR_NEGATIVE_FIELD_COUNT '+c)
        d[c]=y.astype(int)
    d['n']=d['v']+d['w']  # target M-species only: never add l or s
    return d


def load_boechera(individual_path: Path, plot_path: Path) -> tuple[pd.DataFrame,pd.DataFrame,dict]:
    original_md5(individual_path,BOE_INDIV_MD5)
    original_md5(plot_path,BOE_PLOT_MD5)
    lines=Path(individual_path).read_text(encoding='utf-8-sig').splitlines()
    head=lines[0].split()
    required=['Site_name','plot','plantID','Flower_color','Leaf_number','Leaf_area_removed','Elevation','Number_of_fruits']
    if head!=required:
        raise ValueError('BOECHERA_INDIVIDUAL_SCHEMA_CHANGED')
    rows=[]; joined=[]
    for lineno,s in enumerate(lines[1:],2):
        parts=s.split()
        if len(parts)==9:
            if not(parts[2].endswith('.') and parts[3][0].isdigit()):
                raise ValueError('UNEXPECTED_MALFORMED_PLANTID '+str(lineno))
            parts=parts[:2]+[parts[2]+' '+parts[3]]+parts[4:]
            joined.append(lineno)
        if len(parts)!=8:
            raise ValueError('UNEXPECTED_FIELD_COUNT '+str(lineno))
        rows.append(parts)
    d=pd.DataFrame(rows,columns=required)
    if len(d)!=547 or len(joined)!=3 or d.duplicated(['Site_name','plot','plantID']).any():
        raise ValueError('INDIVIDUAL_KEY_OR_SOURCE_COUNT_FAIL')
    if not d.Flower_color.isin(['0','1']).all():
        raise ValueError('UNKNOWN_BOECHERA_COLOUR_CODE')
    d['Flower_color']=d.Flower_color.astype(int)
    plots=pd.read_csv(plot_path,sep=r'\s+',engine='python')
    keys=['Site_name','plot']
    if len(plots)!=26 or plots.duplicated(keys).any():
        raise ValueError('INVALID_BOECHERA_PLOT_KEYS')
    for c in ['Boechera_flowered','Boechera_pigmented']:
        q=pd.to_numeric(plots[c],errors='raise')
        if not q.notna().all() or not (q>=0).all() or not np.equal(q,np.floor(q)).all():
            raise ValueError('BAD_BOECHERA_PLOT_COUNTS')
        plots[c]=q.astype(int)
    if (plots.Boechera_pigmented>plots.Boechera_flowered).any():
        raise ValueError('BOECHERA_COLOUR_EXCEEDS_FLOWERED_COUNT')
    recon=d.groupby(keys,sort=False).Flower_color.agg(size='size',pigmented='sum').reset_index()
    comp=plots[keys+['Boechera_flowered','Boechera_pigmented']].merge(recon,on=keys,how='outer',validate='one_to_one',indicator=True)
    ok=(len(comp)==26 and comp['_merge'].eq('both').all() and (comp.Boechera_flowered==comp['size']).all()
        and (comp.Boechera_pigmented==comp['pigmented']).all())
    if not ok:
        raise ValueError('PLOT_INDIVIDUAL_COUNTS_DO_NOT_RECONCILE')
    plots=plots.copy()
    plots['n']=plots.Boechera_flowered
    plots['v']=plots.Boechera_pigmented
    qa={'n_individuals':len(d),'n_plots':len(plots),'n_sites':plots.Site_name.nunique(),
        'joined_plantid_typos':len(joined),'joined_typo_line_numbers':joined,
        'all_26_plots_match_original_individual_counts':bool(ok),
        'n_missing_foliar_herbivory':int(d.Leaf_area_removed.isin(['NA','']).sum()),
        'n_missing_fruits':int(d.Number_of_fruits.isin(['NA','']).sum())}
    return d,plots,qa


def hold_reason(N:int, m:int, V:int, W:int) -> str|None:
    if N<20:return 'HOLD_LOW_N'
    if m<2:return 'HOLD_LT2_OCCUPIED_QUADRATS'
    if V==0 or W==0:return 'HOLD_MONOMORPHIC'
    if V<5 or W<5:return 'HOLD_LT5_EACH_COLOUR'
    return None


def select_strata(d: pd.DataFrame,species: str,latest_year_only=False) -> list[dict]:
    rows=[]
    if species=='sisyrinchium':
        groups=d.groupby(['year','population'],sort=True)
    elif species=='boechera':
        groups=d.groupby(['Site_name'],sort=True)
    else:raise ValueError('UNKNOWN_SPECIES')
    for key,g in groups:
        if species=='sisyrinchium':
            year,pop=key
            site=pop
            ids=g['line'].astype(str)+':'+g['quad'].astype(str)
        else:
            year='single_census_unreported_year'
            site=key[0] if isinstance(key,tuple) else key
            pop=site
            ids=g['plot'].astype(str)
        n=g['n'].to_numpy(dtype=int)
        v=g['v'].to_numpy(dtype=int)
        record={'year':year,'site':site,'population':pop,'siteyear':year+'|'+site,
                'n':n,'v':v,'plot_ids':ids.tolist(),
                'n_source_quadrats':int(len(g)),
                'n_zero_target_quadrats':int(np.count_nonzero(n==0)),
                'other_species_total':int(g[['l','s']].sum().sum()) if species=='sisyrinchium' else None}
        rows.append(record)
    if latest_year_only and species=='sisyrinchium':
        maxyear={}
        for i,a in enumerate(rows):
            if a['site'] not in maxyear or a['year']>maxyear[a['site']]['year']:
                maxyear[a['site']]=a
        rows=sorted(maxyear.values(),key=lambda a:a['siteyear'])
    return rows


def within_site(n:np.ndarray,v:np.ndarray) -> dict[str,float|int]:
    n=np.asarray(n,dtype=int)
    v=np.asarray(v,dtype=int)
    if n.shape!=v.shape or (n<0).any() or (v<0).any() or (v>n).any():
        raise ValueError('IMPOSSIBLE_QUADRAT_COUNTS')
    zero=int(np.sum(n==0))
    yes=n>0
    nn=n[yes]
    vv=v[yes]
    N=int(nn.sum());V=int(vv.sum());W=N-V;m=len(nn)
    if N<=1 or V==0 or W==0 or m<2:
        raise ValueError('MONOCHROMATIC_OR_NO_WITHIN_SITE_VARIANCE')
    p=V/N
    numerator=float(np.sum((vv-nn*p)**2/nn))
    B=numerator/(N*p*(1-p))
    null_exp=(m-1)/(N-1)
    return {'N':N,'V':V,'W':W,'p_pigmented':p,
            'n_occupied':m,'n_zero_target':zero,'B':B,'expected_B':null_exp,
            'delta':B-null_exp,'observed_mixed':int(np.count_nonzero((vv>0)&(vv<nn))),
            'mixed_fraction':float(np.mean((vv>0)&(vv<nn))),
            'expected_mixed':float(np.mean(expected_mixed(nn,V)))}


def expected_mixed(n:np.ndarray,V:int)->np.ndarray:
    n=np.asarray(n,dtype=int)
    N=int(n.sum())
    if (n<0).any() or V<0 or V>N:raise ValueError('IMPOSSIBLE_COUNTS')
    if N==0:raise ValueError('EMPTY_INVENTORY')
    # Without-replacement, exactly preserving each quadrat's plant opportunity.
    return np.clip(1.-hypergeom.pmf(0,N,V,n)-hypergeom.pmf(n,N,V,n),0.,1.)


def permute_conditional(n:np.ndarray,v:np.ndarray,R:int,seed=None,rng=None) -> dict[str,np.ndarray]:
    if rng is None:rng=np.random.default_rng(seed)
    n=np.asarray(n,dtype=int)
    v=np.asarray(v,dtype=int)
    occupied=n>0
    nn=n[occupied];vv=v[occupied]
    N=int(nn.sum());V=int(vv.sum());m=len(nn)
    if N<=1 or V==0 or V==N or m<=1:raise ValueError('INELIGIBLE_PERMUTATION')
    x=rng.multivariate_hypergeometric(nn,V,size=R,method='marginals')
    if not np.all(x.sum(axis=1)==V) or not np.all((x<=nn)&(x>=0)):
        raise AssertionError('PERMUTATION_INVENTORY_NOT_CONSERVED')
    p=V/N
    B=(((x-nn*p)**2)/nn).sum(axis=1)/(N*p*(1-p))
    mixed=((x>0)&(x<nn)).sum(axis=1)/m
    return {'B':B,'mixed_fraction':mixed}


def summarize_species(records:list[dict])->dict:
    if not records:return {'n_sites':0,'mean_delta':None,'status':'HOLD_NO_ELIGIBLE_SITE'}
    return {'n_sites':len(records),'mean_delta':float(np.mean([x['delta'] for x in records])),
        'mean_B':float(np.mean([x['B'] for x in records])) if all('B' in x for x in records) else None}


def analyze_species(d:pd.DataFrame,species:str,R:int=R,latest_year_only=False,min_occupied=2)->tuple[dict,list[dict],pd.DataFrame,dict]:
    strata=select_strata(d,species,latest_year_only=latest_year_only)
    rng=np.random.default_rng(SEEDS[species] + (10 if latest_year_only else 0) + (1 if min_occupied==3 else 0))
    metric_rows=[];good=[];simulation_B=[];simulation_mix=[]
    full_ledger=[]
    for q in strata:
        n,v=q['n'],q['v'];V=int(v.sum());N=int(n.sum());W=N-V;m=int(np.count_nonzero(n>0))
        status=hold_reason(N,m,V,W)
        if status is None and m<min_occupied:status='HOLD_SENSITIVITY_LT3_QUADRATS'
        row={'siteyear':q['siteyear'],'site':q['site'],'year':q['year'],
             'n_quadrats_total':len(n),'n_occupied':m,'n_zero_target_quadrats':int(np.sum(n==0)),
             'N_target':N,'V_pigmented':V,'W_white':W,'status':status or 'ELIGIBLE',
             'other_species_total':q['other_species_total']}
        for idd,nn,vv in zip(q['plot_ids'],n,v):
            full_ledger.append({'species':species,'site':q['site'],'year':q['year'],
                'quadrat_id':str(idd),'n_target':int(nn),'v_pigmented':int(vv),'w_white':int(nn-vv),
                'status':'NO_TARGET_PLANTS' if nn==0 else 'TARGET_PLANTS_PRESENT'})
        if status is None:
            obs=within_site(n,v)
            row.update(obs)
            nn=n[n>0];vv=v[n>0]
            sims=permute_conditional(nn,vv,R,rng=rng)
            simulation_B.append(sims['B']-obs['expected_B'])
            simulation_mix.append(sims['mixed_fraction'])
            good.append(row)
            row['site_null_B_95pct']=list(map(float,np.quantile(sims['B'],[0.025,0.975])))
            row['site_p_upper']=(1+np.count_nonzero(sims['B']>=obs['B']-1e-12))/(R+1)
            row['site_p_less_mixed']=(1+np.count_nonzero(sims['mixed_fraction']<=obs['mixed_fraction']+1e-12))/(R+1)
        metric_rows.append(row)
    summary={'species':species,'support':{'n_source_siteyear':len(strata),'n_eligible':len(good),
            'n_held':len(strata)-len(good),'n_quadrats_source':sum(x['n_quadrats_total'] for x in metric_rows),
            'n_total_target_plants_source':sum(x['N_target'] for x in metric_rows),
            'n_target_plants_eligible':sum(x['N_target'] for x in good)},
            'inference':'CONDITIONAL_RANDOM_LABEL_BENCHMARK_NOT_CAUSAL','latest_year_only':latest_year_only,'min_occupied':min_occupied}
    null=pd.DataFrame({'replicate':np.arange(1,R+1)})
    if good:
        observed=np.mean([r['delta'] for r in good])
        null['delta']=np.mean(np.vstack(simulation_B),axis=0)
        mix_obs=float(np.mean([r['mixed_fraction'] for r in good]))
        null['mixed_fraction']=np.mean(np.vstack(simulation_mix),axis=0)
        summ={'mean_B_observed':float(np.mean([r['B'] for r in good])),
              'mean_B_expected_finite_random':float(np.mean([r['expected_B'] for r in good])),
              'excess_between_quadrat_diversity_fraction':float(observed),
              'excess_percentage_points':float(100*observed),
              'null_mean_excess':float(null.delta.mean()),
              'null_excess_CI95':list(map(float,np.quantile(null.delta,[0.025,0.975]))),
              'permutation_upper_tail_p':float((1+np.count_nonzero(null.delta.to_numpy()>=observed-1e-12))/(R+1)),
              'observed_site_equal_mixed_quadrat_fraction':mix_obs,
              'expected_site_equal_mixed_fraction':float(np.mean([r['expected_mixed'] for r in good])),
              'null_mean_site_equal_mixed_fraction':float(null.mixed_fraction.mean()),
              'permutation_lower_tail_mixing_p':float((1+np.count_nonzero(null.mixed_fraction.to_numpy()<=mix_obs+1e-12))/(R+1))}
        summary.update(status='COMPUTED_DESCRIPTIVE_AND_CONDITIONAL_NULL',results=summ)
    else:
        summary.update(status='HOLD_NO_ELIGIBLE_POPULATION',results=None)
    return summary,metric_rows,pd.DataFrame(full_ledger),null


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--source-dir',type=Path,default=Path(__file__).resolve().parents[1]/'source')
    parser.add_argument('--outdir',type=Path,default=Path(__file__).resolve().parents[1]/'results')
    parser.add_argument('--reps',type=int,default=R)
    args=parser.parse_args()
    if args.reps<10:raise SystemExit('need at least 10 Monte Carlo draws')
    args.outdir.mkdir(parents=True,exist_ok=True)
    sis=load_sisyrinchium(args.source_dir/'Sisyrinchium_data.zip')
    individual,boe,qa=load_boechera(args.source_dir/'Boechera_natural_populations_individual.txt',
                args.source_dir/'Boechera_natural_populations_plot.txt')
    studies={'sisyrinchium':sis,'boechera':boe}
    main_summary={};ledger=[];poprows=[];null_all=[];sensitivities={}
    for name,d in studies.items():
        summary,rows,quads,null=analyze_species(d,name,R=args.reps)
        main_summary[name]=summary
        ledger.append(quads)
        for rr in rows:poprows.append({'species':name,**rr})
        null_all.append(null.assign(species=name))
        cases=[('at_least_3_occupied',False,3)]
        if name=='sisyrinchium':cases.append(('latest_year_per_named_site',True,2))
        for tag,latest,minocc in cases:
            part,_,_,_=analyze_species(d,name,R=args.reps,latest_year_only=latest,min_occupied=minocc)
            sensitivities[name+'__'+tag]=part
    qa_studies={'boechera':qa,'sisyrinchium':{'n_original_quadrats':len(sis),
        'n_original_siteyear_strata':sis.groupby(['year','population']).ngroups,
        'count_v':int(sis.v.sum()),'count_w':int(sis.w.sum()),
        'count_l_other_species':int(sis.l.sum()),'count_s_other_species':int(sis.s.sum()),
        'n_quadrats_without_M':int((sis.n==0).sum()),'n_missing_humidity':int(sis.humidity.isna().sum())}}
    result={'schema':'fcp_two_independent_field_color_spatial_partition_v1',
        'date_jst':'2026-10-13','protocol':'docs/FCP_TWO_SPECIES_FIELD_MORPH_SPATIAL_PARTITION_PROTOCOL_20261013.md',
        'source_dois':SOURCE_DOIS,'source_sha256':{
           'Sisyrinchium_data.zip':hashlib.sha256((args.source_dir/'Sisyrinchium_data.zip').read_bytes()).hexdigest(),
           'Boechera_natural_populations_individual.txt':hashlib.sha256((args.source_dir/'Boechera_natural_populations_individual.txt').read_bytes()).hexdigest(),
           'Boechera_natural_populations_plot.txt':hashlib.sha256((args.source_dir/'Boechera_natural_populations_plot.txt').read_bytes()).hexdigest()},
        'source_quality':qa_studies,'monte_carlo_reps':args.reps,
        'siteyear_eligibility':{'N_min':20,'V_min':5,'W_min':5,'occupied_quadrats_min':2},
        'source_exposure':'Retrospective exploratory. Site-morph prevalence and raw source field schemas were inspected before registration; color-clumping p values were not.',
        'primary':main_summary,'sensitivities':sensitivities,
        'nonclaims':['The field studies are independent datasets, but populations and quadrats are nested; source sites not random independent global evolutionary replicates.',
          'Not genetic Fst, fitness or selection; null only random placement of recorded colour types across fixed quadrat opportunity within each site-year.',
          'Cannot directly compare the cross-species index as a universal spatial effect due to different quadrat scales, site selection and source ecology.',
          'Original published studies already discuss pigment/abiotic stress and interaction gradients; this reanalysis does not newly discover them.',
          'No original field population estimates from selected Digitalis flower-trait subsamples.']}
    (args.outdir/'result.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    pd.DataFrame(poprows).to_csv(args.outdir/'siteyear_metrics.csv',index=False)
    pd.concat(ledger,ignore_index=True).to_csv(args.outdir/'quadrat_count_ledger.csv',index=False)
    pd.concat(null_all,ignore_index=True).to_csv(args.outdir/'conditional_null_replicates.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    print(json.dumps({name:{'support':x['support'],'results':x['results']} for name,x in main_summary.items()},indent=2,ensure_ascii=False),flush=True)

if __name__=='__main__':main()
