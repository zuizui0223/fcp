#!/usr/bin/env python3
"""Out-of-year geographic transfer of field-counted Linanthus flower colour.

Four fixed, broad census sections (north/south x east/west of a ravine) and
1992–1998 observations. Source-exposed retrospective single-study diagnosis;
no plant-/section-/year-level independence assumptions for p values.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from collections import defaultdict

SOURCE_SHA256='fa408c7ad6e521a8410252d903b024a057b1ee4f7e27406927e77fcc22bf4ac7'
SECTIONS=('north_east','north_west','south_east','south_west')
MODELS=('GLOBAL','SIDE','SECTION')


def load_rows(path:Path)->list[dict]:
    path=Path(path)
    h=hashlib.sha256(path.read_bytes()).hexdigest()
    if h!=SOURCE_SHA256:
        raise ValueError('Census source SHA256 mismatch: '+h)
    with path.open(newline='',encoding='utf-8') as f:
        reader=csv.DictReader(f)
        required={'year','transect','side','blue','white','reported_white_pct'}
        if not required.issubset(reader.fieldnames or []):raise ValueError('Missing source columns')
        rows=[]
        for row in reader:
            year=int(row['year']);transect=row['transect'];side=row['side']
            blue=int(row['blue']);white=int(row['white'])
            if blue<0 or white<0 or side not in ('east','west') or transect not in ('north','south'):
                raise ValueError('Invalid census value')
            n=blue+white; printed=row['reported_white_pct']; pct=float(printed) if printed.strip() else None
            p=white/n if n else None
            if (n==0)!=(pct is None): raise ValueError('Missing colour percentage inconsistent with no plants')
            rows.append({'year':year,'transect':transect,'side':side,'blue':blue,'white':white,'n':n,
                'p':p,'reported_white_pct':pct,'section':transect+'_'+side})
    keys={(r['year'],r['section']) for r in rows}
    if len(rows)!=28 or len(keys)!=28 or {r['year'] for r in rows}!=set(range(1992,1999)):
        raise ValueError('Frozen source structure drift')
    if sorted({r['section'] for r in rows})!=sorted(SECTIONS):
        raise ValueError('Not all four fixed sections available')
    if sum(r['n'] for r in rows)!=100500:raise ValueError('Original total census did not agree')
    return sorted(rows,key=lambda x:(x['year'],x['section']))


def eligible_years(rows:list[dict], min_n:int=20)->list[int]:
    if min_n<1:raise ValueError('Positive minimum needed')
    years=sorted({r['year'] for r in rows})
    return [y for y in years if len([r for r in rows if r['year']==y])==4
        and {r['section'] for r in rows if r['year']==y}==set(SECTIONS)
        and min(r['n'] for r in rows if r['year']==y)>=min_n]


def compute_variance(rows:list[dict], years:list[int])->dict:
    sub=sorted([r for r in rows if r['year'] in years], key=lambda x:(x['year'],x['section']))
    if not years or len(sub)!=4*len(years) or any(x['p'] is None for x in sub):
        raise ValueError('Non-rectangular complete four-section census')
    p=sum(r['p'] for r in sub)/len(sub)
    site={s:sum(r['p'] for r in sub if r['section']==s)/len(years) for s in SECTIONS}
    yy={y:sum(r['p'] for r in sub if r['year']==y)/4 for y in years}
    vt=sum((r['p']-p)**2 for r in sub)/len(sub)
    vs=sum((x-p)**2 for x in site.values())/4
    vy=sum((x-p)**2 for x in yy.values())/len(years)
    vi=sum((r['p']-site[r['section']]-yy[r['year']]+p)**2 for r in sub)/len(sub)
    if not math.isclose(vt,vs+vy+vi,rel_tol=1e-10,abs_tol=1e-14):
        raise AssertionError('Space, time, and interaction must reassemble')
    return {'n_years':len(years),'n_section_years':len(sub), 'grand_white_fraction_equal_section_year':p,
        'total_variance':vt,'spatial_variance':vs,'temporal_variance':vy,'interaction_variance':vi,
        'spatial_fraction':vs/vt if vt else None,
        'temporal_fraction':vy/vt if vt else None,
        'interaction_fraction':vi/vt if vt else None,
        'section_white_fraction_means':site,'year_white_fraction_means':{str(y):yy[y] for y in years},
        'mean_naive_binomial_fraction_sampling_variance':sum(r['p']*(1-r['p'])/r['n'] for r in sub)/len(sub)}


def cv_predictions(rows:list[dict], years:list[int])->list[dict]:
    sub=sorted([r for r in rows if r['year'] in years],key=lambda r:(r['year'],r['section']))
    if len(sub)!=4*len(years) or len(years)<2:raise ValueError('Need >=2 complete years')
    output=[]
    for target in years:
        train=[r for r in sub if r['year']!=target]
        if len(train)!=4*(len(years)-1):raise AssertionError('Year-leakage control')
        fullmean=sum(r['p'] for r in train)/len(train)
        side={s:sum(r['p'] for r in train if r['side']==s)/(2*(len(years)-1)) for s in ('east','west')}
        section={s:sum(r['p'] for r in train if r['section']==s)/(len(years)-1) for s in SECTIONS}
        for observation in [r for r in sub if r['year']==target]:
            preds={'GLOBAL':fullmean,'SIDE':side[observation['side']],'SECTION':section[observation['section']]}
            for model in MODELS:
                pred=preds[model]
                output.append({'year':target,'section':observation['section'],'side':observation['side'],
                    'transect':observation['transect'],'n':observation['n'],'white':observation['white'],
                    'blue':observation['blue'],'observed_white_fraction':observation['p'],'model':model,
                    'prediction':pred,'error_sq':(observation['p']-pred)**2,
                    'n_training_section_years':len(train),'n_training_years':len(years)-1})
    if len(output)!=len(years)*4*len(MODELS):raise AssertionError('Prediction coverage fail')
    return output


def group_annual(rows,years):
    result=[]
    for year in years:
        g=[r for r in rows if r['year']==year]
        side_stats={}
        for side in ('west','east'):
            d=[x for x in g if x['side']==side]
            side_stats[side]={
              'counted_plants':sum(x['n'] for x in d),
              'white_count':sum(x['white'] for x in d),
              'pooled_white_fraction':sum(x['white'] for x in d)/sum(x['n'] for x in d),
              'equal_transect_white_fraction':sum(x['p'] for x in d)/2,
            }
        result.append({'year':year,'counted_plants':sum(r['n'] for r in g),
          'west':side_stats['west'],'east':side_stats['east'],
          'east_minus_west_equal_transect':side_stats['east']['equal_transect_white_fraction']-side_stats['west']['equal_transect_white_fraction'],
          'four_sections_each_has_both_colours':all(x['blue']>0 and x['white']>0 for x in g)})
    return result


def report_case(rows,years):
    v=compute_variance(rows,years)
    pred=cv_predictions(rows,years)
    per_year=[]
    for year in years:
        mse={name:sum(x['error_sq'] for x in pred if x['year']==year and x['model']==name)/4 for name in MODELS}
        per_year.append({'year':year,**{'mse_'+name.lower():mse[name] for name in MODELS},
            'side_vs_global_relative_reduction':(mse['GLOBAL']-mse['SIDE'])/mse['GLOBAL'] if mse['GLOBAL']>0 else None,
            'section_increment_vs_side':mse['SIDE']-mse['SECTION']})
    mean_mse={m:sum(y['mse_'+m.lower()] for y in per_year)/len(years) for m in MODELS}
    return {'years':years,'variance_partition':v,'cv_predictions':pred,'heldout_years':per_year,
        'mean_frequency_MSE':mean_mse,
        'relative_error_reduction_from_side':(mean_mse['GLOBAL']-mean_mse['SIDE'])/mean_mse['GLOBAL'] if mean_mse['GLOBAL'] else None,
        'extra_reduction_from_section':mean_mse['SIDE']-mean_mse['SECTION'],
        'annual_side_contrasts':group_annual(rows,years)}


def run(source)->dict:
    source=Path(source)
    rows=load_rows(source)
    main=eligible_years(rows,20)
    sensitivity=eligible_years(rows,100)
    def miss(x):
        printed=x['reported_white_pct'];calc=100*x['p'] if x['p'] is not None else None
        return None if printed is None or abs(printed-calc)<=0.051 else {
            'year':x['year'],'section':x['section'],'reported_white_pct':printed,'actual_white_pct':calc,
            'source_blue':x['blue'],'source_white':x['white'],'source_n':x['n']}
    mismatches=[m for row in rows if (m:=miss(row)) is not None]
    counts={str(y):sum(r['n'] for r in rows if r['year']==y) for y in sorted({x['year'] for x in rows})}
    excluded=[{'year':y,'n_plant_records':counts[str(y)],
       'n_empty_sections':sum(r['n']==0 for r in rows if r['year']==y),
       'status':'HOLD_NO_PLANTS' if counts[str(y)]==0 else 'HOLD_MIN_SECTION_N'}
      for y in range(1992,1999) if y not in main]
    primary=report_case(rows,main);strong=report_case(rows,sensitivity)
    nby=[counts[str(y)] for y in main]
    return {'schema':'fcp_linanthus_out_of_year_field_frequency_spacetime_v1',
      'date_jst':'2026-10-13','role':'post_publication_exploratory_one_site_field_benchmark',
      'source_publication':'Schemske & Bierzychudek 2007, Evolution 61:2528-2543, doi:10.1111/j.1558-5646.2007.00219.x, Appendix',
      'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
      'n_original_rows':len(rows), 'n_counted_plant_records_all_years':sum(r['n'] for r in rows),
      'n_fixed_sections':4,'n_all_years':7,
      'count_by_year':counts,'excluded_years':excluded,
      'printed_percent_mismatches':mismatches,
      'n_years_main':len(main),'n_section_years_main':4*len(main),
      'counted_plants_in_main':sum(nby),'annual_abundance_max_over_min_main':max(nby)/min(nby),
      'primary':{
         'eligible_min_section_n':20,
         'out_of_year_global_mse':primary['mean_frequency_MSE']['GLOBAL'],
         'out_of_year_side_mse':primary['mean_frequency_MSE']['SIDE'],
         'out_of_year_section_mse':primary['mean_frequency_MSE']['SECTION'],
         **primary},
      'sensitivity_min100':{'eligible_min_section_n':100,**strong},
      'claims':['Stable observed west-versus-east frequency placement over evaluable years at one fixed site.',
        'A heldout-year geographic-side prior greatly reduces four-sector fraction prediction error relative to one global prior.'],
      'nonclaims':['NOT a global biological general law, independent species replication or multiple independent landscapes.',
        'NOT phylogenetic or population genetic F_ST, individual genotype frequencies, or the original 2m-quadrat-level diversity.',
        'No causal support for selection vs drift vs gene flow or specific moisture effects in this reanalysis.',
        'No independent p-value because only 5 correlated repeated-census years in one site.',
        'These frequency-variance fractions are not Shannon or Gini-Simpson diversity allocation.',
        'LOYO is a source-exposed retrospective diagnostic, not genuinely prospective forecasting.'],
      'original_main_H1H2_and_pr172_results_unchanged':True}


def cli():
    ap=argparse.ArgumentParser()
    ap.add_argument('--source',type=Path,default=Path(__file__).resolve().parents[1]/'source/linanthus_counts.csv')
    ap.add_argument('--outdir',type=Path,default=Path(__file__).resolve().parents[1]/'results')
    x=ap.parse_args();z=run(x.source);x.outdir.mkdir(parents=True,exist_ok=True)
    (x.outdir/'result.json').write_text(json.dumps(z,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf8')
    for name,rows in [('heldout_section_model_predictions.csv',z['primary']['cv_predictions']),
                      ('yearly_side_summary.csv',z['primary']['annual_side_contrasts'])]:
        if name.startswith('yearly'):
            rows=[{'year':r['year'],'counted_plants':r['counted_plants'],
                   'west_n':r['west']['counted_plants'],'east_n':r['east']['counted_plants'],
                   'west_white_fraction_count_weighted':r['west']['pooled_white_fraction'],
                   'east_white_fraction_count_weighted':r['east']['pooled_white_fraction'],
                   'west_white_fraction_equal_transect':r['west']['equal_transect_white_fraction'],
                   'east_white_fraction_equal_transect':r['east']['equal_transect_white_fraction'],
                   'east_minus_west_equal_transect':r['east_minus_west_equal_transect']} for r in rows]
        with (x.outdir/name).open('w',newline='',encoding='utf8') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print(json.dumps({'main_years':z['primary']['years'],
                      'counted':z['counted_plants_in_main'],
                      'variance_space_percent':100*z['primary']['variance_partition']['spatial_fraction'],
                      'heldout_year_relative_gain_percent':100*z['primary']['relative_error_reduction_from_side'],
                      'sensitivity_space_percent':100*z['sensitivity_min100']['variance_partition']['spatial_fraction'],
                      'status':'SOURCE_EXPOSED_DESCRIPTIVE_REANALYSIS_NO_CAUSAL_CLAIM'},ensure_ascii=False))

if __name__=='__main__':cli()
