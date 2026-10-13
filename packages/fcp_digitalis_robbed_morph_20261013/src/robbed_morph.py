#!/usr/bin/env python3
"""Source-locked Digitalis per-(individual or patch) colour/robbed-flower comparison.

The source is a *selected morphological sample*, not a population frequency census.
All inference is conditional on five fixed surveyed sites and recorded plant/patch
clusters. A yes/no robbed-flower flag does not identify a robbing animal or
establish nectar, fitness or selection effects.
"""
from __future__ import annotations
import argparse, hashlib, json, math
from pathlib import Path
import numpy as np
import pandas as pd

EXPECTED_SIZE=345652
EXPECTED_MD5='a1200e0197b282ee126700ae234eb95e'
EXPECTED_SHA256='723c7a03b6ffa010e6093996d8203b6770ab1ee1bc6ddc21ec9c3d5b47d3b8c9'
PRIMARY_SITES=('G_1','H_1','H_3','GPS-143','GPS-160')
EXTRA_SITE='H_2'
VALID_MORPHS=frozenset(('violet','pink','white'))
PRIMARY_COLOURS=('violet','pink')
KEYS=('COUNTRY','Year','GPS-code','Plant_ID')
N_PERM=9999
N_BOOT=1999
SEED_PERM=2026101301
SEED_BOOT=2026101302
SEED_LENGTH=2026101303
SEED_DATE=2026101304


def verify_source(path: Path):
    b=Path(path).read_bytes()
    if len(b)!=EXPECTED_SIZE or hashlib.md5(b).hexdigest()!=EXPECTED_MD5:
        raise ValueError('source MD5 or file-size mismatch: do not analyse substitute source')
    if hashlib.sha256(b).hexdigest()!=EXPECTED_SHA256:
        raise ValueError('source SHA256 mismatch')
    return {'file_size_bytes':len(b),'MD5':EXPECTED_MD5,'SHA256':EXPECTED_SHA256,'encoding':'cp1252'}


def load_source(path:Path):
    receipt=verify_source(path)
    d=pd.read_csv(path,encoding='cp1252',dtype={'Plant_ID':str,'Plant':str,'Patch':str},low_memory=False)
    if len(d)!=1323 or d.shape[1]!=54:
        raise ValueError('source shape differs from frozen population')
    return d,receipt


def group_plant_units(frame: pd.DataFrame) -> pd.DataFrame:
    required=(*KEYS,'Flower_morph','Robbed_flower','Length_flower','Date')
    missing=set(required)-set(frame.columns)
    if missing:raise ValueError('missing columns: '+str(sorted(missing)))
    if frame[list(KEYS)].isna().any().any() or any(frame[k].astype(str).str.strip().eq('').any() for k in KEYS):
        raise ValueError('missing plant key')
    results=[]
    for key,g in frame.groupby(list(KEYS),sort=True,dropna=False):
        morphset=set(g['Flower_morph'].dropna().astype(str).str.strip().str.lower())-{''}
        if len(morphset)>1:
            morph='CONFLICT';error='CONFLICTING_MORPH'
        elif not morphset:
            morph='MISSING';error='MISSING_MORPH'
        else:
            morph=next(iter(morphset));error='' if morph in VALID_MORPHS else 'UNKNOWN_MORPH'
        robbed=g.Robbed_flower.astype('string').str.lower().str.strip()
        yes=int(robbed.eq('yes').sum());no=int(robbed.eq('no').sum())
        if yes+no==0 and error=='':error='NO_ASSESSED_FLOWERS'
        v=pd.to_numeric(g.Length_flower,errors='coerce').dropna()
        datekeys=tuple(sorted(set(g.Date.dropna().astype(str).str.strip())-{''}))
        row=dict(zip(KEYS,key))
        row.update(morph=morph,exclusion_reason=error,n_flower_rows=int(len(g)),n_yes=yes,n_no=no,
                   n_unassessed_flower_rows=int(len(g)-yes-no),n_assessed_flower_rows=yes+no,
                   robbed_fraction=float(yes/(yes+no)) if yes+no else float('nan'),
                   length_plant_mean=float(v.mean()) if len(v) else float('nan'),n_length=int(len(v)),
                   n_dates=len(datekeys),single_date=datekeys[0] if len(datekeys)==1 else '',
                   days_seen=';'.join(datekeys))
        results.append(row)
    out=pd.DataFrame(results)
    if out.duplicated(list(KEYS)).any():raise ValueError('duplicate grouped plant/patch keys')
    if int(out.n_flower_rows.sum())!=len(frame):raise ValueError('original flower rows not conserved')
    return out


def select_eligible(plants,sites,endpoint='robbed_fraction',morphs=PRIMARY_COLOURS):
    x=plants[(plants.Year.astype(str)=='2024') & (plants['GPS-code'].isin(sites))
             & (plants.morph.isin(morphs)) & plants.exclusion_reason.eq('')].copy()
    return x.loc[np.isfinite(pd.to_numeric(x[endpoint],errors='coerce'))]


def site_contrasts(plants,sites,min_each=5,endpoint='robbed_fraction'):
    x=select_eligible(plants,sites,endpoint)
    records=[]
    for site in sites:
        rows=x[x['GPS-code']==site]
        v=rows.loc[rows.morph=='violet',endpoint].to_numpy(float)
        p=rows.loc[rows.morph=='pink',endpoint].to_numpy(float)
        if len(v)<min_each or len(p)<min_each:
            raise ValueError(f'HOLD_INSUFFICIENT_MORPH_SUPPORT {site} violet={len(v)} pink={len(p)} min={min_each}')
        records.append({'site':site,'country':str(rows.COUNTRY.iloc[0]),
                        'n_violet_groups':int(len(v)),'n_pink_groups':int(len(p)),
                        'violet_mean':float(v.mean()),'pink_mean':float(p.mean()),
                        'difference_violet_minus_pink':float(v.mean()-p.mean()),
                        'n_inspected_flower_rows':int(rows.n_assessed_flower_rows.sum())})
    return records


def _strata_arrays(plants,sites,min_each,endpoint):
    x=select_eligible(plants,sites,endpoint)
    out=[]
    for site in sites:
        d=x[x['GPS-code']==site]
        v=d.loc[d.morph=='violet',endpoint].to_numpy(float)
        p=d.loc[d.morph=='pink',endpoint].to_numpy(float)
        if len(v)<min_each or len(p)<min_each:
            raise ValueError(f'HOLD_INSUFFICIENT_MORPH_SUPPORT {site}')
        out.append((site,v,p))
    return out


def permutation_main(plants,sites,reps=N_PERM,seed=SEED_PERM,min_each=5,endpoint='robbed_fraction'):
    strata=_strata_arrays(plants,sites,min_each,endpoint)
    point=float(np.mean([v.mean()-p.mean() for _,v,p in strata]))
    rng=np.random.default_rng(seed); draws=np.empty(reps,float)
    for k in range(reps):
        terms=[]
        for _,v,p in strata:
            allv=np.concatenate([v,p]); n=len(v)
            ix=rng.permutation(len(allv))
            terms.append(allv[ix[:n]].mean()-allv[ix[n:]].mean())
        draws[k]=np.mean(terms)
    exceed=int(np.count_nonzero(abs(draws)>=abs(point)-1e-12))
    return {'site_equal_estimate':point,'p_two_sided':(exceed+1)/(reps+1),
            'n_permutations':reps,'permutation_mean':float(draws.mean()),
            'permutation_sd':float(draws.std(ddof=1)) if reps>1 else 0,
            'abs_null_q95':float(np.quantile(abs(draws),0.95)),
            'n_sites':len(strata),'n_eligible_plant_patch_units':int(sum(len(v)+len(p) for _,v,p in strata)),
            'seed':seed,'test_scope':'observed sites and sampled plant/patch units; no causal assignment'}


def bootstrap_main(plants,sites,reps=N_BOOT,seed=SEED_BOOT,min_each=5,endpoint='robbed_fraction'):
    strata=_strata_arrays(plants,sites,min_each,endpoint)
    rng=np.random.default_rng(seed);out=np.empty(reps,float)
    for k in range(reps):
        out[k]=float(np.mean([rng.choice(v,size=len(v),replace=True).mean()-rng.choice(p,size=len(p),replace=True).mean() for _,v,p in strata]))
    return {'bootstrap_ci95':np.quantile(out,[.025,.975]).tolist(), 'n_bootstraps':reps,'seed':seed,
            'resampling_unit':'observed plant/patch groups within site and morph; five study sites fixed'}


def date_overlap_sensitivity(plants,sites,min_sites=3,min_plants=50,reps=N_PERM,seed=SEED_DATE):
    x=select_eligible(plants,sites).copy()
    x=x.loc[x.n_dates.eq(1)].copy()
    parts=[];site_count={}
    for (site,date),g in x.groupby(['GPS-code','single_date'],sort=True):
        if site not in sites:continue
        v=g.loc[g.morph=='violet','robbed_fraction'].to_numpy(float)
        p=g.loc[g.morph=='pink','robbed_fraction'].to_numpy(float)
        if len(v)>=2 and len(p)>=2:
            parts.append((site,date,v,p));site_count[site]=site_count.get(site,0)+len(v)+len(p)
    supported=[s for s in sites if s in site_count]
    n_total=sum(site_count.values())
    record={'n_supported_sites':len(supported),'n_supported_groups':n_total,
            'site_day_strata':[{ 'site':s,'date':day,'n_violet':len(v),'n_pink':len(p)} for s,day,v,p in parts],
            'excluded_mixed_date_group_count':int(select_eligible(plants,sites).n_dates.gt(1).sum()),
            'min_sites_required':min_sites,'min_plants_required':min_plants}
    if len(supported)<min_sites or n_total<min_plants:
        return record|{'status':'HOLD_INSUFFICIENT_DATE_OVERLAP'}
    def point_of_groups(group_terms):
        return float(np.mean([np.mean(group_terms[s]) for s in supported]))
    observed={s:[] for s in supported}
    for site,date,v,p in parts: observed[site].append(float(v.mean()-p.mean()))
    point=point_of_groups(observed)
    rng=np.random.default_rng(seed)
    sims=[]
    for _ in range(reps):
        tmp={s:[] for s in supported}
        for site,date,v,p in parts:
            allv=np.r_[v,p];ind=rng.permutation(len(allv));n=len(v)
            tmp[site].append(float(allv[ind[:n]].mean()-allv[ind[n:]].mean()))
        sims.append(point_of_groups(tmp))
    exceed=sum(abs(v)>=abs(point)-1e-12 for v in sims)
    return record|{'status':'EXPLORATORY_DATE_STRATIFIED','site_equal_estimate':point,
                   'p_two_sided':(exceed+1)/(reps+1),'n_permutations':reps,'seed':seed,
                   'limitation':'Complete plant groups only within overlapping site-day strata; selected time strata are not a random field sample'}



def descriptive_variance_partition(plants,sites):
    """Plant-weighted ANOVA arithmetic; not a mechanism, CI or site-random test."""
    a=plants[(plants.Year.astype(str)=='2024') & plants['GPS-code'].isin(sites)
             & plants.exclusion_reason.eq('') & np.isfinite(plants.robbed_fraction)].copy()
    if len(a)<2:raise ValueError('too few source units for descriptive variance')
    grand=float(a.robbed_fraction.mean())
    total=float(((a.robbed_fraction-grand)**2).mean())
    if total<=0:return {'status':'HOLD_ZERO_TOTAL_VARIANCE','n_plant_patch_units':len(a)}
    site=a.groupby('GPS-code').robbed_fraction.agg(['size','mean'])
    site_between=sum(float(row['size'])*(float(row['mean'])-grand)**2 for _,row in site.iterrows())/len(a)
    combined=a.groupby(['GPS-code','morph']).robbed_fraction.agg(['size','mean']).reset_index()
    combined=combined.merge(site['mean'].rename('site_mean'),left_on='GPS-code',right_index=True)
    morph_within=float((combined['size']*(combined['mean']-combined['site_mean'])**2).sum()/len(a))
    remain=total-site_between-morph_within
    if remain < -1e-12 or abs((site_between+morph_within+remain)-total)>1e-12:
        raise ValueError('within and between variance failed to reassemble')
    return {'status':'DESCRIPTIVE_POSTHOC','n_plant_patch_units':len(a),'total_plant_robbed_fraction_variance':total,
            'site_component_fraction':site_between/total,
            'morph_given_site_component_fraction':morph_within/total,
            'within_site_morph_residual_fraction':remain/total,
            'unit':'One observed plant/patch per row, site frequencies weighted by number of observed groups',
            'hard_boundary':'Variance decomposition of recorded flags only; not explanatory causal R2, not generalization to sites or regions'}


def run(source_path:Path,outdir:Path,reps=N_PERM,boot=N_BOOT):
    raw,source=load_source(source_path)
    plants=group_plant_units(raw)
    outdir.mkdir(parents=True,exist_ok=True)
    plants.to_csv(outdir/'plant_patch_ledger.csv',index=False,float_format='%.12g')
    sites=list(PRIMARY_SITES)
    descriptive=plants.groupby(['Year','GPS-code','exclusion_reason'],dropna=False).agg(n_grouped=('Plant_ID','size'),n_flower_rows=('n_flower_rows','sum')).reset_index().to_dict('records')
    main_sites=site_contrasts(plants,sites)
    primary=permutation_main(plants,sites,reps=reps)
    primary.update(bootstrap_main(plants,sites,reps=boot))
    day=date_overlap_sensitivity(plants,sites,reps=reps)
    len_plants=plants[plants.exclusion_reason.eq('')&plants.n_length.gt(0)].copy()
    length_sites=site_contrasts(len_plants,sites,endpoint='length_plant_mean')
    length=permutation_main(len_plants,sites,reps=1999,seed=SEED_LENGTH,endpoint='length_plant_mean')
    extra=site_contrasts(plants,[EXTRA_SITE],min_each=5)
    observed=plants[(plants.Year.astype(str)=='2024')&plants['GPS-code'].isin(sites)&plants.exclusion_reason.eq('')]
    variance=descriptive_variance_partition(plants,sites)
    site_colour=observed.groupby(['GPS-code','morph']).agg(n_grouped=('Plant_ID','size'),
              n_assessed=('n_assessed_flower_rows','sum'),n_yes=('n_yes','sum'),n_no=('n_no','sum'),
              group_mean_robbed_fraction=('robbed_fraction','mean')).reset_index().to_dict('records')
    all_sites=observed.groupby(['GPS-code','COUNTRY']).agg(n=('Plant_ID','size'),n_yes=('n_yes','sum'),n_no=('n_no','sum'),
                                                  mean_group_rob=('robbed_fraction','mean')).reset_index().to_dict('records')
    status={'source':source,'n_source_flower_rows':len(raw),'n_2023_source_rows':int(raw.Year.eq(2023).sum()),
            'n_2024_source_rows':int(raw.Year.eq(2024).sum()),'n_2024_primary_site_flower_rows':int(raw[raw.Year.eq(2024)&raw['GPS-code'].isin(sites)].shape[0]),
            'n_2024_H2_excluded_flower_rows':int(raw[raw.Year.eq(2024)&raw['GPS-code'].eq('H_2')].shape[0]),
            'n_grouped_units_all_years':len(plants),
            'n_grouped_2024_five_sites':int(plants.Year.eq(2024).mul(plants['GPS-code'].isin(sites)).sum()),
            'n_conflicting_morph_groups_2024_main':int(((plants.Year==2024)&(plants['GPS-code'].isin(sites))&plants.exclusion_reason.eq('CONFLICTING_MORPH')).sum()),
            'n_no_morph_groups_2024_main':int(((plants.Year==2024)&(plants['GPS-code'].isin(sites))&plants.exclusion_reason.eq('MISSING_MORPH')).sum()),
            'n_no_robbed_flag_groups_2024_main':int(((plants.Year==2024)&(plants['GPS-code'].isin(sites))&plants.exclusion_reason.eq('NO_ASSESSED_FLOWERS')).sum()),
            'n_non_yes_no_robbery_rows_2024_main':int((~raw.loc[raw.Year.eq(2024)&raw['GPS-code'].isin(sites),'Robbed_flower'].astype(str).str.lower().isin(['yes','no'])).sum())}
    output={'schema':'fcp_digitalis_fixed_pop_morph_robbed_flower_20261013_v1',
         'date_jst':'2026-10-13','protocol_path':'docs/FCP_DIGITALIS_ROBBED_FLOWER_MORPH_PROTOCOL_20261013.md',
         'original_article_doi':'10.1002/ajb2.70186','source_figshare_doi':'10.6084/m9.figshare.29095451.v2',
         'source_file_id':54632993,'scope':'2024 source-selected individual/patch measurements from five fixed study sites, NOT an ecological census or new sampling',
         'quality':status,'all_grouping_counts':descriptive,
         'primary_endpoint':'site-equal mean of (violet minus pink) mean robbed-flower yes fraction per surveyed plant/patch',
         'main_site_contrasts':main_sites,'main_paired_colour_reference':primary,
         'date_blocked_sensitivity':day,'descriptive_posthoc_variance_partition':variance,'secondary_flower_length_mm':{'site_contrasts':length_sites,'within_site_reference':length},
         'extra_H2_excluded_site':extra[0], 'all_colour_descriptive_site_summary':site_colour,
         'all_site_observation_context':all_sites,
         'explicit_limits':['The source is a selected phenotypic survey, not a randomly sampled flower-color population census.',
          'Individual/patch source units are not necessarily biological genets, especially in G_1.',
          'Some site/year observation protocols differ; excluded 2023 source observations cannot be used for out-of-year tests.',
          'Within-site permutation is a conditional reference assuming exchangeability of observed morph labels; not causal treatment assignment.',
          'Five source sites fixed, not a random sample of biogeographic regions; country comparison only descriptive and confounded with survey dates.',
          'Robbed_flower is a source binary observational flag; no verified insect culprit, nectar reward or fitness costs.',
          'Colour conflict groups excluded; missing robbed flags not converted to negatives.',
          'Source author reported no consistent flower morphology differences among colour morphs; this does not establish a new general mechanism.',
          'No genetic polymorphism mechanisms, local adaptation, pollinator preference or population frequency inferred.'],
         'previous_main_results_changed':False}
    (outdir/'result.json').write_text(json.dumps(output,indent=2,ensure_ascii=False)+'\n')
    pd.DataFrame(main_sites).to_csv(outdir/'site_primary_contrasts.csv',index=False)
    return output


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--outdir',type=Path,required=True)
    p.add_argument('--permutations',type=int,default=N_PERM);p.add_argument('--bootstraps',type=int,default=N_BOOT)
    a=p.parse_args()
    result=run(a.source,a.outdir,reps=a.permutations,boot=a.bootstraps)
    print(json.dumps({'n_sites':len(result['main_site_contrasts']),
          'site_contrasts':result['main_site_contrasts'],
          'primary':result['main_paired_colour_reference'],
          'date':result['date_blocked_sensitivity'],
          'secondary_length':result['secondary_flower_length_mm']['within_site_reference'],
          'quality':result['quality']},ensure_ascii=False),flush=True)
