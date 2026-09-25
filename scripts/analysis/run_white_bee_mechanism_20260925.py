#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, re
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import spearmanr

MORPHS=['white','yellow_orange','red_pink','blue_purple']

def norm(x):
    return re.sub(r'\s+',' ',str(x).strip()).lower()

def bool_series(s):
    if s.dtype==bool: return s.fillna(False)
    return s.astype(str).str.strip().str.lower().isin({'true','1','yes','y'})

def z(s):
    x=pd.to_numeric(s,errors='coerce').astype(float)
    sd=float(x.std(ddof=0))
    if not np.isfinite(sd) or sd<=0:
        return pd.Series(np.nan,index=x.index)
    return (x-float(x.mean()))/sd

def build_bee_metrics(path,chunksize=150000):
    records=defaultdict(int)
    bees=defaultdict(set)
    dbs=defaultdict(set)
    family={}
    cols=['plant_species','bee_species','plant_family','database']
    for chunk in pd.read_csv(path,usecols=cols,chunksize=chunksize,low_memory=False,
                             encoding='utf-8',encoding_errors='replace'):
        chunk=chunk.dropna(subset=['plant_species','bee_species'])
        for r in chunk.itertuples(index=False):
            k=norm(r.plant_species)
            if not k: continue
            records[k]+=1
            bees[k].add(str(r.bee_species).strip())
            if pd.notna(r.database): dbs[k].add(str(r.database).strip())
            if k not in family and pd.notna(r.plant_family):
                family[k]=str(r.plant_family).strip()
    rows=[]
    for k,n in records.items():
        rows.append({
            'species_key':k,
            'n_bee_records':int(n),
            'n_bee_species':int(len(bees[k])),
            'n_bee_databases':int(len(dbs[k])),
            'plant_family':family.get(k,'unknown')
        })
    m=pd.DataFrame(rows)
    cal=m[(m.n_bee_records>=5)&(m.n_bee_species>0)].copy()
    cal['log_records']=np.log1p(cal.n_bee_records)
    cal['log_databases']=np.log1p(cal.n_bee_databases)
    cal['log_bees']=np.log1p(cal.n_bee_species)
    X=sm.add_constant(cal[['log_records','log_databases']],has_constant='add')
    fit=sm.OLS(cal.log_bees,X).fit()
    cal['bee_partner_breadth_resid']=cal.log_bees-fit.predict(X)
    mu=float(cal.bee_partner_breadth_resid.mean())
    sd=float(cal.bee_partner_breadth_resid.std(ddof=0))
    cal['bee_partner_breadth_resid_z']=(cal.bee_partner_breadth_resid-mu)/sd
    m=m.merge(cal[['species_key','bee_partner_breadth_resid_z']],on='species_key',how='left',validate='one_to_one')
    calibration={
      'n_plants':int(len(cal)),
      'formula':'log1p(n_bee_species) ~ log1p(n_bee_records) + log1p(n_bee_databases)',
      'params':{k:float(v) for k,v in fit.params.items()},
      'r_squared':float(fit.rsquared),
      'residual_mean':mu,
      'residual_sd':sd
    }
    return m,calibration

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--measured',required=True)
    p.add_argument('--technical-table',required=True)
    p.add_argument('--join-key',required=True)
    p.add_argument('--high-clip-ids',required=True)
    p.add_argument('--bee-curated',required=True)
    p.add_argument('--outdir',required=True)
    p.add_argument('--permutations',type=int,default=9999)
    p.add_argument('--seed',type=int,default=20260925)
    a=p.parse_args()
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)

    bee,calibration=build_bee_metrics(a.bee_curated)

    use=['photo_id','species','morph','global_classifiable']
    d=pd.read_csv(a.measured,usecols=use)
    d=d[bool_series(d.global_classifiable)&d.morph.isin(MORPHS)].copy()
    d['white']=(d.morph=='white').astype(int)

    tech=pd.read_csv(a.technical_table,compression='gzip',dtype={'measurement_id':str})
    join=pd.read_csv(a.join_key,dtype={'measurement_id':str})
    high=pd.read_csv(a.high_clip_ids,dtype={'measurement_id':str})
    tech=tech.merge(join,on='measurement_id',how='left',validate='one_to_one')
    tech['high_clip']=tech.measurement_id.astype(str).isin(set(high.measurement_id.astype(str)))
    d=d.merge(tech[['photo_id','near_clip_fraction','high_clip']],on='photo_id',how='left',validate='one_to_one')
    d=d[d.near_clip_fraction.notna() & ~d.high_clip.fillna(False)].copy()

    sp=d.groupby('species').white.agg(['sum','count']).reset_index()
    sp['n_white']=sp['sum'].astype(int)
    sp['n_nonwhite']=(sp['count']-sp['sum']).astype(int)
    sp['white_fraction']=sp.n_white/sp['count']
    sp['species_key']=sp.species.map(norm)
    sp=sp[(sp['count']>=40)&(sp.n_white>=5)&(sp.n_nonwhite>=5)].copy()
    x=sp.merge(bee,on='species_key',how='inner',validate='one_to_one')
    x=x[(x.n_bee_records>=5)&x.bee_partner_breadth_resid_z.notna()].copy()

    if len(x)<100:
        result={'schema':'fcp_white_bee_mechanism_v1','status':'not_estimable',
                'eligible_species':int(len(x)),'minimum_species':100,
                'reason':'fewer than 100 species after frozen colour and bee-metric gates'}
        (out/'result.json').write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result,indent=2)); return

    rho=float(spearmanr(x.bee_partner_breadth_resid_z,x.white_fraction).statistic)
    raw_rho=float(spearmanr(np.log1p(x.n_bee_species),x.white_fraction).statistic)
    effort_rho=float(spearmanr(x.bee_partner_breadth_resid_z,np.log1p(x.n_bee_records)).statistic)

    rng=np.random.default_rng(a.seed)
    obs=x.white_fraction.to_numpy(float)
    pred=x.bee_partner_breadth_resid_z.to_numpy(float)
    null=np.empty(a.permutations,float)
    for i in range(a.permutations):
        null[i]=spearmanr(rng.permutation(pred),obs).statistic
    perm_p=float((1+np.sum(null<=rho))/(a.permutations+1))

    x['log_records_z']=z(np.log1p(x.n_bee_records))
    x['log_databases_z']=z(np.log1p(x.n_bee_databases))
    X=sm.add_constant(x[['bee_partner_breadth_resid_z','log_records_z','log_databases_z']],has_constant='add')
    fit=sm.GLM(x.white_fraction,X,family=sm.families.Binomial()).fit(
        cov_type='cluster',cov_kwds={'groups':x.plant_family.astype(str)})
    beta=float(fit.params['bee_partner_breadth_resid_z'])
    p_cluster=float(fit.pvalues['bee_partner_breadth_resid_z'])
    se=float(fit.bse['bee_partner_breadth_resid_z'])

    lofo=[]
    for fam in sorted(x.plant_family.astype(str).unique()):
        q=x[x.plant_family.astype(str)!=fam].copy()
        if len(q)<100: continue
        XX=sm.add_constant(q[['bee_partner_breadth_resid_z','log_records_z','log_databases_z']],has_constant='add')
        try:
            ff=sm.GLM(q.white_fraction,XX,family=sm.families.Binomial()).fit()
            lofo.append(float(ff.params['bee_partner_breadth_resid_z']))
        except Exception:
            pass

    gate=bool(rho<0 and perm_p<0.05 and beta<0 and p_cluster<0.05)
    x.sort_values('species').to_csv(out/'species_bee_white_analysis.csv',index=False)
    pd.DataFrame({'rho_null':null}).to_csv(out/'spearman_permutation_null.csv',index=False)
    result={
      'schema':'fcp_white_bee_mechanism_v1',
      'status':'complete',
      'role':'prospectively_frozen_after_outcome_blind_coverage_before_bee_colour_opening',
      'eligible_species':int(len(x)),
      'represented_families':int(x.plant_family.nunique()),
      'bee_metric_calibration':calibration,
      'primary_prediction':'greater effort-adjusted bee partner breadth -> lower white fraction',
      'primary_spearman_rho':rho,
      'primary_one_sided_permutation_p':perm_p,
      'permutations':int(a.permutations),
      'raw_log_bee_richness_rho':raw_rho,
      'bee_metric_vs_log_records_rho':effort_rho,
      'fractional_glm':{
        'beta_bee_breadth':beta,
        'se_cluster':se,
        'p_family_cluster':p_cluster,
        'odds_multiplier_per_residual_SD':float(np.exp(beta)),
        'ci_low':float(np.exp(beta-1.96*se)),
        'ci_high':float(np.exp(beta+1.96*se))
      },
      'leave_one_family_out':{
        'valid_refits':len(lofo),
        'beta_min':float(np.min(lofo)) if lofo else None,
        'beta_max':float(np.max(lofo)) if lofo else None,
        'fraction_negative':float(np.mean(np.array(lofo)<0)) if lofo else None
      },
      'mechanism_gate_pass':gate,
      'verdict':'BEE_BREADTH_SUPPORTED' if gate else 'BEE_BREADTH_NOT_SUPPORTED_UNDER_THIS_TEST',
      'hard_nonclaims':[
        'does not establish pigmented-to-white evolutionary transition direction',
        'does not measure morph-specific bee preference or pollination effectiveness',
        'does not test moth, butterfly, fly, bird, or bat mechanisms',
        'cross-species association cannot by itself establish causality'
      ]
    }
    (out/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
