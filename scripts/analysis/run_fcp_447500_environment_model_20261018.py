#!/usr/bin/env python3
"""Two-stage Firth logistic / REML model of within-species photo colour.

This is an approximate two-stage analysis, NOT a one-stage GLMM. Responses are
white=1 versus the other three classified morphs=0. No palette-derived exposure
control is used. Seven-variable complete-case selection matches the earlier panel.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import json
import hashlib
from pathlib import Path
import platform
import sys

import numpy as np
import pandas as pd
import scipy
from scipy.optimize import minimize, minimize_scalar, brentq
from scipy.special import expit
from scipy.stats import norm, t, chi2

from scripts.analysis.summarize_fcp_447500_within_species_environment_20261018 import (
    COHORTS, MORPHS, VARIABLES, PHOTO_SHA, ZIP_SHA, demand, sha,
    read_and_sample, json_safe)

PREDICTORS=('BIO5','BIO14','SRADmean')
METHOD_SOURCES=[
 'https://search.r-project.org/CRAN/refmans/logistf/html/logistf-package.html',
 'https://wviechtb.github.io/metafor/reference/misc-recs.html',
 'https://wviechtb.github.io/metafor/reference/profile.html']


def firth_objective(beta, x, y):
    """Negative Jeffreys-penalized binomial log likelihood and analytic gradient."""
    eta=x@beta; p=expit(eta); w=p*(1-p)
    info=x.T@(w[:,None]*x)
    sign,logdet=np.linalg.slogdet(info)
    if sign<=0: return np.inf, np.full(x.shape[1],np.nan)
    try: inv=np.linalg.inv(info)
    except np.linalg.LinAlgError: return np.inf, np.full(x.shape[1],np.nan)
    leverage=w*np.einsum('ij,jk,ik->i',x,inv,x)
    ll=np.sum(y*eta-np.logaddexp(0,eta))+.5*logdet
    score=x.T@(y-p+leverage*(.5-p))
    return -float(ll),-score


def fit_firth(x,y):
    """Fit with explicit rank/convergence checks; no automatic variable dropping."""
    x=np.asarray(x,float); y=np.asarray(y,float)
    if np.linalg.matrix_rank(x)<x.shape[1]: return {'status':'rank_deficient'}
    condition=float(np.linalg.cond(x))
    if not np.isfinite(condition) or condition>1e8:
        return {'status':'ill_conditioned_design','design_condition':condition}
    start=np.zeros(x.shape[1]); start[0]=np.log((y.sum()+.5)/(len(y)-y.sum()+.5))
    result=minimize(firth_objective,start,args=(x,y),jac=True,method='BFGS',
                    options={'gtol':1e-7,'maxiter':400})
    objective,gradient=firth_objective(result.x,x,y)
    # BFGS can report precision loss when the score already satisfies tolerance.
    converged=bool(np.isfinite(objective) and np.isfinite(result.x).all()
                   and np.max(np.abs(gradient))<1e-5)
    p=expit(x@result.x); w=p*(1-p); info=x.T@(w[:,None]*x)
    try: cov=np.linalg.inv(info)
    except np.linalg.LinAlgError: return {'status':'singular_information'}
    eigen=np.linalg.eigvalsh(cov)
    if not converged or np.any(eigen<=0) or not np.isfinite(cov).all():
        return {'status':'nonconverged','max_abs_score':float(np.max(np.abs(gradient))),
                'optimizer_message':str(result.message),'iterations':int(result.nit)}
    h=w*np.einsum('ij,jk,ik->i',x,cov,x)
    return {'status':'ok','beta':result.x,'cov':cov,'design_condition':condition,
            'information_condition':float(np.linalg.cond(info)),
            'max_leverage':float(h.max()),'high_leverage_gt_0p8':bool(h.max()>.8),'max_abs_score':float(np.max(np.abs(gradient))),
            'iterations':int(result.nit),'optimizer_success':bool(result.success),
            'optimizer_message':str(result.message),'penalized_log_likelihood':-objective,
            'min_fitted_probability':float(p.min()),'max_fitted_probability':float(p.max())}


def reml_meta(y,v):
    """Marginal normal-normal REML; modified HKSJ mean and profile tau2 CI."""
    y=np.asarray(y,float); v=np.asarray(v,float); k=len(y)
    if k<3 or not np.isfinite(y).all() or not (np.isfinite(v)&(v>0)).all():
        raise ValueError('meta-analysis needs >=3 finite estimates and positive variances')
    def objective(tau2):
        vv=v+tau2; w=1/vv; mu=np.dot(w,y)/w.sum()
        return .5*(np.log(vv).sum()+np.log(w.sum())+np.dot(w,(y-mu)**2))
    upper=max(1.,float(np.var(y))*4)
    opt=minimize_scalar(objective,bounds=(0,upper),method='bounded',options={'xatol':1e-10})
    while opt.x>.95*upper:
        upper*=4
        if upper>1e8: raise ValueError('REML optimum unbounded')
        opt=minimize_scalar(objective,bounds=(0,upper),method='bounded',options={'xatol':1e-10})
    tau2=float(opt.x) if objective(opt.x)<objective(0) else 0.
    w=1/(v+tau2); mu=float(np.dot(w,y)/w.sum()); q=float(np.dot(w,(y-mu)**2)/(k-1))
    muvar=max(1.,q)/w.sum(); crit=t.ppf(.975,k-1)
    cutoff=objective(tau2)+chi2.ppf(.95,1)/2
    lo=0. if objective(0)<=cutoff else brentq(lambda z:objective(z)-cutoff,0,tau2)
    hi=max(upper,tau2*2+1)
    while objective(hi)<cutoff: hi*=2
    hi=brentq(lambda z:objective(z)-cutoff,tau2,hi)
    eb_weight=tau2/(tau2+v)
    eb_mean=mu+eb_weight*(y-mu)
    eb_var=tau2*v/(tau2+v)
    # Conditional EB intervals deliberately exclude estimated hyperparameter uncertainty.
    sign_prob=norm.cdf(eb_mean/np.sqrt(eb_var)) if tau2>0 else np.full(k,np.nan)
    fixed_equal=float(y.mean()); fixed_se=float(np.sqrt(v.sum())/k)
    pred=t.ppf(.975,k-2)*np.sqrt(tau2+muvar)
    out={'n_species':k,'reml_mean':mu,'reml_mean_se_mhksj':float(np.sqrt(muvar)),
         'reml_mean_ci_low':float(mu-crit*np.sqrt(muvar)),
         'reml_mean_ci_high':float(mu+crit*np.sqrt(muvar)),
         'tau2':tau2,'tau':float(np.sqrt(tau2)),
         'tau2_profile_ci_low':float(lo),'tau2_profile_ci_high':float(hi),
         'working_new_species_prediction_low':float(mu-pred),
         'working_new_species_prediction_high':float(mu+pred),
         'hksj_q':q,'hksj_q_used':max(1.,q),
         'equal_species_unshrunk_mean':fixed_equal,
         'equal_species_unshrunk_se_independence':fixed_se,
         'equal_species_unshrunk_ci_low':fixed_equal-norm.ppf(.975)*fixed_se,
         'equal_species_unshrunk_ci_high':fixed_equal+norm.ppf(.975)*fixed_se,
         'equal_species_eb_mean_conditional':float(eb_mean.mean()),
         'selected_panel_mean_eb_positive_probability_conditional':float(sign_prob.mean()) if tau2>0 else None,
         'working_normal_latent_positive_fraction':float(norm.cdf(mu/np.sqrt(tau2))) if tau2>0 else None,
         'observed_estimate_positive_fraction':float(np.mean(y>0)),
         'reml_optimizer_success':bool(opt.success)}
    return out,eb_mean,eb_var,sign_prob


def eligible(g,minimum_each=5):
    nw=int(g.white.sum()); return len(g)>=40 and min(nw,len(g)-nw)>=minimum_each


def thin_photos(g,mode):
    """Response-blind deterministic thinning, with missing observer IDs excluded."""
    ordered=g.assign(_thin_order=g.photo_id.map(lambda v: hashlib.sha256(('20261018:'+str(v)).encode()).hexdigest())).sort_values('_thin_order',kind='stable')
    if mode=='observer_thinned':
        ordered=ordered.loc[ordered.observer_id.notna()].copy()
        return ordered.drop_duplicates('observer_id')
    if mode=='spatial_thinned':
        ordered=ordered.copy()
        ordered['_latcell']=np.floor(ordered.latitude).astype(int)
        ordered['_loncell']=np.floor(ordered.longitude).astype(int)
        return ordered.drop_duplicates(['_latcell','_loncell'])
    return ordered


def run_models(panel):
    records=[]; diagnostics=[]; selection=[]
    for index,(tid,base) in enumerate(panel.groupby('inat_taxon_id',sort=True)):
        native=base[list(PREDICTORS)].to_numpy(float)
        scales=native.std(axis=0,ddof=0)
        for mode in ('primary','observer_thinned','spatial_thinned','geographic_trend'):
            g=thin_photos(base,mode)
            selected=eligible(g)
            entry={'inat_taxon_id':int(tid),'species':str(base.species.iloc[0]),
                   'source_cohort':str(base.fcp_source_cohort.iloc[0]),'model':mode,
                   'n_primary':len(base),'n_retained':len(g),'n_white':int(g.white.sum()),
                   'n_nonwhite':int((~g.white).sum()),'eligible_after_thinning':selected,
                   'n_primary_observers':int(base.observer_id.nunique()),
                   'n_primary_missing_observer':int(base.observer_id.isna().sum())}
            if not selected:
                entry['status']='ineligible_after_thinning'; selection.append(entry); continue
            threshold=32*np.finfo(float).eps*np.maximum(1,np.abs(native).max(axis=0))
            if np.any(scales<=threshold):
                entry['status']='zero_or_near_zero_primary_predictor_sd'; selection.append(entry); continue
            z=(g[list(PREDICTORS)].to_numpy(float)-g[list(PREDICTORS)].mean().to_numpy())/scales
            x=np.column_stack([np.ones(len(g)),z])
            if mode=='geographic_trend':
                radians=np.deg2rad(g.longitude.to_numpy(float))
                resultant=np.abs(np.exp(1j*radians).mean())
                origin=np.angle(np.exp(1j*radians).mean())
                delta=np.rad2deg(np.angle(np.exp(1j*(radians-origin))))
                geography=np.column_stack([g.latitude.to_numpy(float),delta])
                geographic_sd=geography.std(axis=0)
                entry['longitude_circular_resultant']=float(resultant)
                entry['unwrapped_longitude_range_degrees']=float(np.ptp(delta))
                if resultant<.1 or np.ptp(delta)>180 or np.any(geographic_sd<1e-8):
                    entry['status']='geographic_basis_not_identifiable'; selection.append(entry); continue
                geography=(geography-geography.mean(axis=0))/geographic_sd
                x=np.column_stack([x,geography])
            result=fit_firth(x,g.white.to_numpy(float))
            entry['status']=result['status']; selection.append(entry)
            diagnostics.append({**entry,**{k:v for k,v in result.items() if k not in ('beta','cov')}})
            if result['status']!='ok': continue
            beta=result['beta']; cov=result['cov']
            for j,variable in enumerate(PREDICTORS,1):
                se=float(np.sqrt(cov[j,j]))
                row={**entry,'variable':variable,'slope':float(beta[j]),'se_approx':se,
                     'variance_approx':float(cov[j,j]),'ci_low_wald':float(beta[j]-1.95996398454*se),
                     'ci_high_wald':float(beta[j]+1.95996398454*se),
                     'primary_within_species_sd_native':float(scales[j-1]),
                     'slope_per_native_unit':float(beta[j]/scales[j-1]),'intercept':float(beta[0]),
                     'primary_predictor_min':float(native[:,j-1].min()),
                     'primary_predictor_max':float(native[:,j-1].max())}
                for jj,vv in enumerate(PREDICTORS,1): row['cov_with_'+vv]=float(cov[j,jj])
                records.append(row)
        if (index+1)%250==0: print('FIT_SPECIES',index+1,flush=True)
    return pd.DataFrame(records),pd.DataFrame(diagnostics),pd.DataFrame(selection)


def summarize(slopes):
    summaries=[]; eb=[]
    for model,model_rows in slopes.groupby('model',sort=False):
        for cohort in ('pooled',*COHORTS):
            subset=model_rows if cohort=='pooled' else model_rows.loc[model_rows.source_cohort.eq(cohort)]
            for variable,g in subset.groupby('variable',sort=False):
                if len(g)<3: continue
                result,m,v,prob=reml_meta(g.slope,g.variance_approx)
                summaries.append({'model':model,'group':cohort,'variable':variable,**result})
                if cohort=='pooled':
                    for tid,mean,var,pr in zip(g.inat_taxon_id,m,v,prob):
                        eb.append({'model':model,'variable':variable,'inat_taxon_id':int(tid),
                                   'eb_mean_conditional':float(mean),'eb_variance_conditional':float(var),
                                   'eb_positive_probability_conditional':float(pr)})
    # Matched sets isolate thinning effects from the loss of eligible species.
    primary=slopes.loc[slopes.model.eq('primary')]
    comparisons=[]
    for model in ('observer_thinned','spatial_thinned','geographic_trend'):
        s=slopes.loc[slopes.model.eq(model)]
        for variable,g in s.groupby('variable',sort=False):
            b=primary.loc[primary.variable.eq(variable)].set_index('inat_taxon_id')
            g=g.set_index('inat_taxon_id'); shared=b.index.intersection(g.index)
            if len(shared)<3: continue
            b=b.loc[shared]; g=g.loc[shared]
            out,_,_,_=reml_meta(b.slope,b.variance_approx)
            summaries.append({'model':model+'_matched_primary','group':'pooled','variable':variable,**out})
            comparisons.append({'sensitivity':model,'variable':variable,'n_matched_species':len(shared),
                                'matched_primary_equal_mean':float(b.slope.mean()),
                                'thinned_equal_mean':float(g.slope.mean()),
                                'paired_equal_mean_change':float((g.slope-b.slope).mean()),
                                'paired_median_change':float((g.slope-b.slope).median()),
                                'estimate_sign_agreement':float(np.mean((g.slope>0)==(b.slope>0)))})
    return pd.DataFrame(summaries),pd.DataFrame(eb),pd.DataFrame(comparisons)


def load_panel(photos,root):
    demand(sha(photos)==PHOTO_SHA,'source photograph SHA mismatch')
    for name,digest in ZIP_SHA.items(): demand(sha(root/name)==digest,name+' SHA mismatch')
    cols=['inat_taxon_id','species','photo_id','fcp_source_cohort','morph','global_classifiable',
          'latitude','longitude','observer_id','observation_id','observed_on','positional_accuracy_m']
    full=pd.read_csv(photos,usecols=cols,low_memory=False)
    demand(len(full)==447500 and full.inat_taxon_id.nunique()==4475 and full.photo_id.is_unique,'source identities')
    classified=full.global_classifiable.fillna(False).astype(str).str.lower().isin(['true','1','yes','y'])
    demand(classified.sum()==227107 and full.loc[classified].morph.isin(MORPHS).all(),'classified source differs')
    p=full.loc[classified].copy().reset_index(drop=True); meta=[]
    for number in (1,4,5,12,14,15):
        vals,m=read_and_sample(root/'bio'/f'wc2.1_10m_bio_{number}.tif',p.latitude,p.longitude,root/'wc2.1_10m_bio.zip')
        p['BIO'+str(number)]=vals;meta.append(m)
    monthly=[]
    for month in range(1,13):
        vals,m=read_and_sample(root/'srad'/f'wc2.1_10m_srad_{month:02d}.tif',p.latitude,p.longitude,root/'wc2.1_10m_srad.zip')
        monthly.append(vals);meta.append(m)
    p['SRADmean']=np.stack(monthly,axis=1).mean(axis=1)
    cc=np.isfinite(p[list(VARIABLES)]).all(axis=1); missing=int((~cc).sum());p=p.loc[cc].copy()
    p['white']=p.morph.eq('white')
    panel=p.groupby('inat_taxon_id',group_keys=False).filter(eligible)
    demand(len(panel)==86114 and panel.inat_taxon_id.nunique()==1574,'primary complete-case panel changed')
    return panel,meta,missing


def genus_bootstrap(slopes, repetitions=1000, seed=20261018):
    # Pragmatic taxonomy-based dependence check; genus string is NOT a phylogeny.
    primary=slopes.loc[slopes.model.eq('primary')].copy()
    primary['genus']=primary.species.str.split().str[0]
    out=[]
    for variable,g in primary.groupby('variable',sort=False):
        aggregate=g.groupby('genus').slope.agg(['sum','count'])
        rng=np.random.default_rng(seed)
        draws=rng.integers(0,len(aggregate),size=(repetitions,len(aggregate)))
        values=aggregate['sum'].to_numpy()[draws].sum(axis=1)/aggregate['count'].to_numpy()[draws].sum(axis=1)
        out.append({'variable':variable,'n_species':len(g),'n_genera':len(aggregate),
                    'repetitions':repetitions,'seed':seed,'equal_species_mean':float(g.slope.mean()),
                    'genus_cluster_bootstrap_sd':float(values.std(ddof=1)),
                    'genus_cluster_percentile_low':float(np.quantile(values,.025)),
                    'genus_cluster_percentile_high':float(np.quantile(values,.975))})
    return pd.DataFrame(out)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--photos',type=Path,required=True)
    parser.add_argument('--worldclim-root',type=Path,required=True)
    parser.add_argument('--outdir',type=Path,required=True)
    args=parser.parse_args();started=datetime.now(timezone.utc).isoformat()
    panel,rasters,missing=load_panel(args.photos,args.worldclim_root)
    print('AUTHENTICATED_PANEL',panel.inat_taxon_id.nunique(),len(panel),flush=True)
    slopes,diagnostics,selection=run_models(panel)
    summary,eb,matched=summarize(slopes)
    # >=10 each uses identical first-stage fits on the predefined stricter panel.
    stricter_ids=[tid for tid,g in panel.groupby('inat_taxon_id') if eligible(g,minimum_each=10)]
    strict=slopes.loc[slopes.model.eq('primary') & slopes.inat_taxon_id.isin(stricter_ids)].copy()
    strict['model']='primary_ge10';strict_summary,strict_eb,_=summarize(strict)
    summary=pd.concat([summary,strict_summary],ignore_index=True);eb=pd.concat([eb,strict_eb],ignore_index=True)
    args.outdir.mkdir(parents=True,exist_ok=True)
    frames={'species_slopes.csv':slopes,'fit_diagnostics.csv':diagnostics,'selection_audit.csv':selection,
            'meta_summary.csv':summary,'species_empirical_bayes.csv':eb,'matched_sensitivity.csv':matched,
            'genus_cluster_bootstrap.csv':genus_bootstrap(slopes)}
    for name,frame in frames.items(): frame.to_csv(args.outdir/name,index=False)
    profile={'schema':'fcp_447500_within_species_two_stage_environment_model_v1',
      'analysis_type':'post_hoc_two_stage_approximate_firth_logistic_then_marginal_normal_REML_not_GLMM',
      'n_primary_eligible_species':1574,'n_primary_photos':len(panel),'n_source_photos':447500,
      'n_classified_photos':227107,'n_missing_seven_variable_climate':missing,
      'primary_by_cohort_species':panel.groupby('fcp_source_cohort').inat_taxon_id.nunique().to_dict(),
      'primary_unique_observation_ids':int(panel.observation_id.nunique()),
      'primary_missing_observer_ids':int(panel.observer_id.isna().sum()),
      'n_ge10_species':len(stricter_ids),'n_ge10_photos':int(panel.inat_taxon_id.isin(stricter_ids).sum()),
      'high_leverage_fit_counts':diagnostics.groupby('model').high_leverage_gt_0p8.sum().astype(int).to_dict(),
      'selection_counts':selection.groupby(['model','status']).size().unstack(fill_value=0).to_dict(orient='index'),
      'method':{'response':'white=1; yellow_orange/red_pink/blue_purple=0',
        'predictors':list(PREDICTORS),'formula':'logit P(white_si)=alpha_s + sum_k beta_sk * centered_environment_sik / SD_primary_sk',
        'primary_selection':'seven-variable complete cases; >=40 photographs and >=5 each response',
        'stage_one':'Firth Jeffreys penalized likelihood, intercept plus 3 joint slopes per species',
        'sampling_covariance':'inverse unpenalized Fisher information at Firth estimate; approximate Wald',
        'scale':'each species primary-panel population SD (ddof=0); sensitivity centered anew, same primary SD',
        'stage_two':'separate marginal normal-normal REML for each joint climate coefficient; sampling variance included',
        'meta_covariance':'three slope covariances saved; between-axis covariance not modeled',
        'hypermean_ci':'modified HKSJ: t(K-1), max(1,q)/sum(w)',
        'tau2_ci':'restricted likelihood profile, chi-square(1) cutoff; approximate at zero boundary',
        'prediction_interval':'working normal new species slope, t(K-2)*sqrt(tau2+modified_HKSJ_mean_variance)',
        'equal_species_mean':'arithmetic mean of unshrunk fitted slopes; finite-panel SE sqrt(sum sampling variances)/K under independence',
        'EB':'marginal empirical Bayes conditional on REML mu/tau; ignores hyperparameter uncertainty; sign fractions suppressed when tau2=0 because degenerate plugin probabilities ignore mean/tau uncertainty',
        'observer_sensitivity':'one minimum SHA256(seed:photo_id) per species/observer, missing observers excluded; same eligibility applied',
        'spatial_sensitivity':'one minimum SHA256(seed:photo_id) per species/1-degree lat-lon cell; same eligibility applied; not minimum separation',
        'geographic_trend':'add standardized within-species latitude and circular-unwrapped longitude; reject resultant<0.1/range>180 degrees/zero SD; fit full 6-column design with rank and condition gates',
        'genus_cluster_bootstrap':'1000 resamples of genus clusters, all species per sampled genus; arithmetic mean slopes; descriptive superpopulation stability, not fixed-panel sampling CI or full phylogenetic correction',
        'matched_comparison':'full primary fit restricted to successfully fitted sensitivity species',
        'ge10_sensitivity':'primary fits subset with >=10 white and >=10 chromatic',
        'method_sources':METHOD_SOURCES},
      'limitations':[
        'Approximate two-stage normal sampling model, not an exact one-stage hierarchical binomial model',
        'Firth Fisher/Wald sampling uncertainty is asymptotic and can be poor for sparse or highly collinear species',
        'Primary sampling variances assume independent photos; thinning reduces repeated sampling but does not eliminate spatial or observer confounding',
        'Equal-species finite-panel intervals assume independence of species; no phylogenetic covariance is included',
        'Thinning uses response-blind SHA256 seed 20261018; retained panels differ and matched baselines are provided',
        'Geographic trend is a limited linear adjustment, not a spatial covariance model; thinning is not confounding control. No seasonal/date, cultivation, measurement-quality, background or exposure adjustment; no all-cohort near_clip data',
        'Climate is 1970-2000 long-term normals, solar radiation is not direct UV, and SD units differ physically between species',
        'Latent sign probabilities and new-species prediction intervals depend on normality/exchangeability; no actual unknown-species or genus holdout validation',
        'Joint coefficients are conditional associations, not causal importance, genetic reaction norms, pigment mechanisms or fitness effects'],
      'summary':summary.to_dict(orient='records')}
    (args.outdir/'profile.json').write_text(json.dumps(json_safe(profile),indent=2,sort_keys=True,allow_nan=False)+'\n')
    receipt={'started_utc':started,'completed_utc':datetime.now(timezone.utc).isoformat(),'command':sys.argv,
       'python':sys.version,'platform':platform.platform(),'packages':{'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__},
       'photo_sha256':PHOTO_SHA,'worldclim_zip_sha256':ZIP_SHA,'rasters':rasters,
       'script_sha256':sha(Path(__file__)),
       'imported_sampling_script_sha256':sha(Path(__file__).with_name('summarize_fcp_447500_within_species_environment_20261018.py')),
       'output_sha256':{name:sha(args.outdir/name) for name in (*frames,'profile.json')}}
    (args.outdir/'execution_receipt.json').write_text(json.dumps(json_safe(receipt),indent=2,sort_keys=True,allow_nan=False)+'\n')
    print(summary.loc[summary.group.eq('pooled'),['model','variable','n_species','reml_mean','reml_mean_ci_low','reml_mean_ci_high','tau','equal_species_unshrunk_mean']].to_string(index=False),flush=True)
    print('FCP_ENVIRONMENT_TWO_STAGE_COMPLETE',flush=True)

if __name__=='__main__': main()
