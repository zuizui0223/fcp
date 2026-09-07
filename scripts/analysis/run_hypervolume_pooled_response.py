#!/usr/bin/env python3
"""Uncertainty-aware sharing, synthetic labels on frozen actual coordinates.

Specification d4ba51f5cb5dabf9080a2d58358f49191a673cf5 precedes outcomes.
No real colours, climate rasters, or true-response parameters enter predictors.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import log_ndtr, logsumexp
import run_hypervolume_actual_support_blocked as support
import run_hypervolume_support_response_diagnostic as previous
import run_hypervolume_transfer_benchmark as base

SEED = 2026090704
REPS = 250
SPEC = 'd4ba51f5cb5dabf9080a2d58358f49191a673cf5'
K = 48
THETA = np.arange(K) * np.pi / K
AXES = np.stack([np.cos(THETA), np.sin(THETA)], axis=1)
FRACTIONS = np.array([0., .25, .5, .75, 1.])
SLOPES = {'saturating_probit': np.array([0., .5, 1., 2.]),
          'linear_logit': np.array([0., .5, 1., 2., 4., 8.])}
MODELS = tuple(name+'_'+strategy for name in SLOPES
               for strategy in ('joint_train', 'single_source_average')) + ('qda_equal',)
NUISANCE = tuple(support.NUISANCE)
POSITIVE = tuple(previous.POSITIVE)
PARENT_HASHES = {
    'run_hypervolume_actual_support_blocked.py': 'dbedd48765816f6a2c8678585c62c52add32bdbd75a3fab59a66f06d8084891e',
    'run_hypervolume_transfer_benchmark.py': 'b5312f75128992f99b45898a6e5190aab1a86e18ad59573d26ca5fae66d8219b',
    'run_hypervolume_support_response_diagnostic.py': '93585ec4e1cee4aca7e316d93590694f2f9f1bbb073d6b7cf8567bcb3ea98808',
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def seed_for(*parts: object) -> int:
    return int.from_bytes(hashlib.sha256('|'.join(map(str, (SEED,)+parts)).encode()).digest()[:8], 'little')


def schedule(geo: support.Geometry, stage: str, replicate: int) -> np.ndarray:
    rng = np.random.default_rng(seed_for('schedule', stage, replicate))
    priority = rng.random(len(geo.df))
    result = np.empty((4, 40, 20), dtype=np.int32)
    for fold in range(4):
        chosen = np.r_[rng.choice(geo.training_eligible[fold], 20, replace=False),
                       rng.choice(geo.test_eligible[fold], 20, replace=False)]
        for j, sid in enumerate(chosen):
            pool = geo.pools[(1, fold, int(sid))]
            result[fold, j] = pool[np.argsort(priority[pool], kind='stable')[:20]]
    return result


def labels_for(geo: support.Geometry, stage: str, world: str,
               amplitude: float, fraction: float, replicate: int) -> np.ndarray:
    # Return labels only: no oracle-response object is exposed to predictors.
    rng = np.random.default_rng(seed_for('labels', stage, world, amplitude, fraction, replicate))
    ng = base.unit_vectors(rng, (369, 3))
    ne = base.unit_vectors(rng, (369, 2))
    if world == 'environmental_shared':
        chosen = rng.permutation(369)[:round(369*fraction)]
        ne[chosen] = base.unit_vectors(rng, (1, 2))[0]
    pg = np.einsum('nd,nd->n', geo.g, ng[geo.sid])
    pe = np.einsum('nd,nd->n', geo.e, ne[geo.sid])
    if world.startswith('environmental'):
        projection = pe
    elif world == 'geographic_specific':
        projection = pg
    elif world == 'mixed_specific':
        projection = np.where((rng.random(369)<.5)[geo.sid], pg, pe)
    elif world == 'no_structure':
        projection = np.zeros(len(geo.df))
    else:
        raise ValueError('unknown synthetic world: '+world)
    signs = rng.choice([-1., 1.], 369)
    return amplitude*signs[geo.sid]*np.tanh(projection/.5) + rng.normal(size=len(geo.df)) > 0


def likelihood_features(x: np.ndarray, link: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Colour-blind log-probabilities; x shape (fold,species,photo,2)."""
    if x.ndim != 4 or x.shape[-1] != 2 or not np.isfinite(x).all():
        raise ValueError('finite (fold,species,photo,2) coordinates required')
    if link not in SLOPES:
        raise ValueError('unknown link')
    projection = np.einsum('fsnd,kd->fskn', x, AXES)
    if link == 'saturating_probit':
        latent = np.tanh(projection[..., None, :]/.5)*SLOPES[link][None, None, None, :, None]
        lp, lq = log_ndtr(latent), log_ndtr(-latent)
    else:
        latent = projection[..., None, :]*SLOPES[link][None, None, None, :, None]
        lp, lq = -np.logaddexp(0., -latent), -np.logaddexp(0., latent)
    return lp-lq, lq.sum(-1), (lp+lq).sum(-1)


def axis_log_likelihood(y: np.ndarray, features) -> np.ndarray:
    """Integrate ONE sign and strength for the entire species sample.

    y=(world,fold,species,photo), result=(world,fold,species,axis).
    """
    odds, zero, both = features
    if y.ndim != 4 or not np.isin(y, [0, 1]).all():
        raise ValueError('binary four-dimensional labels required')
    if y.shape[1:3] != odds.shape[:2] or y.shape[-1] != odds.shape[-1]:
        raise ValueError('label/feature shape mismatch')
    ll = np.einsum('wfsn,fskan->wfska', y.astype(float), odds, optimize=True)+zero
    signs = np.logaddexp(ll, both-ll)-np.log(2.)
    return logsumexp(signs, axis=-1)-np.log(odds.shape[-2])


def normalized_log_evidence(ll: np.ndarray) -> np.ndarray:
    if ll.ndim < 2 or ll.shape[-1] != K or not np.isfinite(ll).all():
        raise ValueError('finite axis likelihoods required')
    return ll-logsumexp(ll, axis=-1, keepdims=True)+np.log(K)


def train_posterior(train_lr: np.ndarray) -> np.ndarray:
    """Only training likelihood ratios are accepted; result (...,fraction,axis)."""
    if train_lr.shape[-1] != K or train_lr.shape[-2] == 0 or not np.isfinite(train_lr).all():
        raise ValueError('nonempty finite training likelihood ratios required')
    lf = np.full(5, -np.inf); lnf = np.full(5, -np.inf)
    np.log(FRACTIONS, out=lf, where=FRACTIONS>0)
    np.log(1-FRACTIONS, out=lnf, where=FRACTIONS<1)
    mix = np.logaddexp(lnf[:, None], lf[:, None]+train_lr[..., :, None, :])
    lp = mix.sum(axis=-3)
    return np.exp(lp-logsumexp(lp, axis=(-2,-1), keepdims=True))


def predictive_gain(posterior: np.ndarray, target_lr: np.ndarray) -> np.ndarray:
    """Per-target predictive gain; never updates the training posterior."""
    weighted_axis = np.einsum('...fk,f->...k', posterior, FRACTIONS)
    independent = np.einsum('...fk,f->...', posterior, 1-FRACTIONS)
    ratio = independent[..., None]+np.einsum('...k,...tk->...t', weighted_axis, np.exp(target_lr))
    return np.log(np.maximum(ratio, np.finfo(float).tiny))


def marginal_scores(lr: np.ndarray, n_train: int = 20, n_photos: int = 20):
    if lr.ndim != 4 or not 0<n_train<lr.shape[-2] or n_photos<=0:
        raise ValueError('invalid sharing score dimensions')
    tr, te = lr[..., :n_train, :], lr[..., n_train:, :]
    posterior = train_posterior(tr)
    pooled = predictive_gain(posterior, te).mean(-1)/n_photos
    separate = np.exp(tr).mean(-2)/K
    single = np.log(np.maximum(np.einsum('...k,...tk->...t', separate, np.exp(te)), np.finfo(float).tiny)).mean(-1)/n_photos
    mean_f = np.einsum('...fk,f->...', posterior, FRACTIONS)
    return {'joint_train': pooled, 'single_source_average': single}, mean_f


def qda_scores(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    n_worlds, folds, ns, n = y.shape
    xb = np.broadcast_to(x, (n_worlds,)+x.shape).reshape(-1,ns,n,2)
    yb = y.reshape(-1,ns,n)
    mu,cov,valid = base.fit_classes(xb,yb)
    pred = previous.gaussian_log_odds(xb[:,20:],mu[:,:20],cov[:,:20])>0
    ari = base.adjusted_rand_binary(yb[:,None,20:],pred)
    pairs = valid[:,:20,None]&valid[:,None,20:]
    return np.where(pairs,ari,0.).mean(axis=(1,2)).reshape(n_worlds,folds)


def summarize(output: Path) -> dict:
    df = pd.read_csv(output/'scores.csv')
    key=['stage','world','amplitude','shared_fraction','replicate','model']
    if len(df)!=30000 or df.duplicated(key).any():
        raise RuntimeError('score census failure')
    if not np.isfinite(df.select_dtypes(include='number')).all().all():
        raise RuntimeError('nonfinite scores')
    fold_cols=[f'score_fold{f}' for f in range(4)]
    if np.max(np.abs(df[fold_cols].mean(axis=1)-df.score))>1e-12:
        raise RuntimeError('fold averaging mismatch')
    quantiles=[]; thresholds=[]; rates=[]; paired=[]
    for model in MODELS:
        cuts=[]
        for w,a,f in NUISANCE:
            q=df[(df.stage=='calibration')&(df.model==model)&(df.world==w)&(df.amplitude==a)&(df.shared_fraction==f)]
            if len(q)!=REPS: raise RuntimeError('calibration arm census')
            cut=float(np.quantile(q.score,.975,method='higher'));cuts.append(cut)
            quantiles.append(dict(model=model,world=w,amplitude=a,quantile_975=cut))
        threshold=max(cuts)
        thresholds.append(dict(model=model,threshold=threshold))
        for w,a,f in NUISANCE+POSITIVE:
            q=df[(df.stage=='evaluation')&(df.model==model)&(df.world==w)&(df.amplitude==a)&(df.shared_fraction==f)]
            if len(q)!=REPS: raise RuntimeError('evaluation arm census')
            count=int((q.score>threshold).sum());lo,hi=base.wilson(count,REPS)
            rates.append(dict(model=model=model,world=w,amplitude=a,shared_fraction=f,count=count,
                              replicates=REPS,rate=count/REPS,wilson95_low=lo,wilson95_high=hi,
                              threshold=threshold,mean_score=float(q.score.mean()),
                              mean_train_fraction=float(q.train_fraction.mean())))
    cutmap={r['model']:r['threshold'] for r in thresholds}
    for w,a,f in POSITIVE:
        q=df[(df.stage=='evaluation')&(df.world==w)&(df.amplitude==a)&(df.shared_fraction==f)]
        wide=q.pivot(index='replicate',columns='model',values='score')
        comparisons=[(m,'qda_equal') for m in MODELS[:-1]]+[(link+'_joint_train',link+'_single_source_average') for link in SLOPES]
        for m,ref in comparisons:
            pred=(wide[m]>cutmap[m]).to_numpy();before=(wide[ref]>cutmap[ref]).to_numpy()
            diff=pred.astype(float)-before.astype(float)
            paired.append(dict(world=w,amplitude=a,shared_fraction=f,model=m,reference=ref,
                               gained=int((pred&~before).sum()),lost=int((~pred&before).sum()),
                               paired_difference=float(diff.mean()),mc_se=float(diff.std(ddof=1)/np.sqrt(REPS))))
    for name,records in [('calibration_quantiles',quantiles),('thresholds',thresholds),('decision_rates',rates),('paired_detections',paired)]:
        pd.DataFrame(records).to_csv(output/(name+'.csv'),index=False,lineterminator='\n')
    summary={'protocol':'hypervolume-pooled-response-evidence-v1','status':'complete_synthetic_run_pending_independent_verification',
             'specification_commit':SPEC,'seed':SEED,'unique_worlds':6000,'score_records':len(df),
             'models':list(MODELS),'full_sharing':[r for r in rates if r['world']=='environmental_shared' and r['shared_fraction']==1],
             'biological_colour_values_read':False,'real_climate_used':False,'true_parameters_used_in_predictors':False,
             'parent_decisions_modified':False,'empirical_inference_opened':False}
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary


def run(geometry: Path, output: Path):
    if output.exists(): raise FileExistsError('Use a fresh output; no scientific overwrite.')
    for name,expected in PARENT_HASHES.items():
        if digest(Path(__file__).with_name(name))!=expected:
            raise RuntimeError('parent runner drift: '+name)
    geo=support.Geometry(geometry)
    output.mkdir(parents=True)
    records=[];posterior_records=[];schedules={}
    for stage in ('calibration','evaluation'):
        scenarios=NUISANCE+(POSITIVE if stage=='evaluation' else ())
        schedules[stage]=np.stack([schedule(geo,stage,r) for r in range(REPS)])
        for r,idx in enumerate(schedules[stage]):
            x=geo.e[idx]
            y=np.stack([labels_for(geo,stage,w,a,f,r)[idx] for w,a,f in scenarios])
            scores={};train_fraction={}
            for link in SLOPES:
                ll=axis_log_likelihood(y,likelihood_features(x,link))
                scored,pf=marginal_scores(normalized_log_evidence(ll))
                for strategy,value in scored.items():
                    scores[link+'_'+strategy]=value
                    train_fraction[link+'_'+strategy]=pf if strategy=='joint_train' else np.full_like(pf,-1.)
            scores['qda_equal']=qda_scores(x,y)
            train_fraction['qda_equal']=np.full((len(scenarios),4),-1.)
            for ai,(w,a,f) in enumerate(scenarios):
                for model in MODELS:
                    values=scores[model][ai]
                    row=dict(stage=stage,world=w,amplitude=a,shared_fraction=f,replicate=r,model=model,
                             score=float(values.mean()),train_fraction=float(train_fraction[model][ai].mean()))
                    row.update({f'score_fold{k}':float(values[k]) for k in range(4)})
                    records.append(row)
            if (r+1)%25==0:
                print(json.dumps(dict(stage=stage,replicates=r+1,total=REPS,score_rows=len(records))),flush=True)
    np.savez_compressed(output/'schedules.npz',**schedules)
    pd.DataFrame(records).to_csv(output/'scores.csv',index=False,lineterminator='\n')
    pd.DataFrame(geo.capacity).to_csv(output/'fold_capacity.csv',index=False,lineterminator='\n')
    summary=summarize(output)
    provenance={'specification_commit':SPEC,'geometry_sha256':digest(geometry),'runner_sha256':digest(Path(__file__)),
                'parent_hashes':PARENT_HASHES,'python':platform.python_version(),'numpy':np.__version__,
                'rows':len(geo.df),'retained_rows':int(geo.mask.sum()),'eligible_species':len(geo.species),
                'species_role_overlap':len(geo.train_sid&geo.test_sid),'all_prior_decisions_unchanged':True}
    (output/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    checksums={p.name:digest(p) for p in output.iterdir() if p.is_file()}
    (output/'checksums.json').write_text(json.dumps(checksums,indent=2)+'\n')
    print(json.dumps({'status':summary['status'],'worlds':6000,'records':len(records)}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--geometry',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    run(args.geometry,args.output)
