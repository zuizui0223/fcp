#!/usr/bin/env python3
"""Support/response diagnostics on actual coordinates and synthetic colours only.

Specification committed BEFORE fresh outcomes: ee3f603f75f29e9fe60cce59d9d4d47552f8365e.
All alternatives are reported; none changes a frozen biological decision.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.stats import chi2
import run_hypervolume_actual_support_blocked as parent
import run_hypervolume_transfer_benchmark as base

SEED = 2026090703
REPS = 250
SPEC = 'ee3f603f75f29e9fe60cce59d9d4d47552f8365e'
MODELS = ('qda_equal', 'qda_empirical_prior', 'conditional_logit', 'oracle_response')
WEIGHTS = ('equal', 'overlap')
NUISANCE = parent.NUISANCE
POSITIVE = [('environmental_shared', 1., f) for f in (.25, .5, 1.)] + [('environmental_shared', 2., 1.)]
PARENT_HASHES = {
    'run_hypervolume_actual_support_blocked.py': 'dbedd48765816f6a2c8678585c62c52add32bdbd75a3fab59a66f06d8084891e',
    'run_hypervolume_transfer_benchmark.py': 'b5312f75128992f99b45898a6e5190aab1a86e18ad59573d26ca5fae66d8219b',
}


def digest(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def seed_for(*parts: object) -> int:
    return int.from_bytes(hashlib.sha256('|'.join(map(str, (SEED,)+parts)).encode()).digest()[:8], 'little')


def schedule(geo: parent.Geometry, stage: str, r: int) -> np.ndarray:
    rng = np.random.default_rng(seed_for('schedule', stage, r))
    priority = rng.random(len(geo.df))
    result = np.empty((4, 40, 20), dtype=np.int32)
    for fold in range(4):
        chosen = np.r_[rng.choice(geo.training_eligible[fold], 20, replace=False), rng.choice(geo.test_eligible[fold], 20, replace=False)]
        for j, sid in enumerate(chosen):
            pool = geo.pools[(1, fold, int(sid))]  # retained, buffered region holdout
            result[fold, j] = pool[np.argsort(priority[pool], kind='stable')[:20]]
    return result


def labels_and_truth(geo: parent.Geometry, stage: str, world: str, amplitude: float, fraction: float, r: int):
    rng = np.random.default_rng(seed_for('labels', stage, world, amplitude, fraction, r))
    ng = base.unit_vectors(rng, (369, 3)); ne = base.unit_vectors(rng, (369, 2))
    if world == 'environmental_shared':
        shared = rng.permutation(369)[:round(369*fraction)]
        ne[shared] = base.unit_vectors(rng, (1, 2))[0]
    elif world not in ('geographic_specific', 'environmental_specific', 'mixed_specific', 'no_structure'):
        raise ValueError(world)
    if world.startswith('environmental'): use_e = np.ones(369, dtype=bool)
    elif world == 'mixed_specific': use_e = rng.random(369) < .5
    else: use_e = np.zeros(369, dtype=bool)
    pg = np.einsum('nd,nd->n', geo.g, ng[geo.sid])
    pe = np.einsum('nd,nd->n', geo.e, ne[geo.sid])
    projection = np.where(use_e[geo.sid], pe, pg)
    if world == 'no_structure': projection = np.zeros(len(geo.df))
    signs = rng.choice([-1., 1.], 369)
    y = amplitude*signs[geo.sid]*np.tanh(projection/.5) + rng.normal(size=len(geo.df)) > 0
    return y, ng, ne, signs, use_e


def gaussian_log_odds(target: np.ndarray, mu: np.ndarray, cov: np.ndarray) -> np.ndarray:
    """Batched exact Gaussian density log ratio, with no prior term."""
    b, nt, n, d = target.shape
    inv = np.linalg.inv(cov); ld = np.linalg.slogdet(cov)[1]
    im = np.einsum('btkij,btkj->btki', inv, mu)
    quad = -.5*(inv[:, :, 1]-inv[:, :, 0]).reshape(b, -1, d*d)
    linear = im[:, :, 1]-im[:, :, 0]
    const = -.5*((mu[:, :, 1]*im[:, :, 1]).sum(-1)-(mu[:, :, 0]*im[:, :, 0]).sum(-1)+ld[:, :, 1]-ld[:, :, 0])
    coeff = np.concatenate((quad, linear, const[..., None]), axis=-1)
    x = target.reshape(b, nt*n, d)
    xx = (x[..., :, None]*x[..., None, :]).reshape(b, nt*n, d*d)
    features = np.concatenate((xx, x, np.ones((b, nt*n, 1))), axis=-1)
    return (coeff @ features.swapaxes(-1, -2)).reshape(b, -1, nt, n)


def overlap_weights(x: np.ndarray) -> np.ndarray:
    """Colour-blind fraction of target points in each source's 95% ellipsoid."""
    source, target = x[:, :20], x[:, 20:]
    mu = source.mean(-2); dx = source-mu[..., None, :]
    cov = base.shrink_covariance(np.einsum('btni,btnj->btij', dx, dx)/19)
    delta = target[:, None]-mu[:, :, None, None]
    distance = np.einsum('batni,baij,batnj->batn', delta, np.linalg.inv(cov), delta, optimize=True)
    return (distance <= chi2.ppf(.95, 2)).mean(-1)


def logistic_fit(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, float]:
    """Penalized conditional logit; hyperparameters fixed in the contract."""
    z = np.concatenate((np.ones(x.shape[:-1]+(1,)), x), axis=-1)
    beta = np.zeros(x.shape[:-2]+(3,)); penalty = np.array([1e-6, 1., 1.])
    gradient_max = float('inf')
    for _ in range(40):
        prob = expit(np.einsum('...ni,...i->...n', z, beta))
        grad = np.einsum('...ni,...n->...i', z, prob-y) + penalty*beta
        gradient_max = float(np.abs(grad).max())
        if gradient_max < 1e-8: break
        hess = np.einsum('...ni,...nj,...n->...ij', z, z, prob*(1-prob)) + np.diag(penalty)
        step = np.linalg.solve(hess, grad[..., None])[..., 0]
        # Fixed numerical safeguard, not a fitted model parameter.
        step /= np.maximum(1., np.max(np.abs(step), axis=-1, keepdims=True)/5.)
        beta -= step
    prob = expit(np.einsum('...ni,...i->...n', z, beta))
    grad = np.einsum('...ni,...n->...i', z, prob-y) + penalty*beta
    gradient_max = float(np.abs(grad).max())
    if not np.isfinite(beta).all() or gradient_max > 1e-6:
        raise RuntimeError(f'logit convergence failed: gradient={gradient_max}')
    return beta, gradient_max


def aggregate(ari: np.ndarray, valid: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Every target retains equal mass, including zero-support targets."""
    total = weights.sum(axis=1)
    numerator = (weights*np.where(valid, ari, 0.)).sum(axis=1)
    target_score = np.divide(numerator, total, out=np.zeros_like(numerator), where=total > 0)
    return target_score.mean(axis=1)


def score_batch(x: np.ndarray, g: np.ndarray, y: np.ndarray, truth, weights: np.ndarray):
    mu, cov, valid = base.fit_classes(x, y)
    odds = gaussian_log_odds(x[:, 20:], mu[:, :20], cov[:, :20])
    n1 = y[:, :20].sum(-1).astype(float)
    prior = np.log(np.maximum(n1, 1)/np.maximum(20-n1, 1))
    beta, grad = logistic_fit(x[:, :20], y[:, :20])
    logit = np.einsum('btd,bund->btun', beta[..., 1:], x[:, 20:]) + beta[..., 0, None, None]
    ng, ne, signs, use_e, structured = truth
    oracle_e = np.einsum('btd,bund->btun', ne, x[:, 20:])
    oracle_g = np.einsum('btd,bund->btun', ng, g[:, 20:])
    oracle = np.where(use_e[..., None, None], oracle_e, oracle_g)*signs[..., None, None]
    prediction = [odds > 0, odds+prior[..., None, None] > 0, logit > 0, (oracle > 0)&structured[:, None, None, None]]
    pairs = valid[:, :20, None] & valid[:, None, 20:]
    result = {}
    for name, pred in zip(MODELS, prediction):
        ari = base.adjusted_rand_binary(y[:, None, 20:], pred)
        info = pairs & pred.any(-1) & (~pred).any(-1)
        for weight_name, w in [('equal', np.ones_like(weights)), ('overlap', weights)]:
            result[(name, weight_name)] = {
                'transfer': aggregate(ari, pairs, w),
                'informative_weight_fraction': aggregate(info.astype(float), np.ones_like(pairs), w),
                'valid_pair_weight_fraction': aggregate(pairs.astype(float), np.ones_like(pairs), w),
            }
    return result, grad


def write_table(path: Path, rows) -> None:
    pd.DataFrame(rows).to_csv(path, index=False, lineterminator='\n')


def summarize(output: Path) -> dict:
    df = pd.read_csv(output/'scores.csv')
    key = ['stage','world','amplitude','shared_fraction','replicate','model','weight']
    if len(df) != 48000 or df.duplicated(key).any(): raise ValueError('score census failure')
    if not np.isfinite(df.select_dtypes(include='number').to_numpy()).all(): raise ValueError('nonfinite scores')
    quantiles, thresholds, rates = [], [], []
    for model in MODELS:
        for weight in WEIGHTS:
            all_cal = df[(df.stage=='calibration')&(df.model==model)&(df.weight==weight)]
            cuts = []
            for world, amp, frac in NUISANCE:
                q = all_cal[(all_cal.world==world)&(all_cal.amplitude==amp)]
                if len(q)!=REPS: raise ValueError('calibration counts')
                value = float(np.quantile(q.transfer, .975, method='higher'))
                quantiles.append(dict(model=model,weight=weight,world=world,amplitude=amp,threshold=value))
                cuts.append((amp, value))
            for scheme, known in [('pooled_unknown_strength', None), ('oracle_known_strength_1', 1.), ('oracle_known_strength_2', 2.)]:
                cut = max(value for amp,value in cuts if known is None or amp in (0.,known))
                thresholds.append(dict(model=model,weight=weight,calibration=scheme,threshold=cut))
                for world, amp, frac in NUISANCE+POSITIVE:
                    q = df[(df.stage=='evaluation')&(df.model==model)&(df.weight==weight)&(df.world==world)&(df.amplitude==amp)&(df.shared_fraction==frac)]
                    if len(q)!=REPS: raise ValueError('evaluation counts')
                    k = int((q.transfer>cut).sum()); lo,hi=base.wilson(k,REPS)
                    rates.append(dict(model=model,weight=weight,calibration=scheme,world=world,amplitude=amp,shared_fraction=frac,count=k,replicates=REPS,rate=k/REPS,wilson95_low=lo,wilson95_high=hi,threshold=cut,mean_transfer=float(q.transfer.mean()),informative_weight_fraction=float(q.informative_weight_fraction.mean())))
    write_table(output/'calibration_quantiles.csv',quantiles)
    write_table(output/'thresholds.csv',thresholds)
    write_table(output/'decision_rates.csv',rates)
    contrasts=[]
    for world,amp,frac in NUISANCE+POSITIVE:
        q=df[(df.stage=='evaluation')&(df.world==world)&(df.amplitude==amp)&(df.shared_fraction==frac)]
        wide=q.pivot(index='replicate',columns=['model','weight'],values='transfer')
        comparisons=[(m,'equal','qda_equal','equal') for m in MODELS[1:]]+[(m,'overlap',m,'equal') for m in MODELS]
        for m,w,ref,rw in comparisons:
            diff=(wide[(m,w)]-wide[(ref,rw)]).to_numpy(); se=float(diff.std(ddof=1)/np.sqrt(REPS)); mean=float(diff.mean())
            contrasts.append(dict(world=world,amplitude=amp,shared_fraction=frac,model=m,weight=w,reference_model=ref,reference_weight=rw,paired_mean=mean,mc_se=se,mc95_low=mean-1.96*se,mc95_high=mean+1.96*se))
    write_table(output/'paired_contrasts.csv',contrasts)
    summary={'protocol':'hypervolume-support-response-diagnostic-v1','status':'complete_synthetic_diagnostic','specification_commit':SPEC,'seed':SEED,'unique_worlds':6000,'score_records':len(df),'fold_score_records':len(df)*4,'real_climate_used':False,'biological_colour_values_read':False,'empirical_decisions_modified':False,'unknown_strength_full_sharing':[r for r in rates if r['calibration']=='pooled_unknown_strength' and r['world']=='environmental_shared' and r['shared_fraction']==1.], 'claim_ceiling':'Conditional synthetic recoverability only, not real environmental sharing or optimal-power bounds; all eight diagnostic designs retained.'}
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary


def run(geometry_path: Path, output: Path, batch: int) -> None:
    if output.exists(): raise FileExistsError('Use a fresh output directory; no favourable overwrite.')
    if batch<1 or batch>25: raise ValueError('batch must be 1..25')
    for name, expected in PARENT_HASHES.items():
        if digest(Path(__file__).with_name(name))!=expected: raise ValueError('parent code drift: '+name)
    geo=parent.Geometry(geometry_path)
    output.mkdir(parents=True)
    schedules={stage:np.stack([schedule(geo,stage,r) for r in range(REPS)]) for stage in ('calibration','evaluation')}
    np.savez_compressed(output/'schedules.npz',**schedules)
    overlap={}; support=[]
    for stage,idx in schedules.items():
        w=overlap_weights(geo.e[idx].reshape(-1,40,20,2)).reshape(REPS,4,20,20); overlap[stage]=w
        total=w.sum(2);sq=(w*w).sum(2)
        ess=np.divide(total*total,sq,out=np.zeros_like(total),where=sq>0)
        for r in range(REPS):
            for fold in range(4):
                support.append(dict(stage=stage,replicate=r,fold=fold,mean_pair_overlap=float(w[r,fold].mean()),nonzero_pair_fraction=float((w[r,fold]>0).mean()),target_coverage_fraction=float((total[r,fold]>0).mean()),effective_source_count_mean=float(ess[r,fold].mean()),effective_source_count_median=float(np.median(ess[r,fold]))))
    write_table(output/'support_geometry.csv',support)
    records=[]; max_gradient=0.
    for stage in ('calibration','evaluation'):
        for world,amp,frac in NUISANCE+(POSITIVE if stage=='evaluation' else []):
            for start in range(0,REPS,batch):
                stop=min(REPS,start+batch); idx=schedules[stage][start:stop]; b=stop-start
                assignments=[labels_and_truth(geo,stage,world,amp,frac,r) for r in range(start,stop)]
                y=np.stack([z[0][idx[j]] for j,z in enumerate(assignments)]).reshape(-1,40,20)
                train_ids=geo.sid[idx[:,:,:20,0]]
                truth=[np.stack([z[k][train_ids[j]] for j,z in enumerate(assignments)]).reshape((-1,20)+(assignments[0][k].shape[1:])) for k in (1,2,3,4)]
                truth.append(np.full(b*4,world!='no_structure',dtype=bool))
                x=geo.e[idx].reshape(-1,40,20,2); g=geo.g[idx].reshape(-1,40,20,3)
                scores,grad=score_batch(x,g,y,truth,overlap[stage][start:stop].reshape(-1,20,20))
                max_gradient=max(max_gradient,grad)
                for (model,weight),values in scores.items():
                    reshaped={k:v.reshape(b,4) for k,v in values.items()}
                    for j,r in enumerate(range(start,stop)):
                        rec=dict(stage=stage,world=world,amplitude=amp,shared_fraction=frac,replicate=r,model=model,weight=weight)
                        for k,v in reshaped.items():
                            rec[k]=float(v[j].mean())
                            for fold in range(4):rec[f'{k}_fold{fold}']=float(v[j,fold])
                        records.append(rec)
            print(json.dumps(dict(stage=stage,world=world,amplitude=amp,fraction=frac,completed=REPS)),flush=True)
    write_table(output/'scores.csv',records)
    result=summarize(output)
    result['logit_max_absolute_gradient']=max_gradient
    result['python_version']=platform.python_version();result['numpy_version']=np.__version__
    result['runner_sha256']=digest(Path(__file__))
    result['parent_runner_sha256']=PARENT_HASHES
    result['geometry_sha256']=digest(geometry_path)
    (output/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    (output/'checksums.json').write_text(json.dumps({p.name:digest(p) for p in sorted(output.iterdir()) if p.is_file()},indent=2)+'\n')
    print('COMPLETE: all 6000 synthetic worlds; no empirical decision changed.',flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--geometry',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--batch',type=int,default=10)
    args=ap.parse_args();run(args.geometry,args.output,args.batch)
