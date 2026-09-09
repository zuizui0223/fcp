#!/usr/bin/env python3
"""Frozen synthetic transfer benchmark on actual geometry; no real colours read.

Protocol committed before synthetic outcomes at db9e34175a22c8f21bde8f018db69b26dde7cc58.
Requires the coordinate-only export and its audit.json, plus the unchanged v1
Gaussian functions beside this file. All four folds contribute to every score.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import platform
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
import run_hypervolume_transfer_benchmark as base

SPECIFICATION_COMMIT = 'db9e34175a22c8f21bde8f018db69b26dde7cc58'
EXPORT_SHA = '5fa0c45a3f54f7441e2e372fa4ff4e813548c4cc44896a12a5ba8e97fb5af743'
SEED, RADIUS_KM, BUFFER_KM = 2026090702, 6371.0088, 500.0
N_TRAIN, N_TEST, N_PHOTOS, N_REPS = 20, 20, 20, 250
FRAMES = ('retained', 'complete_same_species')
MODES = ('species_only', 'species_and_region')
DESIGNS = tuple((f,m) for f in FRAMES for m in MODES)
NUISANCE = [('no_structure',0.,0.)] + [(w,a,0.) for w in ('geographic_specific','environmental_specific','mixed_specific') for a in (.5,1.,2.)]
POSITIVE = [(w,1.,f) for w in ('geographic_shared','environmental_shared') for f in (.25,.5,1.)]


def seed_for(*parts: object) -> int:
    return int.from_bytes(hashlib.sha256('|'.join(map(str,(SEED,)+parts)).encode()).digest()[:8],'little')


def sector_id(longitude: np.ndarray) -> np.ndarray:
    return np.floor(((np.asarray(longitude)+180.)%360.)/90.).astype(int)


def xyz_of(latitude: np.ndarray, longitude: np.ndarray) -> np.ndarray:
    lat,lon=np.deg2rad(latitude),np.deg2rad(longitude)
    return np.stack((np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)),axis=-1)


def chord_to_km(chord: np.ndarray) -> np.ndarray:
    return 2*RADIUS_KM*np.arcsin(np.clip(np.asarray(chord)/2,0,1))


class Geometry:
    """Only photo identity, species, coordinates and fixed measurement mask."""
    def __init__(self, path: Path):
        self.source_path=path
        meta=json.loads(path.with_name('audit.json').read_text())
        if hashlib.sha256(path.read_bytes()).hexdigest()!=EXPORT_SHA or meta['export_sha256']!=EXPORT_SHA or meta['source_sha256']!=base.SOURCE_SHA256:
            raise ValueError('coordinate export lineage mismatch')
        d=pd.read_csv(path)
        expected=['photo_id','species','latitude','longitude','global_classifiable']
        if list(d.columns)!=expected:
            raise ValueError('only the five frozen geometry/mask columns are accepted')
        if len(d)!=50000 or d.species.nunique()!=500 or d.photo_id.duplicated().any():
            raise ValueError('source counts/identities drifted')
        values=d.global_classifiable.astype(str).str.lower()
        if not values.isin(['true','false','1','0']).all():
            raise ValueError('unrecognized mask value')
        d['global_classifiable']=values.isin(['true','1'])
        if d.global_classifiable.sum()!=25377:
            raise ValueError('source mask count drifted')
        counts=d.loc[d.global_classifiable].groupby('species').size()
        species=sorted(counts[counts>=40].index)
        if len(species)!=369:
            raise ValueError('eligible species drifted')
        self.species=species
        d=d[d.species.isin(species)].sort_values(['species','photo_id'],kind='stable').reset_index(drop=True)
        if len(d)!=36900 or d.global_classifiable.sum()!=21424:
            raise ValueError('eligible geometry counts drifted')
        if not np.isfinite(d[['latitude','longitude']].to_numpy()).all() or not d.latitude.between(-90,90).all() or not d.longitude.between(-180,180).all():
            raise ValueError('invalid coordinates')
        self.df=d
        self.sid=pd.Categorical(d.species,categories=species).codes.astype(np.int32)
        self.xyz=xyz_of(d.latitude.to_numpy(),d.longitude.to_numpy())
        self.g=np.sqrt(3.)*self.xyz
        self.e=base.environment(self.xyz)
        self.mask=d.global_classifiable.to_numpy()
        self.sector=sector_id(d.longitude.to_numpy())
        order=sorted(range(369),key=lambda i:hashlib.sha256(('2026090702|species_split|'+species[i]).encode()).hexdigest())
        self.train_sid=set(order[:184]); self.test_sid=set(order[184:])
        self.pools={}; self.capacity=[]; self.species_capacity=[]
        self.training_eligible=[]; self.test_eligible=[]
        for fold in range(4):
            inside=self.sector==fold
            dist=chord_to_km(cKDTree(self.xyz[inside]).query(self.xyz)[0])
            away=(~inside)&(dist>=BUFFER_KM)
            tr,te=[],[]
            for i,sp in enumerate(species):
                species_mask=self.sid==i
                ntr=int(np.count_nonzero(species_mask&self.mask&away))
                nte=int(np.count_nonzero(species_mask&self.mask&inside))
                self.species_capacity.append({'fold':fold,'species':sp,'role':'training' if i in self.train_sid else 'evaluation','retained_outside_buffer':ntr,'retained_inside_sector':nte})
                if i in self.train_sid and ntr>=N_PHOTOS: tr.append(i)
                if i in self.test_sid and nte>=N_PHOTOS: te.append(i)
                for di,(frame,mode) in enumerate(DESIGNS):
                    fm=self.mask if frame=='retained' else np.ones(len(d),dtype=bool)
                    loc=(away if mode=='species_and_region' else np.ones(len(d),dtype=bool)) if i in self.train_sid else inside
                    self.pools[(di,fold,i)]=np.flatnonzero(species_mask&fm&loc).astype(np.int32)
            if len(tr)<N_TRAIN or len(te)<N_TEST:
                raise ValueError(f'not_evaluable fold {fold}: training={len(tr)},evaluation={len(te)}')
            self.training_eligible.append(np.asarray(tr));self.test_eligible.append(np.asarray(te))
            usable_away=away&np.isin(self.sid,tr)
            row={'fold':fold,'longitude_low':-180+90*fold,'longitude_high':-90+90*fold,'training_species':len(tr),'evaluation_species':len(te),'retained_training_photos':int(np.count_nonzero(usable_away&self.mask)),'retained_evaluation_photos':int(np.count_nonzero(inside&self.mask&np.isin(self.sid,te))),'minimum_training_to_sector_photo_km':float(dist[usable_away].min())}
            self.capacity.append(row)
        if [r['training_species'] for r in self.capacity]!=[133,135,141,165] or [r['evaluation_species'] for r in self.capacity]!=[60,62,59,26]:
            raise ValueError('metadata-only preflight and execution disagree')

    def schedule(self, stage: str, replicate: int) -> np.ndarray:
        """Same species, photo priority and evaluation targets across four designs."""
        rng=np.random.default_rng(seed_for('schedule',stage,replicate))
        priorities=rng.random(len(self.df))
        result=np.empty((4,4,N_TRAIN+N_TEST,N_PHOTOS),dtype=np.int32)
        for fold in range(4):
            chosen=np.r_[rng.choice(self.training_eligible[fold],N_TRAIN,replace=False),rng.choice(self.test_eligible[fold],N_TEST,replace=False)]
            for di in range(4):
                for j,i in enumerate(chosen):
                    pool=self.pools[(di,fold,int(i))]
                    result[di,fold,j]=pool[np.argsort(priorities[pool],kind='stable')[:N_PHOTOS]]
        return result


def labels_for(geo: Geometry, stage: str, world: str, amplitude: float, fraction: float, replicate: int) -> np.ndarray:
    """One coherent global label assignment per original photo per synthetic world."""
    rng=np.random.default_rng(seed_for('labels',stage,world,amplitude,fraction,replicate))
    ng=base.unit_vectors(rng,(369,3));ne=base.unit_vectors(rng,(369,2))
    if world.endswith('_shared'):
        chosen=rng.permutation(369)[:round(369*fraction)]
        if world=='geographic_shared':ng[chosen]=base.unit_vectors(rng,(1,3))[0]
        else:ne[chosen]=base.unit_vectors(rng,(1,2))[0]
    pg=np.einsum('nd,nd->n',geo.g,ng[geo.sid]);pe=np.einsum('nd,nd->n',geo.e,ne[geo.sid])
    if world.startswith('geographic'): projection=pg
    elif world.startswith('environmental'): projection=pe
    elif world=='mixed_specific': projection=np.where((rng.random(369)<.5)[geo.sid],pg,pe)
    elif world=='no_structure': projection=np.zeros(len(geo.df))
    else: raise ValueError(world)
    signs=rng.choice([-1.,1.],369)[geo.sid]
    return amplitude*signs*np.tanh(projection/.5)+rng.normal(size=len(geo.df))>0


def predict_qda(x: np.ndarray, mu: np.ndarray, cov: np.ndarray) -> np.ndarray:
    """Exact quadratic Gaussian log-odds evaluated as a batched feature product.

    Shapes: x=(B,target,photos,D), mu=(B,train,2,D), cov=(B,train,2,D,D).
    Algebraically equal to the parent direct Mahalanobis expression.
    """
    b,nt,np_,d=x.shape
    inv=np.linalg.inv(cov);logdet=np.linalg.slogdet(cov)[1]
    invmu=np.einsum('btkij,btkj->btki',inv,mu)
    quad=-.5*(inv[:,:,1]-inv[:,:,0]).reshape(b,-1,d*d)
    linear=invmu[:,:,1]-invmu[:,:,0]
    const=-.5*((mu[:,:,1]*invmu[:,:,1]).sum(-1)-(mu[:,:,0]*invmu[:,:,0]).sum(-1)+logdet[:,:,1]-logdet[:,:,0])
    coefficients=np.concatenate((quad,linear,const[...,None]),axis=-1)
    flat=x.reshape(b,nt*np_,d)
    xx=(flat[..., :,None]*flat[...,None,:]).reshape(b,nt*np_,d*d)
    features=np.concatenate((xx,flat,np.ones((b,nt*np_,1))),axis=-1)
    return (coefficients@features.swapaxes(-1,-2)>0).reshape(b,-1,nt,np_)


def score_domain(x: np.ndarray, labels: np.ndarray) -> dict[str,np.ndarray]:
    if x.shape[1:]!=(N_TRAIN+N_TEST,N_PHOTOS,x.shape[-1]) or x.shape[:-1]!=labels.shape:
        raise ValueError('unexpected score dimensions')
    mu,cov,valid=base.fit_classes(x,labels)
    predicted=predict_qda(x[:,N_TRAIN:],mu[:,:N_TRAIN],cov[:,:N_TRAIN])
    ari=base.adjusted_rand_binary(labels[:,None,N_TRAIN:],predicted)
    pairs=valid[:,:N_TRAIN,None]&valid[:,None,N_TRAIN:]
    nonconstant=(predicted.any(-1)&(~predicted).any(-1))
    return {'transfer':np.where(pairs,ari,0.).mean((1,2)),
            'separation':np.where(valid,1-base.affinity(mu,cov),0.).mean(1),
            'valid_train_fraction':valid[:,:N_TRAIN].mean(1),
            'valid_test_fraction':valid[:,N_TRAIN:].mean(1),
            'informative_pair_fraction':(pairs&nonconstant).mean((1,2))}


def write_rows(path: Path, rows: list[dict]) -> None:
    if not rows: raise ValueError('no rows')
    path.parent.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(rows).to_csv(path,index=False,lineterminator='\n')


def prepare(geo: Geometry, output: Path) -> dict[str,np.ndarray]:
    schedules={}
    for stage in ('calibration','evaluation'):
        schedules[stage]=np.stack([geo.schedule(stage,r) for r in range(N_REPS)])
    np.savez_compressed(output/'schedules.npz',**schedules)
    write_rows(output/'fold_capacity.csv',geo.capacity)
    write_rows(output/'species_capacity.csv',geo.species_capacity)
    geo.df.to_csv(output/'geometry_index.csv',index=False,lineterminator='\n')
    inclusions=[]
    for stage,s in schedules.items():
        for fold in range(4):
            species=geo.sid[s[:,0,fold,:,0]]
            for role,sub in [('training',species[:,:N_TRAIN]),('evaluation',species[:,N_TRAIN:])]:
                counts=np.bincount(sub.ravel(),minlength=369)
                inclusions.extend({'stage':stage,'fold':fold,'role':role,'species':geo.species[i],'inclusions':int(counts[i])} for i in range(369))
    write_rows(output/'species_inclusion_counts.csv',inclusions)
    photo_inclusions=[]
    for stage,s in schedules.items():
        for di,(frame,mode) in enumerate(DESIGNS):
            counts=np.bincount(s[:,di].ravel(),minlength=len(geo.df))
            for i in np.flatnonzero(counts):
                photo_inclusions.append({'stage':stage,'frame':frame,'mode':mode,'photo_id':geo.df.photo_id.iloc[i],'inclusions':int(counts[i])})
    write_rows(output/'photo_inclusion_counts.csv',photo_inclusions)
    audit={'protocol':'hypervolume-actual-support-blocked-transfer-v1','specification_commit':SPECIFICATION_COMMIT,'geometry_sha256':EXPORT_SHA,'species':369,'complete_rows':36900,'retained_rows':21424,'training_species':184,'evaluation_species':185,'species_overlap':len(geo.train_sid&geo.test_sid),'buffer_km':BUFFER_KM,'folds':geo.capacity,'biological_colour_values_read':False,'real_environment_used':False,'parent_decisions_modified':False,'selection_note':'The fixed classifiability mask is measurement-derived; not assumed random missingness.'}
    (output/'geometry_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    return schedules


def arm(geo: Geometry, schedules: dict[str,np.ndarray], output: Path, stage: str, world: str, amplitude: float, fraction: float, batch: int) -> Path:
    path=output/f'{stage}__{world}__a{amplitude:g}__f{fraction:g}.csv'
    if path.exists(): raise FileExistsError(path)
    rows=[]
    for start in range(0,N_REPS,batch):
        stop=min(start+batch,N_REPS); indices=schedules[stage][start:stop]
        labels=np.stack([labels_for(geo,stage,world,amplitude,fraction,r)[indices[r-start]] for r in range(start,stop)])
        y=labels.reshape(-1,N_TRAIN+N_TEST,N_PHOTOS)
        for representation,coordinate in [('geographic',geo.g),('environmental',geo.e)]:
            x=coordinate[indices].reshape(-1,N_TRAIN+N_TEST,N_PHOTOS,coordinate.shape[-1])
            scores={k:v.reshape(stop-start,4,4) for k,v in score_domain(x,y).items()}
            for r in range(start,stop):
                for di,(frame,mode) in enumerate(DESIGNS):
                    rec={'stage':stage,'world':world,'amplitude':amplitude,'shared_fraction':fraction,'replicate':r,'frame':frame,'mode':mode,'representation':representation}
                    for key,values in scores.items():
                        z=values[r-start,di];rec[key]=float(z.mean())
                        for f in range(4):rec[f'{key}_fold{f}']=float(z[f])
                    rows.append(rec)
    write_rows(path,rows)
    print(json.dumps({'stage':stage,'world':world,'amplitude':amplitude,'fraction':fraction,'worlds':N_REPS}),flush=True)
    return path


def summarize(output: Path) -> dict:
    paths=sorted(output.glob('*__*__a*__f*.csv'))
    if len(paths)!=26: raise ValueError(f'expected 26 arms, got {len(paths)}')
    df=pd.concat([pd.read_csv(p) for p in paths],ignore_index=True)
    key=['stage','world','amplitude','shared_fraction','replicate','frame','mode','representation']
    if len(df)!=52000 or df.duplicated(key).any() or not np.isfinite(df.select_dtypes(include='number').to_numpy()).all():
        raise ValueError('missing, duplicated or nonfinite score records')
    thresholds=[]; decisions=[]; flags=[]
    for frame,mode in DESIGNS:
        cut={}
        for domain in ('geographic','environmental'):
            q=df[(df.stage=='calibration')&(df.frame==frame)&(df['mode']==mode)&(df.representation==domain)]
            vals=[]
            for w,a,f in NUISANCE:
                z=q[(q.world==w)&(q.amplitude==a)].transfer.to_numpy()
                if len(z)!=N_REPS:raise ValueError('incomplete nuisance calibration')
                t=float(np.quantile(z,.975,method='higher'));vals.append(t)
                thresholds.append({'frame':frame,'mode':mode,'representation':domain,'nuisance':w,'amplitude':a,'threshold':t})
            cut[domain]=max(vals)
        for w,a,f in NUISANCE+POSITIVE:
            z=df[(df.stage=='evaluation')&(df.frame==frame)&(df['mode']==mode)&(df.world==w)&(df.amplitude==a)&(df.shared_fraction==f)]
            p=z.pivot(index='replicate',columns='representation',values='transfer')
            if len(p)!=N_REPS or p.isna().any().any(): raise ValueError('incomplete evaluation')
            g=p.geographic.to_numpy()>cut['geographic'];e=p.environmental.to_numpy()>cut['environmental']
            masks={'geography_only':g&~e,'environment_only':e&~g,'both':g&e,'neither':~g&~e,'any':g|e,'geographic_positive':g,'environmental_positive':e}
            for name,mask in masks.items():
                k=int(mask.sum());lo,hi=base.wilson(k,N_REPS)
                decisions.append({'frame':frame,'mode':mode,'world':w,'amplitude':a,'shared_fraction':f,'decision':name,'count':k,'replicates':N_REPS,'rate':k/N_REPS,'wilson95_low':lo,'wilson95_high':hi})
            for r in range(N_REPS):flags.append({'frame':frame,'mode':mode,'world':w,'amplitude':a,'shared_fraction':f,'replicate':r,'geographic_positive':bool(g[r]),'environmental_positive':bool(e[r])})
    write_rows(output/'calibration_thresholds.csv',thresholds);write_rows(output/'decision_rates.csv',decisions);write_rows(output/'evaluation_decisions.csv',flags)
    agg=['stage','frame','mode','world','amplitude','shared_fraction','representation']
    df.groupby(agg).mean(numeric_only=True).drop(columns='replicate').reset_index().to_csv(output/'mean_statistics.csv',index=False)
    contrasts=[]
    evaluation=df[df.stage=='evaluation']
    for (w,a,f,domain),group in evaluation.groupby(['world','amplitude','shared_fraction','representation']):
        p=group.pivot(index='replicate',columns=['frame','mode'],values='transfer')
        definitions=[('complete_minus_retained_'+m,('complete_same_species',m),('retained',m)) for m in MODES]+[('species_only_minus_blocked_'+fr,(fr,'species_only'),(fr,'species_and_region')) for fr in FRAMES]
        for label,plus,minus in definitions:
            delta=(p[plus]-p[minus]).to_numpy();se=float(delta.std(ddof=1)/np.sqrt(N_REPS));mean=float(delta.mean())
            contrasts.append({'world':w,'amplitude':a,'shared_fraction':f,'representation':domain,'contrast':label,'mean_paired_difference':mean,'paired_mc_se':se,'normal_mc95_low':mean-1.96*se,'normal_mc95_high':mean+1.96*se,'replicates':N_REPS})
    write_rows(output/'paired_score_contrasts.csv',contrasts)
    d=pd.DataFrame(decisions)
    nuisance=d[(~d.world.str.endswith('_shared'))&(d.decision=='any')]
    maxrow=nuisance.loc[nuisance.rate.idxmax()].to_dict()
    result={'protocol':'hypervolume-actual-support-blocked-transfer-v1','status':'complete_actual_support_synthetic_benchmark','specification_commit':SPECIFICATION_COMMIT,'seed':SEED,'unique_synthetic_worlds':6500,'designs':4,'folds':4,'fold_world_design_evaluations':104000,'score_records':52000,'calibration_replicates_per_arm':250,'evaluation_replicates_per_arm':250,'maximum_recorded_nuisance_any_positive':maxrow,'decisions':decisions,'biological_colour_values_read':False,'real_environment_used':False,'parent_decisions_modified':False,'runner_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'parent_runner_sha256':hashlib.sha256(Path(base.__file__).read_bytes()).hexdigest(),'geometry_sha256':EXPORT_SHA,'python':platform.python_version(),'numpy':np.__version__,'uncertainty':'Wilson intervals omit calibration-threshold uncertainty; rates conditional on frozen coordinates/mask and the finite synthetic nuisance family. No global error guarantee across designs.','claim_ceiling':'Synthetic response recoverability under actual support and species-and-region separation; not real climate, biological sharing prevalence, causal missingness or causal environment-vs-geography identification.'}
    (output/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    manifest={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in output.iterdir() if p.is_file() and p.name!='checksums.json'}
    (output/'checksums.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='decisions'},indent=2),flush=True)
    return result


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--geometry',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--batch-size',type=int,default=4)
    args=parser.parse_args()
    if not 1<=args.batch_size<=16: raise ValueError('batch size must be 1..16')
    if args.output_dir.exists(): raise FileExistsError('use a fresh output directory; no overwrites or favourable reruns')
    args.output_dir.mkdir(parents=True)
    geo=Geometry(args.geometry);schedules=prepare(geo,args.output_dir)
    print(json.dumps({'stage':'geometry_verified','folds':geo.capacity}),flush=True)
    for stage,worlds in [('calibration',NUISANCE),('evaluation',NUISANCE+POSITIVE)]:
        for w,a,f in worlds:arm(geo,schedules,args.output_dir,stage,w,a,f,args.batch_size)
    summarize(args.output_dir)
    return 0

if __name__=='__main__':raise SystemExit(main())
