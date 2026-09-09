#!/usr/bin/env python3
"""Outcome-blind support audit for the real CHELSA coordinate export.

No flower-colour values are accepted. This audit only determines whether the
already-fixed species/region holdout design remains feasible after block-level
environmental missingness and describes the geometry available for later
synthetic method qualification.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import chi2

EXPORT_SHA = "33ab7c56df9bcd5a1f88ae97346ab7108d4a12cdffb4168359448b308ef5a4eb"
SEED = 20260907060
BUFFER_KM = 500.0
RADIUS_KM = 6371.0088
BLOCKS = {
    "thermal_regime": ["env_bio01", "env_bio04", "env_bio05", "env_bio06"],
    "water_balance": ["env_bio12", "env_bio14", "env_bio15", "env_cmi_mean"],
    "atmospheric_energy_dryness": ["env_vpd_mean", "env_rsds_mean", "env_gdd5", "env_sfcWind_mean"],
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def xyz_of(lat_deg: np.ndarray, lon_deg: np.ndarray) -> np.ndarray:
    lat, lon = np.deg2rad(lat_deg), np.deg2rad(lon_deg)
    return np.stack((np.cos(lat)*np.cos(lon), np.cos(lat)*np.sin(lon), np.sin(lat)), axis=1)


def chord_to_km(chord: np.ndarray) -> np.ndarray:
    return 2*RADIUS_KM*np.arcsin(np.clip(chord/2, 0, 1))


def shrink_covariance(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mu=x.mean(axis=0); dx=x-mu
    cov=dx.T@dx/max(len(x)-1,1)
    iso=np.trace(cov)/x.shape[1]
    cov=.8*cov+.2*iso*np.eye(x.shape[1])+1e-8*np.eye(x.shape[1])
    return mu,cov


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",type=Path,required=True)
    ap.add_argument("--output-json",type=Path,required=True)
    ap.add_argument("--output-capacity",type=Path,required=True)
    ap.add_argument("--output-support",type=Path,required=True)
    args=ap.parse_args()
    if sha(args.input)!=EXPORT_SHA: raise RuntimeError("real-climate export checksum drift")
    d=pd.read_csv(args.input)
    required=["photo_id","species","latitude","longitude","global_classifiable",*[x for v in BLOCKS.values() for x in v]]
    if list(d.columns)!=required or len(d)!=50000 or d.photo_id.duplicated().any(): raise RuntimeError("input shape/column drift")
    values=d.global_classifiable.astype(str).str.lower()
    if not values.isin(["true","false","1","0"]).all(): raise RuntimeError("bad classifiability mask")
    mask=values.isin(["true","1"]).to_numpy()
    counts=d.loc[mask].groupby("species").size()
    species=sorted(counts[counts>=40].index.astype(str))
    if len(species)!=369: raise RuntimeError("eligible species drift")
    q=d[d.species.astype(str).isin(species)].sort_values(["species","photo_id"],kind="stable").reset_index(drop=True)
    if len(q)!=36900: raise RuntimeError("eligible complete geometry drift")
    qmask=q.global_classifiable.astype(str).str.lower().isin(["true","1"]).to_numpy()
    if int(qmask.sum())!=21424: raise RuntimeError("eligible retained geometry drift")
    sid=pd.Categorical(q.species.astype(str),categories=species).codes.astype(np.int32)
    xyz=xyz_of(q.latitude.to_numpy(float),q.longitude.to_numpy(float))
    sector=np.floor(((q.longitude.to_numpy(float)+180)%360)/90).astype(int)
    order=sorted(range(369),key=lambda i:hashlib.sha256(("2026090702|species_split|"+species[i]).encode()).hexdigest())
    train_sid=set(order[:184]);test_sid=set(order[184:])

    capacity=[]; block_state={}; geometry_summary={}
    for block,cols in BLOCKS.items():
        x=q[cols].to_numpy(float); complete=np.isfinite(x).all(axis=1); retained=qmask&complete
        xr=x[retained]
        cov=np.cov(xr,rowvar=False,ddof=0); eig=np.linalg.eigvalsh(cov)
        eff=float(eig.sum()**2/np.square(eig).sum())
        corr=np.corrcoef(xr,rowvar=False)
        geometry_summary[block]={
            "variables":cols,"retained_complete_rows":int(retained.sum()),
            "retained_mean":xr.mean(axis=0).tolist(),"retained_sd":xr.std(axis=0).tolist(),
            "covariance_eigenvalues":eig.tolist(),"effective_rank_participation":eff,
            "max_abs_offdiag_correlation":float(np.max(np.abs(corr-np.eye(4)))),
            "correlation":corr.tolist(),
        }
        tr_el=[];te_el=[];pools={}
        for fold in range(4):
            inside=sector==fold
            dist=chord_to_km(cKDTree(xyz[inside]).query(xyz)[0])
            away=(~inside)&(dist>=BUFFER_KM)
            tr=[];te=[]
            for i,sp in enumerate(species):
                sm=sid==i
                trpool=np.flatnonzero(sm&retained&away).astype(np.int32)
                tepool=np.flatnonzero(sm&retained&inside).astype(np.int32)
                pools[(fold,i,"tr")]=trpool;pools[(fold,i,"te")]=tepool
                if i in train_sid and len(trpool)>=20:tr.append(i)
                if i in test_sid and len(tepool)>=20:te.append(i)
            if len(tr)<20 or len(te)<20: raise RuntimeError(f"not evaluable {block} fold {fold}: {len(tr)}/{len(te)}")
            tr_el.append(np.asarray(tr,dtype=int));te_el.append(np.asarray(te,dtype=int))
            capacity.append({"block":block,"fold":fold,"training_species":len(tr),"evaluation_species":len(te),
                             "retained_training_photos":int(np.count_nonzero(retained&away&np.isin(sid,tr))),
                             "retained_evaluation_photos":int(np.count_nonzero(retained&inside&np.isin(sid,te))),
                             "minimum_training_to_sector_photo_km":float(dist[retained&away&np.isin(sid,tr)].min())})
        block_state[block]=(x,tr_el,te_el,pools)

    support=[]
    for bi,(block,(x,tr_el,te_el,pools)) in enumerate(block_state.items()):
        rng=np.random.default_rng(SEED+bi)
        for replicate in range(100):
            priority=rng.random(len(q))
            for fold in range(4):
                tr=rng.choice(tr_el[fold],20,replace=False);te=rng.choice(te_el[fold],20,replace=False)
                sources=[];targets=[]
                for i in tr:
                    pool=pools[(fold,int(i),"tr")]
                    idx=pool[np.argsort(priority[pool],kind="stable")[:20]]
                    sources.append(shrink_covariance(x[idx]))
                for i in te:
                    pool=pools[(fold,int(i),"te")]
                    idx=pool[np.argsort(priority[pool],kind="stable")[:20]]
                    targets.append(x[idx])
                w=np.zeros((20,20),dtype=float)
                for a,(mu,cov) in enumerate(sources):
                    inv=np.linalg.inv(cov)
                    for t,xt in enumerate(targets):
                        delta=xt-mu
                        md=np.einsum("ni,ij,nj->n",delta,inv,delta,optimize=True)
                        w[a,t]=float(np.mean(md<=chi2.ppf(.95,4)))
                total=w.sum(axis=0);sq=np.square(w).sum(axis=0)
                ess=np.divide(total*total,sq,out=np.zeros_like(total),where=sq>0)
                support.append({"block":block,"replicate":replicate,"fold":fold,
                                "mean_source_target_photo_coverage":float(w.mean()),
                                "target_any_source_fraction":float(np.mean(total>0)),
                                "effective_source_count_mean":float(ess.mean()),
                                "effective_source_count_median":float(np.median(ess))})
    cap=pd.DataFrame(capacity);sup=pd.DataFrame(support)
    expected_tr=[133,135,141,165];expected_te=[60,62,59,26]
    for block in BLOCKS:
        z=cap[cap.block==block].sort_values("fold")
        if z.training_species.tolist()!=expected_tr or z.evaluation_species.tolist()!=expected_te:
            raise RuntimeError("block-specific species capacity differs from frozen design: "+block)
    summary=sup.groupby(["block","fold"],as_index=False).agg(
        mean_source_target_photo_coverage=("mean_source_target_photo_coverage","mean"),
        target_any_source_fraction=("target_any_source_fraction","mean"),
        effective_source_count_mean=("effective_source_count_mean","mean"),
    )
    payload={
        "protocol":"hypervolume-real-climate-support-audit-v1","status":"complete_outcome_blind_real_climate_support_audit",
        "seed":SEED,"input_sha256":EXPORT_SHA,"eligible_species":369,"retained_rows":21424,
        "training_species":184,"evaluation_species":185,"buffer_km":BUFFER_KM,
        "block_geometry":geometry_summary,"capacity":capacity,
        "support_summary":summary.to_dict(orient="records"),
        "support_schedules":100,"photos_per_species":20,"sources_per_fold":20,"targets_per_fold":20,
        "biological_colour_values_read":False,"synthetic_colour_labels_generated":False,
        "real_environment_inference_opened":False,"parent_empirical_decisions_modified":False,
        "interpretation":"Outcome-blind feasibility/support description only. Gaussian ellipsoid coverage is a support diagnostic, not ecological niche overlap or a biological result."
    }
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_capacity.parent.mkdir(parents=True,exist_ok=True)
    args.output_support.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(payload,indent=2)+"\n")
    cap.to_csv(args.output_capacity,index=False,lineterminator="\n")
    sup.to_csv(args.output_support,index=False,lineterminator="\n")
    print(json.dumps({"status":payload["status"],"effective_rank":{k:v["effective_rank_participation"] for k,v in geometry_summary.items()},"support_summary":payload["support_summary"]},indent=2))
    return 0

if __name__=="__main__": raise SystemExit(main())
