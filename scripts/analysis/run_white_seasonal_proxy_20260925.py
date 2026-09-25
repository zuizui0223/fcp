#!/usr/bin/env python3
"""Post-hoc seasonal proxy diagnostic for the white-transition mechanism line.

This diagnostic is exploratory. It uses no external climate data. A positive
summer-proximity association is evidence for temporal sorting, not proof that
temperature causes white flowers.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import binomtest

STATES=["white","yellow_orange","red_pink","blue_purple"]
CONTRASTS={
    "white_vs_anthocyanic":["red_pink","blue_purple"],
    "white_vs_all_nonwhite":["yellow_orange","red_pink","blue_purple"],
    "white_vs_yellow_orange":["yellow_orange"],
}

def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for x in iter(lambda:f.read(1024*1024),b""): h.update(x)
    return h.hexdigest()

def bool_series(s):
    if s.dtype==bool: return s.fillna(False)
    return s.astype(str).str.strip().str.lower().isin({"true","1","yes"})

def holm(ps):
    out=[np.nan]*len(ps)
    valid=sorted([(i,float(p)) for i,p in enumerate(ps) if np.isfinite(p)],key=lambda x:x[1])
    m=len(valid); prev=0.0
    for rank,(i,p) in enumerate(valid):
        v=max(prev,min(1.0,(m-rank)*p)); out[i]=v; prev=v
    return out

def high_clip_photo_ids(high_path,join_path):
    high=pd.read_csv(high_path,dtype={"measurement_id":str})
    join=pd.read_csv(join_path,dtype={"measurement_id":str})
    ids=set(high.measurement_id.astype(str))
    out=pd.to_numeric(join.loc[join.measurement_id.astype(str).isin(ids),"photo_id"],errors="coerce").dropna().astype("int64")
    if len(out)!=len(ids): raise SystemExit("incomplete high-clip mapping")
    return set(out.tolist())

def screen(frame,metrics,comparison,min_per_state,nboot,rng):
    d=frame.loc[frame.morph.eq("white")|frame.morph.isin(comparison)].copy()
    d["is_white"]=d.morph.eq("white")
    rows=[]
    for metric in metrics:
        effects=[]
        for species,g in d.groupby("species",sort=True):
            x=pd.to_numeric(g[metric],errors="coerce")
            ok=x.notna(); g=g.loc[ok].copy(); x=x.loc[ok].astype(float)
            nw=int(g.is_white.sum()); nc=int((~g.is_white).sum())
            if nw<min_per_state or nc<min_per_state: continue
            sd=float(x.std(ddof=0))
            if not np.isfinite(sd) or sd<=0: continue
            z=(x-float(x.mean()))/sd
            effects.append(float(z.loc[g.is_white].mean()-z.loc[~g.is_white].mean()))
        a=np.asarray(effects,float)
        boots=np.array([rng.choice(a,size=len(a),replace=True).mean() for _ in range(nboot)])
        pos=int((a>0).sum()); neg=int((a<0).sum())
        rows.append({
            "metric":metric,"n_species":len(a),"mean_effect_sd":float(a.mean()),
            "median_effect_sd":float(np.median(a)),
            "bootstrap_ci_low":float(np.quantile(boots,.025)),
            "bootstrap_ci_high":float(np.quantile(boots,.975)),
            "positive_species":pos,"negative_species":neg,
            "sign_p":float(binomtest(pos,pos+neg,.5).pvalue),
        })
    out=pd.DataFrame(rows); out["sign_p_holm"]=holm(out.sign_p.tolist()); return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--measured",required=True)
    p.add_argument("--high-clip-ids",required=True)
    p.add_argument("--sealed-join-key",required=True)
    p.add_argument("--outdir",required=True)
    p.add_argument("--min-per-state",type=int,default=5)
    p.add_argument("--bootstraps",type=int,default=9999)
    p.add_argument("--seed",type=int,default=20260925)
    a=p.parse_args()
    measured=Path(a.measured); high=Path(a.high_clip_ids); join=Path(a.sealed_join_key)
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    df=pd.read_csv(measured,low_memory=False)
    req={"species","photo_id","latitude","observed_on","morph","global_classifiable"}
    if not req<=set(df): raise SystemExit(f"missing columns {sorted(req-set(df))}")
    keep=bool_series(df.global_classifiable)&df.morph.astype(str).isin(STATES)
    d=df.loc[keep].copy(); d["photo_id"]=pd.to_numeric(d.photo_id,errors="raise").astype("int64")
    hp=high_clip_photo_ids(high,join); d=d.loc[~d.photo_id.isin(hp)].copy()
    d["date"]=pd.to_datetime(d.observed_on,errors="coerce")
    d["doy"]=d.date.dt.dayofyear
    lat=pd.to_numeric(d.latitude,errors="coerce")
    d["abs_lat"]=lat.abs()
    peak=np.where(lat>=0,196,15)
    d["summer_index"]=np.cos(2*np.pi*(d.doy.to_numpy(float)-peak)/365.25)
    d["outside_tropics"]=lat.abs()>=23.5
    d["hemisphere"]=np.where(lat>=0,"N","S")

    rng=np.random.default_rng(a.seed)
    parts=[]
    for name,comparison in CONTRASTS.items():
        x=screen(d,["abs_lat","summer_index"],comparison,a.min_per_state,a.bootstraps,rng)
        x.insert(0,"contrast",name); parts.append(x)
    main=pd.concat(parts,ignore_index=True); main.to_csv(out/"proxy_summary.csv",index=False)

    sens=[]
    for name,comparison in CONTRASTS.items():
        x=screen(d.loc[d.outside_tropics],["summer_index"],comparison,a.min_per_state,a.bootstraps,rng)
        x.insert(0,"scope","outside_tropics"); x.insert(1,"contrast",name); sens.append(x)
        for hem in ["N","S"]:
            x=screen(d.loc[d.hemisphere.eq(hem)],["summer_index"],comparison,a.min_per_state,a.bootstraps,rng)
            x.insert(0,"scope",f"hemisphere_{hem}"); x.insert(1,"contrast",name); sens.append(x)
    sensitivity=pd.concat(sens,ignore_index=True); sensitivity.to_csv(out/"sensitivity_summary.csv",index=False)

    primary=main.loc[(main.contrast=="white_vs_anthocyanic")&(main.metric=="summer_index")].iloc[0]
    result={
        "schema":"white_seasonal_proxy_diagnostic_v1",
        "status":"complete","date_jst":"2026-09-25",
        "inference_status":"posthoc_exploratory_diagnostic_after_warm_side_partial_climate_signal_opened",
        "inputs":{"measured_sha256":sha256(measured),"high_clip_ids_sha256":sha256(high),"sealed_join_key_sha256":sha256(join)},
        "counts":{"post_highclip_classifiable_rows":int(len(d)),"species":int(d.species.nunique()),"high_clip_ids":int(len(hp))},
        "summer_index":"cosine proximity to day 196 in Northern Hemisphere or day 15 in Southern Hemisphere; +1 is local midsummer, -1 local midwinter",
        "primary":{
            "contrast":"white_vs_anthocyanic","n_species":int(primary.n_species),
            "mean_effect_sd":float(primary.mean_effect_sd),
            "bootstrap_ci":[float(primary.bootstrap_ci_low),float(primary.bootstrap_ci_high)],
            "sign_p":float(primary.sign_p),"sign_p_holm_two_proxies":float(primary.sign_p_holm)
        },
        "hard_nonclaims":[
            "calendar summer proximity is not measured temperature",
            "association can reflect morph phenology or correlated geography rather than thermal causation",
            "does not establish pigmented-to-white evolutionary direction",
            "does not test pollinator causation",
            "does not modify the frozen New Phytologist manuscript"
        ]
    }
    (out/"result.json").write_text(json.dumps(result,indent=2)+"\n")
    print(main.to_string(index=False)); print(sensitivity.to_string(index=False)); print(json.dumps(result,indent=2))

if __name__=="__main__": main()
