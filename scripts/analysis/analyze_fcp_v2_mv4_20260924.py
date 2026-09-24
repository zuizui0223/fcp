#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "packages" / "disttrait" / "src"))

from disttrait.spatial import spatial_permutation_null

BIO = ["white","yellow","orange","red","pink","magenta","purple","blue","bronze"]
MORPHS = ["white","yellow_orange","red_pink","blue_purple"]
N_PERM = 999
SEED = 20260923
MIN_CLASSIFIABLE = 40
EARTH_KM = 6371.0088


def as_bool(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.strip().str.casefold().isin({"1","true","yes"})


def fraction_cols(prefix: str = "base") -> list[str]:
    return [f"{prefix}_fraction_{c}" for c in BIO]


def base_classifiable(df: pd.DataFrame) -> pd.Series:
    return (
        df["base_status"].astype(str).eq("classified_four_state_morph")
        & df["base_morph"].astype(str).isin(MORPHS)
    )


def maximum_span_km(lat: np.ndarray, lon: np.ndarray) -> float:
    lat=np.radians(np.asarray(lat,float)); lon=np.radians(np.asarray(lon,float))
    dlat=lat[:,None]-lat[None,:]; dlon=lon[:,None]-lon[None,:]
    a=np.sin(dlat/2.0)**2+np.cos(lat[:,None])*np.cos(lat[None,:])*np.sin(dlon/2.0)**2
    a=np.clip(a,0.0,1.0)
    return float(np.nanmax(2.0*EARTH_KM*np.arcsin(np.sqrt(a))))


def diversity(morph: pd.Series) -> float:
    c=morph.astype(str).value_counts().reindex(MORPHS,fill_value=0).to_numpy(float)
    p=c/c.sum()
    return float(1.0-np.sum(p*p))


def rank_residual(values: np.ndarray, controls: list[np.ndarray]) -> np.ndarray:
    y=rankdata(np.asarray(values,float),method="average")
    cols=[np.ones(len(y),float)]
    for v in controls:
        cols.append(rankdata(np.asarray(v,float),method="average"))
    X=np.column_stack(cols)
    b,*_=np.linalg.lstsq(X,y,rcond=None)
    return y-X@b


def partial_rank(x: np.ndarray,y: np.ndarray,controls: list[np.ndarray]) -> float:
    ex=rank_residual(x,controls); ey=rank_residual(y,controls)
    den=float(np.linalg.norm(ex)*np.linalg.norm(ey))
    return float(ex@ey/den) if den>1e-15 else float("nan")


def species_spatial(
    all_rows: pd.DataFrame,
    analysis_rows: pd.DataFrame,
    *,
    species: str,
    key_suffix: str,
) -> tuple[dict, np.ndarray] | None:
    classifiable=analysis_rows.loc[base_classifiable(analysis_rows)].copy()
    if len(classifiable)<MIN_CLASSIFIABLE:
        return None
    lat=pd.to_numeric(classifiable["latitude"],errors="coerce").to_numpy(float)
    lon=pd.to_numeric(classifiable["longitude"],errors="coerce").to_numpy(float)
    traits=classifiable[fraction_cols()].apply(pd.to_numeric,errors="coerce").to_numpy(float)
    finite=np.isfinite(lat)&np.isfinite(lon)&np.isfinite(traits).all(axis=1)&(traits.sum(axis=1)>0)
    lat=lat[finite]; lon=lon[finite]; traits=traits[finite]
    morph=classifiable.loc[finite,"base_morph"].astype(str)
    if len(lat)<MIN_CLASSIFIABLE:
        return None

    observed,null=spatial_permutation_null(
        lat,lon,traits,
        n_permutations=N_PERM,
        seed=SEED,
        key=f"FCPV2|{species}|{key_suffix}",
    )
    raw_lat=pd.to_numeric(all_rows["latitude"],errors="coerce").to_numpy(float)
    raw_lon=pd.to_numeric(all_rows["longitude"],errors="coerce").to_numpy(float)
    if not (np.isfinite(raw_lat).all() and np.isfinite(raw_lon).all()):
        raise RuntimeError(f"nonfinite frozen geometry for {species}")

    clear_tech_fail=all_rows["base_status"].astype(str).eq("not_evaluable_roi_or_flip_gate")
    row={
        "species":species,
        "panel":str(all_rows["panel"].iloc[0]),
        "n_classifiable":int(len(lat)),
        "D":diversity(morph),
        "spatial_observed_rho":float(observed),
        "log1p_span_all_100":float(np.log1p(maximum_span_km(raw_lat,raw_lon))),
        "technical_failure_rate_all_100":float(clear_tech_fail.mean()),
    }
    return row,null


def summarize(panel: pd.DataFrame,null: np.ndarray) -> dict:
    if len(panel)<3:
        return {"evaluable":False,"species":int(len(panel))}
    d=panel["D"].to_numpy(float)
    y=panel["spatial_observed_rho"].to_numpy(float)
    span=panel["log1p_span_all_100"].to_numpy(float)
    tech=panel["technical_failure_rate_all_100"].to_numpy(float)
    observed=partial_rank(d,y,[span,tech])
    vals=np.asarray(
        [partial_rank(d,null[:,j],[span,tech]) for j in range(null.shape[1])],
        float,
    )
    return {
        "evaluable":bool(np.isfinite(observed) and np.isfinite(vals).all()),
        "species":int(len(panel)),
        "partial_rho_D_spatial_adjusted_span_technical":float(observed),
        "p_upper_geometry_preserving_spatial_null":float((1+np.sum(vals>=observed-1e-15))/(len(vals)+1)),
        "null_mean":float(vals.mean()),
        "null_q025":float(np.quantile(vals,0.025)),
        "null_q975":float(np.quantile(vals,0.975)),
    }


def run_condition(df: pd.DataFrame, condition: str) -> dict:
    rows=[]; nulls=[]
    for species,g in df.groupby("species",sort=True):
        if len(g)!=100:
            raise RuntimeError(f"{species}: frozen 100-row denominator drift")
        if condition=="base":
            analysis=g
        elif condition=="roi_stable_only":
            analysis=g.loc[~as_bool(g["stratum_roi_unstable_fixed"])].copy()
        else:
            raise ValueError(condition)
        result=species_spatial(g,analysis,species=str(species),key_suffix=condition)
        if result is None:
            continue
        row,null=result
        rows.append(row); nulls.append(null)
    table=pd.DataFrame(rows)
    null=np.stack(nulls) if nulls else np.empty((0,N_PERM),float)
    overall=summarize(table,null) if len(table) else {"evaluable":False,"species":0}
    panels={}
    for label in ("P","N"):
        idx=np.flatnonzero(table["panel"].astype(str).to_numpy()==label) if len(table) else np.array([],int)
        panels[label]=summarize(table.iloc[idx].reset_index(drop=True),null[idx]) if len(idx) else {"evaluable":False,"species":0}
    return {
        "condition":condition,
        "eligible_species":int(len(table)),
        "panel_counts":{str(k):int(v) for k,v in table["panel"].value_counts().to_dict().items()} if len(table) else {},
        "pooled":overall,
        "panel_descriptive":panels,
        "species_table":table,
        "null":null,
    }


def main() -> None:
    p=argparse.ArgumentParser()
    p.add_argument("--biological-seal-dir",type=Path,required=True)
    p.add_argument("--output-dir",type=Path,required=True)
    args=p.parse_args()

    summary=json.loads((args.biological_seal_dir/"biological_summary.json").read_text())
    if summary.get("status")!="BIOLOGICAL_PASS_B_SEALED":
        raise RuntimeError("biological Pass B seal not terminal")
    if int(summary.get("rows",-1))!=40000 or int(summary.get("partition_receipts",-1))!=256:
        raise RuntimeError("Pass B census drift")
    if summary.get("source_identity",{}).get("replacement_rows")!=0:
        raise RuntimeError("replacement detected")

    df=pd.read_csv(args.biological_seal_dir/"biological_table.csv.gz",compression="gzip",dtype={"measurement_id":str})
    if len(df)!=40000 or df["species"].nunique()!=400:
        raise RuntimeError("biological table denominator drift")
    exact=df["source_identity_status"].astype(str).eq("exact_source_sha_match")
    work=df.loc[exact].copy()

    base=run_condition(work,"base")
    stable=run_condition(work,"roi_stable_only")

    out=args.output_dir
    out.mkdir(parents=True,exist_ok=True)
    base["species_table"].to_csv(out/"base_species_spatial.csv",index=False)
    stable["species_table"].to_csv(out/"roi_stable_species_spatial.csv",index=False)
    np.savez_compressed(out/"spatial_nulls.npz",base=base["null"],roi_stable=stable["null"])

    def clean(x: dict) -> dict:
        return {k:v for k,v in x.items() if k not in {"species_table","null"}}

    result={
        "schema":"fcp_v2_mv4_result_v1",
        "status":"FCP_V2_MV4_COMPLETE_WITH_FROZEN_CHANNEL_LIMITATION",
        "rows":40000,
        "source_exact_rows":int(exact.sum()),
        "per_species_spatial_null_replicates":N_PERM,
        "base":clean(base),
        "roi_stable_only":clean(stable),
        "matched_flower_minus_background":{
            "evaluable":False,
            "status":"NOT_EVALUABLE_FROZEN_TECHNICAL_CHANNEL_ABSENT",
            "reason":"The pre-biological v2 technical seal did not persist the background palette counts/fractions required to reproduce the existing matched flower-minus-background JSD-difference statistic. Background Lab summaries are available but are not substituted post hoc.",
            "rescue_attempted":False,
        },
        "claim_boundary":[
            "MV4 reuses the species-specific JSD spatial statistic and vertex-permutation null.",
            "The cross-species association is adjusted for all-100 sampled span and clear ROI/flip technical-failure rate.",
            "No ecological predictor is introduced.",
            "Matched-background exact control is reported not evaluable rather than reconstructed after biological opening.",
        ],
    }
    (out/"result.json").write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2,sort_keys=True,allow_nan=False))


if __name__=="__main__":
    main()
