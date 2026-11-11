#!/usr/bin/env python3
"""FCP within-species geographic photo-colour sorting after sampling-time controls.

Input is the SHA-authenticated immutable 447,500 original measured photos.
A species' observed colour multiset, cells, observer/photo IDs, calendar
month/year and cell sample sizes are fixed. We shuffle complete COARSE photo
labels (a) globally, (b) within month-of-year, (c) within exact year-month.
These are strictly POST-HOC diagnostic nulls, not adaptation/IBD/IBE/fitness
tests. No exact point coordinates or observation identifiers leave this script.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from fcp_global_photo_atlas_20261020 import (
    STATES, check, geographic_units, group_stats, high_purity, read_original,
)
from fcp_equal_area_spatial_null_20261020 import score

SEED = "FCP_WITHIN_SPECIES_MONTH_NULL_20261020"
MODES = ("unrestricted", "calendar_month", "calendar_year_month")
N_REPS = 199
TYPES = ("tv", "white_abs", "nonwhite_residual")
COHORTS = ("discovery", "validation", "prospective", "wave4", "new_U100")


def date_filter(data):
    parsed=pd.to_datetime(data.observed_on,errors="coerce")
    keep=parsed.notna() & parsed.dt.year.between(2000,2026)
    out=data.loc[keep].copy()
    out["month"]=parsed.loc[keep].dt.month.to_numpy(dtype=np.int16)
    out["year_month"]=(parsed.loc[keep].dt.year*12+parsed.loc[keep].dt.month).to_numpy(dtype=np.int32)
    return out,{"classifiable":int(len(data)), "dated_2000_2026":int(len(out)),
                "missing_or_out_of_range":int(len(data)-len(out))}


def make_matched_cell_sets(d):
    full=geographic_units(d)
    all_units=group_stats(full,["inat_taxon_id","cell"],min_n=10,min_observers=3)
    pure=geographic_units(high_purity(d))
    q_units=group_stats(pure,["inat_taxon_id","cell"],min_n=10,min_observers=3)
    candidate=q_units.groupby("inat_taxon_id").cell.nunique()
    candidate=set(candidate.loc[candidate>=2].index)
    keep=q_units.loc[q_units.inat_taxon_id.isin(candidate),["inat_taxon_id","cell"]]
    orig=full.merge(keep,on=["inat_taxon_id","cell"],how="inner",validate="many_to_one")
    good=pure.merge(keep,on=["inat_taxon_id","cell"],how="inner",validate="many_to_one")
    sets={
        "same_taxa_cells_original_quality":orig,
        "same_taxa_cells_high_purity":good,
    }
    check(orig.inat_taxon_id.nunique()==good.inat_taxon_id.nunique() and
          set(orig.inat_taxon_id.unique())==set(good.inat_taxon_id.unique()),
          "quality comparison has different species")
    check(keep.drop_duplicates().shape[0]==len(keep),"nonunique taxon-cell inclusion")
    return sets,{
        "classifiable_in_valid_dated_geography":int(len(full)),
        "strict_purity_dated_geography":int(len(pure)),
        "taxa_matched_at_least_two_cells":int(len(candidate)),
        "shared_taxon_cell_units":int(len(keep)),
        "matched_original_photos":int(len(orig)),
        "matched_pure_photos":int(len(good)),
    }


def stable_seed(*items):
    payload="|".join(map(str,(SEED,*items))).encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8],"little")


def strata_groups(g,mode):
    if mode=="unrestricted":
        return [np.arange(len(g),dtype=int)]
    col="month" if mode=="calendar_month" else "year_month"
    return [np.asarray(ix,dtype=int) for _,ix in
            g.groupby(col,sort=True).indices.items()]


def opportunity_counts(g,blocks):
    cells=g.cell.to_numpy()
    observers=g.observer_id.to_numpy()
    obs_ids=g.observation_id.to_numpy()
    crosscell_positions=0
    n_movable_blocks=0
    independent_pairs=0
    for idx in blocks:
        if len(idx)<2 or len(set(cells[idx]))<2:
            continue
        n_movable_blocks+=1
        crosscell_positions+=len(idx)
        u,v=np.triu_indices(len(idx),1)
        i=idx[u];j=idx[v]
        different=(cells[i]!=cells[j])
        validobs=pd.notna(observers[i]) & pd.notna(observers[j])
        validid=pd.notna(obs_ids[i]) & pd.notna(obs_ids[j])
        independent_pairs+=int(np.count_nonzero(different & validobs & validid
               & (observers[i]!=observers[j]) & (obs_ids[i]!=obs_ids[j])))
    return {
        "n_cross_cell_exchangeable_blocks":n_movable_blocks,
        "n_cross_cell_exchangeable_positions":crosscell_positions,
        "n_cross_cell_distinct_observer_photo_pairs":independent_pairs,
    }


def taxon_test(g,tag,reps=N_REPS):
    g=g.sort_values(["cell","photo_id"],kind="stable").reset_index(drop=True)
    cells,cell_n=np.unique(g.cell.to_numpy(),return_counts=True)
    if len(cells)<2:return None
    check(g.fcp_source_cohort.nunique()==1,"source cohort overlap")
    categories=pd.Categorical(g.morph,categories=STATES).codes
    check((categories>=0).all(),"invalid morph label")
    n=len(g)
    start=np.insert(np.cumsum(cell_n)[:-1],0,0)
    den=cell_n[:,None]
    lookup=np.eye(len(STATES),dtype=np.int16)
    def metric(labels):
        cellcounts=np.add.reduceat(lookup[labels],start,axis=0)
        return np.asarray(score(cellcounts/den),dtype=float)
    observed=metric(categories)
    cohort=str(g.fcp_source_cohort.iloc[0])
    tid=int(g.inat_taxon_id.iloc[0])
    block_by={m:strata_groups(g,m) for m in MODES}
    info={m:opportunity_counts(g,block_by[m]) for m in MODES if m!="unrestricted"}
    # For season-constrained worlds, only strata spanning >1 geographical
    # cell can alter the statistic. Shuffling single-cell strata is redundant.
    arrs={}
    for mode in MODES:
        blocks=[q for q in block_by[mode] if len(q)>1 and
                (mode=="unrestricted" or len(set(g.cell.to_numpy()[q]))>1)]
        rng=np.random.default_rng(stable_seed(tag,tid,mode))
        null=np.empty((reps,3),dtype=float)
        for k in range(reps):
            x=categories.copy()
            for block in blocks:
                x[block]=rng.permutation(categories[block])
            null[k]=metric(x)
        arrs[mode]=null
    row={"taxon_id":tid,"cohort":cohort,"n_photos":n,"n_cells":len(cells),
         "n_observed_white":int(np.count_nonzero(categories==0)),
         "n_observed_nonwhite":int(np.count_nonzero(categories!=0)),
         "n_months_observed":int(g.month.nunique()),
         "n_year_months_observed":int(g.year_month.nunique()),
         **{"observed_"+t:float(v) for t,v in zip(TYPES,observed)}}
    for mode in MODES:
        for i,t in enumerate(TYPES):
            row[f"{mode}_null_{t}"]=float(arrs[mode][:,i].mean())
            row[f"{mode}_excess_{t}"]=float(observed[i]-arrs[mode][:,i].mean())
        if mode!="unrestricted":
            row.update({f"{mode}_{k}":v for k,v in info[mode].items()})
    return row,arrs


def self_test():
    # Two geographic photo-sampled cells. All colour changes are seasonal,
    # with 18 January white / 2 June red in cell A and reverse in cell B.
    # Unrestricted shuffle destroys the spurious spatial seasonal imbalance;
    # calendar-month shuffle leaves it exactly unchanged.
    d=[]
    for cell in ("EA03_01","EA03_02"):
        njan=18 if cell.endswith("01") else 2
        for i in range(20):
            jan=i<njan
            d.append({"inat_taxon_id":1,"fcp_source_cohort":"discovery","cell":cell,
                      "photo_id":f"{cell}:{i}","observer_id":f"o{cell}{i}",
                      "observation_id":f"obs{cell}{i}","month":1 if jan else 6,
                      "year_month":2026*12+(1 if jan else 6),
                      "morph":"white" if jan else "red_pink"})
    row,worlds=taxon_test(pd.DataFrame(d),"fixture",reps=99)
    check(row["calendar_month_n_cross_cell_exchangeable_blocks"]==2,
          "monthly cross-cell support should exist")
    check(abs(row["calendar_month_excess_tv"])<1e-12,
          "purely time-composition-driven geography falsely retained")
    check(abs(row["calendar_year_month_excess_tv"])<1e-12,
          "purely year-month-composition-driven geography falsely retained")
    check(row["unrestricted_excess_tv"]>0.3,
          "unrestricted null must see seasonal spatial sampling")
    # Exact score decomposition also applies after conditioning.
    for mode in MODES:
        check(abs(row[f"{mode}_excess_tv"]-row[f"{mode}_excess_white_abs"]-
                  row[f"{mode}_excess_nonwhite_residual"])<1e-11,
              "additive decomposition broken")
    print("FCP_MONTH_MATCHED_GLOBAL_GRID_SYNTHETIC_PASS")


def summarize(frame,per_species,mode,support_filter):
    if support_filter=="geometry_ge5_independent_pairs" and mode!="unrestricted":
        mask=(frame[f"{mode}_n_cross_cell_distinct_observer_photo_pairs"]>=5).to_numpy()
    elif support_filter=="geometry_ge5_independent_pairs" and mode=="unrestricted":
        # Align unrestricted null to the exact subset for each time mode,
        # only by invoking this function with its actual filtered rows.
        raise RuntimeError("invalid geometry comparison")
    else:
        mask=np.ones(len(frame),bool)
    df=frame.loc[mask]
    if len(df)==0:
        return {"n_species":0,"status":"not_estimable_no_exchangeable_geometry"}
    vals=np.stack([per_species[tid][mode] for tid in df.taxon_id],axis=0)
    observed=df[["observed_"+t for t in TYPES]].mean().to_numpy(float)
    exp=vals.mean(axis=0)
    nullmean=exp.mean(axis=0)
    p=(1+np.count_nonzero(exp>=observed[None,:],axis=0))/(vals.shape[1]+1)
    return {
        "n_species":int(len(df)),
        "cohort_counts":{str(k):int(v) for k,v in df.cohort.value_counts().items()},
        "mean_observed":dict(zip(TYPES,map(float,observed))),
        "mean_null":dict(zip(TYPES,map(float,nullmean))),
        "mean_excess":dict(zip(TYPES,map(float,observed-nullmean))),
        "one_sided_199_world_vertex_p":dict(zip(TYPES,map(float,p))),
        "fraction_taxa_month_geometry_5plus":float(np.mean(
          frame["calendar_month_n_cross_cell_distinct_observer_photo_pairs"]>=5)),
        "fraction_taxa_year_month_geometry_5plus":float(np.mean(
          frame["calendar_year_month_n_cross_cell_distinct_observer_photo_pairs"]>=5)),
        "posthoc_exploratory":True,
    }


def cohort_sensitivity(frame,per_species):
    """Same original nulls evaluated by disjoint taxon cohort; no per-cohort refit."""
    out={}
    for cohort in (*COHORTS,"leave_out_new_U100"):
        if cohort=="leave_out_new_U100":
            sub=frame.loc[frame.cohort.ne("new_U100")]
        else:
            sub=frame.loc[frame.cohort.eq(cohort)]
        if len(sub)<20:
            out[cohort]={"n_species":int(len(sub)),"status":"insufficient_species"}
            continue
        subres={}
        for mode in MODES:
            x=summarize(sub,per_species,mode,"all")
            subres[mode]={
                "n_species":x["n_species"],
                "mean_excess_tv":x["mean_excess"]["tv"],
                "p_upper_tv":x["one_sided_199_world_vertex_p"]["tv"],
                "mean_excess_white_photo_frequency":x["mean_excess"]["white_abs"],
                "p_upper_white":x["one_sided_199_world_vertex_p"]["white_abs"],
            }
        u=subres["unrestricted"]["mean_excess_tv"]
        m=subres["calendar_month"]["mean_excess_tv"]
        ym=subres["calendar_year_month"]["mean_excess_tv"]
        out[cohort]={
          "n_species":int(len(sub)),
          "source_cohort":cohort,
          "n_with_month_crosscell_observer_pair_5plus":int(
              (sub.calendar_month_n_cross_cell_distinct_observer_photo_pairs>=5).sum()),
          "n_with_year_month_crosscell_observer_pair_5plus":int(
              (sub.calendar_year_month_n_cross_cell_distinct_observer_photo_pairs>=5).sum()),
          "worlds":subres,
          "absolute_drop_excess_from_unrestricted_to_month":float(u-m),
          "absolute_drop_excess_from_month_to_yearmonth":float(m-ym),
          "yearmonth_to_unrestricted_excess_ratio":float(ym/u) if abs(u)>1e-12 else None,
          "status":"posthoc_source_disjoint_photo_cohort_comparison",
        }
    return out


def analyze_one(d,tag):
    per_species={}
    rows=[]
    for tid,g in d.groupby("inat_taxon_id",sort=True):
        test=taxon_test(g,tag)
        if test is None:continue
        row,arr=test
        rows.append(row)
        per_species[tid]=arr
    frame=pd.DataFrame(rows)
    check(len(frame)>0,"no geographically eligible species after date matching")
    results={
        "n_species":int(len(frame)),
        "source_cohort_sensitivity":cohort_sensitivity(frame,per_species),
        "all_eligible":{m:summarize(frame,per_species,m,"all")
                        for m in MODES},
        "month_geometry_5plus_independent_cross_cell_pairs":{
            m:summarize(frame.loc[frame.calendar_month_n_cross_cell_distinct_observer_photo_pairs>=5],
                        per_species,m,"all") for m in MODES},
        "year_month_geometry_5plus_independent_cross_cell_pairs":{
            m:summarize(frame.loc[frame.calendar_year_month_n_cross_cell_distinct_observer_photo_pairs>=5],
                        per_species,m,"all") for m in MODES},
        "n_taxa_with_month_geometry_5plus":int(
           (frame.calendar_month_n_cross_cell_distinct_observer_photo_pairs>=5).sum()),
        "n_taxa_with_year_month_geometry_5plus":int(
           (frame.calendar_year_month_n_cross_cell_distinct_observer_photo_pairs>=5).sum()),
        "n_taxa_with_any_month_crosscell_exchangeable_block":int(
           (frame.calendar_month_n_cross_cell_exchangeable_blocks>0).sum()),
        "n_taxa_with_any_year_month_crosscell_exchangeable_block":int(
           (frame.calendar_year_month_n_cross_cell_exchangeable_blocks>0).sum()),
    }
    return results,frame


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--photos",type=Path)
    ap.add_argument("--out",type=Path)
    ap.add_argument("--self-test",action="store_true")
    a=ap.parse_args()
    self_test()
    if a.self_test:return
    check(a.photos is not None and a.out is not None,"photo source/output required")
    read=read_original(a.photos)
    dated,dateledger=date_filter(read)
    subsets,opportunity=make_matched_cell_sets(dated)
    results={}
    tabs=[]
    for tag,subset in subsets.items():
        result,taxon=analyze_one(subset,tag)
        results[tag]=result
        taxon.insert(0,"source_collection",tag)
        tabs.append(taxon)
    out={
        "schema":"fcp_world_geographic_colour_sampling_time_conditioned_v1",
        "source_photo_sha256":"8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6",
        "actual_genuine_original_photo_rows":447500,
        "actual_classifiable_photos":227107,
        "source_method":"same four coarse photographic labels, same taxon-cell support thresholds, identical species and cells across photo-purity treatments",
        "date_filter":dateledger,"opportunity":opportunity,
        "permutations":N_REPS,
        "metrics":"between-cell mean four-state TV; absolute white/rest frequency difference; residual nonwhite TV",
        "nulls":"same exact taxa/cells/observer positions/photographic labels; label swaps unrestricted or confined to calendar month or exact year-month",
        "results":results,
        "warnings":[
          "Post hoc test only; original 50km, original prospective white W and earlier IBD/IBE results unchanged",
          "Month-conditioned label exchangeability captures sampling phenology and genuine temporal colour variation together; does not identify either cause",
          "Year-month nonexchangeability can make null uninformative; opportunity denominators are mandatory",
          "Regional photo-colour sorting is not genetic IBD, phenotype fitness, selection, demographic history or breeding-population FCP",
          "Native occurrence, horticulture, flower-age and independent observer-specific image error remain unverified",
          "Multiple related posthoc statistics; no unqualified p-value discovery"
        ]
    }
    a.out.mkdir(parents=True,exist_ok=True)
    (a.out/"fcp_world_geographic_photo_month_null_result.json").write_text(
      json.dumps(out,indent=2,sort_keys=True,allow_nan=False)+"\n")
    pd.concat(tabs,ignore_index=True).to_csv(
        a.out/"fcp_world_geographic_photo_month_aggregate_taxon_no_localities.csv",index=False)
    print("FCP_REAL_WORLD_MONTH_CONDITIONED_SPATIAL",json.dumps(out,sort_keys=True,allow_nan=False))


if __name__=="__main__":main()
