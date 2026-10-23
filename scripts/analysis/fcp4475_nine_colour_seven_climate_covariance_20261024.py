#!/usr/bin/env python3
"""Source-locked 4475-taxon continuous colour/climate covariance decomposition.

Do not equate four photographic source labels to genetic morphs. For measured
photo palette fractions Y and standardized climate Z, the exact finite-sample
identity is cov(Z,Y)=within-original-photo-label cov + between-label cov.

Between component tracks original categorical label frequency/climate sorting.
Within component tracks continuous 9-hue shifts *conditional on that label*.
Neither separates biologically continuous versus discrete inherited morphs,
and any photo-hue/climate signal can be photographic/cultivation/geographic.
A second analysis orthogonalizes climate Z to observed geographic 3D vectors
and calendar month; if climate is nearly geographic-collinear, no conditional
environment conclusion is possible.

All original photographed taxa retained in the manifest and per-species CSV.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "analysis"))
from fcp_full4475_hypothesis_census_20261021 import (
    TOTAL_TAXA, TOTAL_ROWS, CLASSIFIABLE_ROWS, SOURCE_SHA, BIO_ARCHIVE_SHA,
    load_original, sample_climate, all_coords_ok, check, PHOTO_LABELS
)
from fcp_all4475_seven_abiotic_axes_20261021 import (
    add_srad, SRAD_ZIP_SHA, AXES,
)

SEED = 20261024
PERMS = 199
BOOTS = 999
MIN_CLASS = 40
MIN_CELLS = 5
MIN_DISTINCT_OBSERVERS = 20
MIN_GROUP = 5
RESIDUAL_CLIMATE_VARIANCE_MIN = 0.05
NAX = 7
N_HUES = 9


def essential(ok, msg):
    if not ok:
        raise RuntimeError("FCP_COLOUR_CLIMATE_COV_DECOMPOSITION: "+msg)


def design_geo_month(g):
    la = np.deg2rad(g.latitude.to_numpy(float))
    lo = np.deg2rad(g.longitude.to_numpy(float))
    C = np.cos(la)
    xyz = np.stack((C * np.cos(lo), C * np.sin(lo), np.sin(la)), axis=1)
    xyz -= xyz.mean(axis=0)
    span = np.sqrt(np.mean(np.sum(xyz**2, axis=1)))
    xyz /= max(span, 1e-7)
    month = g.date.dt.month.fillna(0).astype(int).to_numpy()
    # Do not put a separate dummy for every observer: at most two photos
    # contributed per observer to the outcome-blind acquisition sample.
    return np.column_stack((
        np.ones(len(g)), xyz,
        np.eye(13, dtype=float)[month, 1:]
    )), month


def exact_partition(palette, z, codes):
    """z centred per photo array shape n x d, palette shape n x 9."""
    palette = np.asarray(palette, float)
    z = np.asarray(z, float)
    n = len(palette)
    assert len(z) == n and palette.shape[1] == 9
    ymean = palette.mean(axis=0)
    y = palette - ymean
    z = z - z.mean(axis=0)
    total = z.T @ y / n
    between = np.zeros_like(total)
    counts = np.bincount(codes, minlength=4)
    for k in range(4):
        if counts[k]:
            mask = codes == k
            between += np.outer(z[mask].sum(axis=0),
                                palette[mask].mean(axis=0) - ymean) / n
    within = total - between
    essential(np.max(np.abs(total - (between + within))) <= 1e-12,
              "law of total covariance failed")
    return total, between, within


def metric_bundle(parts):
    """shape components x seven climate axes x palette hues."""
    data = np.stack(parts)
    return np.stack((
        np.linalg.norm(data, axis=-1).transpose(1,0),
        data[:, :, 0].transpose(1,0)
    ), axis=0)  # two metrics x seven axes x three components


def month_shuffle(months, rng):
    out = np.arange(len(months))
    for m in np.unique(months):
        idx = np.flatnonzero(months == m)
        if len(idx) > 1:
            out[idx] = rng.permutation(idx)
    return out



def fixed_colour_component_basis(Y, codes):
    """Once per source taxon, create exact colour covariance decomposition basis.

    The three n x 9 matrices are total, photographic-label-between and
    photographic-label-within deviations. These sum exactly to total and
    never depend on climate or randomization.
    """
    Y=np.asarray(Y, dtype=float)
    ymean=Y.mean(axis=0)
    labels=np.zeros_like(Y)
    for k in range(4):
        mask = (codes==k)
        if mask.any():
            labels[mask]=Y[mask].mean(axis=0)
    total=Y-ymean
    between=labels-ymean
    within=Y-labels
    essential(np.max(np.abs(total-between-within))<1e-12,
              "fixed photo-label basis fails exact covariance identity")
    return np.column_stack((total,between,within))


def projected_colour_climate_metrics(basis, climate):
    """Species covariance for raw/georesid climate; output reg x metric x axis x comp.

    Repeated null evaluation needs only 14 x n times n x 27 dense matrix
    multiplication; no repeated group-by or normalization of source hues.
    """
    n=basis.shape[0]
    essential(basis.shape==(n,27) and climate.shape==(2,n,NAX),
              "matrix-based source-partition input shape failed")
    centred=climate-climate.mean(axis=1,keepdims=True)
    product=(centred.transpose(0,2,1)@basis/n).reshape(2,NAX,3,N_HUES)
    return np.stack((np.linalg.norm(product,axis=-1), product[:,:,:,0]),axis=1)


def one_species(g, one_observer=False, perms=PERMS):
    species = str(g.species.iloc[0])
    row = {
        "inat_taxon_id": int(g.inat_taxon_id.iloc[0]),
        "species": species,
        "genus": species.split()[0],
        "source_cohort": str(g.fcp_source_cohort.iloc[0]),
        "n_classified_all_original": len(g),
        "n_geo_climate_complete": None,
        "n_cells_10arcmin": None,
        "n_distinct_observers": None,
        "one_photo_per_observer": bool(one_observer),
        "status": "no_eligible_continuous_climate_estimate",
    }
    if len(g) < MIN_CLASS:
        row["status"] = "fewer_than_40_source_classifiable"
        return row,None,None,None,None
    if one_observer:
        g = (g.dropna(subset=["observer_id"])
             .sort_values(["observer_id","photo_id"],kind="stable")
             .drop_duplicates("observer_id")
             .sort_values("photo_id",kind="stable"))
    geo = all_coords_ok(g)
    valid = geo & np.isfinite(g[list(AXES)].to_numpy(float)).all(axis=1)
    g = g.loc[valid].reset_index(drop=True)
    n = len(g)
    row["n_geo_climate_complete"] = n
    if n < MIN_CLASS:
        row["status"] = "fewer_than_40_geovalid_complete_seven_climate"
        return row,None,None,None,None
    # WorldClim 10-arc-minute lattice: photographs in same cell do not form
    # independent climate replicates, regardless of source-photo count.
    lat = g.latitude.to_numpy(float)
    lon = g.longitude.to_numpy(float)
    loc = np.stack((np.floor((90-lat)*6),np.floor((lon+180)*6)),axis=1).astype(int)
    ncell = len(np.unique(loc, axis=0))
    row["n_cells_10arcmin"] = ncell
    nobs = g.observer_id.dropna().nunique()
    row["n_distinct_observers"] = int(nobs)
    if ncell < MIN_CELLS or nobs < MIN_DISTINCT_OBSERVERS:
        row["status"] = "too_few_climate_grid_cells_or_observers"
        return row,None,None,None,None

    Y = np.stack(g.palette_norm9.to_numpy()).astype(float)
    essential(Y.shape == (n,N_HUES) and np.isfinite(Y).all(),
              "unmeasured nine colour vector")
    essential(np.max(np.abs(Y.sum(axis=1)-1)) < 1e-10,
              "normalized palette mass differed from original")
    labels = pd.Categorical(g.morph,categories=PHOTO_LABELS).codes.astype(int)
    essential(np.all(labels>=0),"unrecognized original coarse photographed colour")
    colourcounts = np.bincount(labels,minlength=4)
    row["n_labels_with_five_images"] = int((colourcounts>=5).sum())
    row["n_labels_observed"] = int((colourcounts>0).sum())
    row["frequency_separation_opportunity"] = bool((colourcounts>=MIN_GROUP).sum()>=2)
    env = g[list(AXES)].to_numpy(float)
    zmean = env.mean(axis=0)
    sd = env.std(axis=0)
    varying = sd > np.maximum(1e-7,2e-10*np.maximum(np.abs(zmean),1))
    sd_safe = np.where(varying,sd,1)
    Z = (env-zmean)/sd_safe
    # Mean covariance with climate residualized on simple spatial pattern and
    # calendar month. This measures remaining data variance, not true causal
    # climate variation independent of all geographical history.
    G,months = design_geo_month(g)
    rank = np.linalg.matrix_rank(G)
    essential(rank <= min(n,G.shape[1]),"geo-month rank > possible")
    zfit = G @ np.linalg.lstsq(G,Z,rcond=1e-8)[0]
    R = Z-zfit
    ratio = np.mean(R**2,axis=0)/np.maximum(np.mean(Z**2,axis=0),1e-12)
    row["geo_month_design_rank"] = int(rank)
    row["climate_residual_fraction_median"] = float(np.median(ratio[varying])) if varying.any() else None
    has_resid = varying & (ratio>=RESIDUAL_CLIMATE_VARIANCE_MIN)
    row["axes_with_noncollinear_geo_residual"] = int(has_resid.sum())
    available = np.stack((varying,has_resid),axis=0)
    # raw vs geography/month-residual standardized source-climate deviations
    both=np.stack((Z,R),axis=0)
    basis=fixed_colour_component_basis(Y,labels)
    observed=projected_colour_climate_metrics(basis,both)
    rng = np.random.default_rng(SEED+int(row["inat_taxon_id"])*53+int(one_observer))
    null=np.empty((perms,2,2,NAX,3),float)
    for j in range(perms):
        perm=month_shuffle(months,rng)
        null[j]=projected_colour_climate_metrics(basis,both[:,perm,:])
    # Within and between magnitudes cannot be interpreted as additive shares
    # because vectors may cancel. Do not take a component ratio as variance
    # explained or as proportion of genetic morph frequency change.
    for regime,label in enumerate(("raw","geo_month_residual")):
        for axis,name in enumerate(AXES):
            row[f"{name}_{label}_eligible"]=bool(available[regime,axis])
            row[f"{name}_{label}_residual_fraction"]=float(ratio[axis]) if varying[axis] else None
            if not available[regime,axis]:
                continue
            for component,j in (("total",0),("between_coarse_label",1),("within_coarse_label",2)):
                if j == 1 and not row["frequency_separation_opportunity"]:
                    # Do not infer photographed between-label climate
                    # sorting from one or extremely sparse label.
                    continue
                for metric,m in (("hue_vector_norm",0),("signed_white_fraction",1)):
                    val=float(observed[regime,m,axis,j])
                    row[f"{name}_{label}_{component}_{metric}"]=val
                    row[f"{name}_{label}_{component}_{metric}_excess_vs_month_null"]=float(
                        val - null[:,regime,m,axis,j].mean())
    row["status"]="available_seven_abiotic_axis_decomposition"
    return row, observed, null, available, n


def fold_species(rows,observed,nulls,eligible,one_observer):
    # All species enter no-source-exclusion ledger; no zeroes for ineligible
    # taxa, including colour-specific 2+ supported class opportunities.
    counts = np.zeros((2,NAX,3),int)
    normObs=np.zeros((2,NAX,3),float)
    normNull=np.zeros((PERMS,2,NAX,3),float)
    signedObs=np.zeros((2,NAX,3),float)
    signedNull=np.zeros((PERMS,2,NAX,3),float)
    used=[]
    for row,obs,nul,ok in zip(rows,observed,nulls,eligible):
        if obs is None:continue
        used.append(row)
        for regime in range(2):
            for ax in range(NAX):
                if not ok[regime,ax]:continue
                for component in range(3):
                    if component == 1 and not row["frequency_separation_opportunity"]:
                        continue
                    counts[regime,ax,component]+=1
                    normObs[regime,ax,component]+=obs[regime,0,ax,component]
                    signedObs[regime,ax,component]+=obs[regime,1,ax,component]
                    normNull[:,regime,ax,component]+=nul[:,regime,0,ax,component]
                    signedNull[:,regime,ax,component]+=nul[:,regime,1,ax,component]
    out={}
    for regime,name in enumerate(("raw","geo_month_residual")):
        out[name]={}
        for ax,axis in enumerate(AXES):
            count=counts[regime,ax]
            n=int(count[0])
            if not n:
                out[name][axis]={"n_species":0,"status":"not_estimable"}
                continue
            results={}
            for j,component in enumerate(("total","between_photo_labels","within_photo_labels")):
                n_component=int(count[j])
                if n_component == 0:
                    results[component]={"n_species":0,"status":"insufficient_colour_label_support"}
                    continue
                ob=normObs[regime,ax,j]/n_component
                nu=normNull[:,regime,ax,j]/n_component
                wobs=signedObs[regime,ax,j]/n_component
                wnull=signedNull[:,regime,ax,j]/n_component
                pnorm=float((1+np.count_nonzero(nu>=ob))/(PERMS+1))
                # True distribution of source photos, signed scalar white
                # loading may be positive or negative; use two-sided null.
                pwhite=float((1+np.count_nonzero(np.abs(wnull-wnull.mean())
                                              >=abs(wobs-wnull.mean())))/(PERMS+1))
                results[component]={
                    "n_species":n_component,
                    "observed_norm_mean":float(ob),
                    "month_null_norm_mean":float(nu.mean()),
                    "hue_vector_norm_excess_over_month_null":float(ob-nu.mean()),
                    "hue_norm_permutation_upper_p_199":pnorm,
                    "observed_signed_white_fraction_change_per_climate_SD":float(wobs),
                    "month_null_signed_white_mean":float(wnull.mean()),
                    "signed_white_two_sided_p_199":pwhite,
                }
            out[name][axis]={"n_species":n,"n_between_label_supported":int(count[1]),"components":results}
    return out


def holm(pvals):
    ps=np.array(pvals,float)
    ix=np.argsort(ps,kind="stable")
    out=np.zeros_like(ps)
    now=0.
    for i,j in enumerate(ix):
        now=max(now,min(1.,(len(ps)-i)*ps[j]))
        out[j]=now
    return [float(x) for x in out]


def add_multiplicity(summary):
    for regime in ("raw","geo_month_residual"):
        for comp in ("total","between_photo_labels","within_photo_labels"):
            for metric in ("hue_norm_permutation_upper_p_199","signed_white_two_sided_p_199"):
                ps=[]
                for ax in AXES:
                    row=summary[regime][ax]
                    ps.append(row.get("components",{}).get(comp,{}).get(metric,1))
                adjusted=holm(ps)
                for ax,p in zip(AXES,adjusted):
                    x=summary[regime][ax]
                    if x["n_species"] and x["components"][comp].get("n_species",0):
                        x["components"][comp][metric+"_holm7"]=p


def bootstrap_genus(per_taxa,summary,one_observer,rep=BOOTS):
    """Exact genus clusters aggregated with vectorized sufficient statistics."""
    d=pd.DataFrame(per_taxa)
    codes,genus_names=pd.factorize(d.genus,sort=True)
    G=len(genus_names)
    rng=np.random.default_rng(SEED+100+int(one_observer))
    draws=rng.integers(G,size=(rep,G))
    results={}
    for regime in ("raw","geo_month_residual"):
        for component in ("total","between_coarse_label","within_coarse_label"):
            for measure in ("hue_vector_norm","signed_white_fraction"):
                cols=[f"{ax}_{regime}_{component}_{measure}_excess_vs_month_null"
                      for ax in AXES]
                vals=d.reindex(columns=cols).to_numpy(float)
                keep=np.isfinite(vals)
                sums=np.zeros((G,NAX))
                nums=np.zeros((G,NAX))
                np.add.at(sums,codes,np.where(keep,vals,0))
                np.add.at(nums,codes,keep.astype(float))
                boot_sum=sums[draws].sum(axis=1)
                boot_n=nums[draws].sum(axis=1)
                estimates=np.divide(boot_sum,boot_n,out=np.full_like(boot_sum,np.nan),where=boot_n>0)
                names={}
                for j,axis in enumerate(AXES):
                    cell=summary[regime][axis].get("components",{}).get(
                        {"total":"total","between_coarse_label":"between_photo_labels",
                         "within_coarse_label":"within_photo_labels"}[component],{})
                    if cell.get("n_species",0)>0:
                        v=estimates[:,j]
                        v=v[np.isfinite(v)]
                        if len(v)>rep*.9:
                            names[axis]={"genus_bootstrap_95ci":[float(x) for x in np.quantile(v,[.025,.975])],
                                         "n_bootstrap_draws":int(len(v))}
                results[regime+"_"+component+"_"+measure]=names
    return results,G


def main():
    p=argparse.ArgumentParser()
    for k in ("photos","worldclim-bio","worldclim-srad","out"):
        p.add_argument("--"+k,type=Path)
    p.add_argument("--self-test",action="store_true")
    p.add_argument("--permutations",type=int,default=PERMS)
    args=p.parse_args()
    if args.self_test:
        y=np.array([[1,0,0,0,0,0,0,0,0],
                    [1,0,0,0,0,0,0,0,0],
                    [0,1,0,0,0,0,0,0,0],
                    [0,1,0,0,0,0,0,0,0]],float)
        z=np.array([[-1,0],[-1,2],[1,1],[1,3]],float)
        t,b,w=exact_partition(y,z,np.array([0,0,1,1]))
        essential(np.max(abs(t-b-w))<1e-14 and abs(w[0,0])<1e-14,
                  "between/within decomposition not exact")
        y2=np.array([[.9,.1]+[0]*7,[.7,.3]+[0]*7,
                     [.1,.9]+[0]*7,[.3,.7]+[0]*7],float)
        t,b,w=exact_partition(y2,z,np.array([0,0,1,1]))
        essential(np.max(abs(t-b-w))<1e-14,
                  "continuous hue shifts lost covariance identity")
        colour_basis=fixed_colour_component_basis(y2,np.array([0,0,1,1]))
        zz=np.zeros((2,4,NAX),float)
        zz[:,:,0]=np.array([-1,0,1,2],float)
        calculated=projected_colour_climate_metrics(colour_basis,zz)
        independent=exact_partition(y2,zz[0],np.array([0,0,1,1]))
        essential(np.allclose(calculated[0,0,0,:],
                  [np.linalg.norm(x[0]) for x in independent],atol=1e-12),
                  "fast source covariance differs from independent exact fixture")
        q=month_shuffle(np.array([2,2,3,3]),np.random.default_rng(SEED))
        essential(set(q[:2])=={0,1} and set(q[2:])=={2,3},
                  "month-conditional permutation leaked among months")
        # Artificial same-month all-photo labels anchored to static geography
        # must be allowed in full Z shuffles but not declared independent
        # month effects when month has no eligible swaps.
        essential(holm([.01,.05,.2])[0]<=.03,"Holm multiplicity broken")
        print("FCP4475_CLIMATE_HUE_COVARIANCE_SELF_TEST_PASS")
        return
    essential(args.permutations==PERMS,"permutation count must remain 199 for CI comparability")
    essential(args.photos and args.worldclim_bio and args.worldclim_srad and args.out,
              "source photo and both WorldClim original archives required")
    source,photos=load_original(args.photos)
    photos=sample_climate(photos,args.worldclim_bio)
    photos=add_srad(photos,args.worldclim_srad)
    ledger={}
    for thinned in (False,True):
        label="observer_thinned" if thinned else "original"
        rows=[];allobs=[];allnull=[];allok=[]
        for ix,(tid,g) in enumerate(photos.groupby("inat_taxon_id",sort=True)):
            row,obs,nul,ok,_ = one_species(g,thinned)
            rows.append(row);allobs.append(obs);allnull.append(nul);allok.append(ok)
            if (ix+1)%1000==0:
                print("FCP_CLIMATE_CONTINUOUS_TAXA",label,ix+1,flush=True)
        essential(len(rows)==TOTAL_TAXA and
                  sum(r["n_classified_all_original"] for r in rows)==CLASSIFIABLE_ROWS,
                  "all 4475/227107 source counts not preserved")
        summary=fold_species(rows,allobs,allnull,allok,thinned)
        add_multiplicity(summary)
        genus_ci,genus_n=bootstrap_genus(rows,summary,thinned)
        table=pd.DataFrame(rows)
        args.out.mkdir(parents=True,exist_ok=True)
        table.to_csv(args.out/f"all4475_{label}_photo9_climate7_per_taxon.csv",index=False)
        ledger[label]={
            "total_species_in_source":len(rows),
            "n_genus_bootstrap":genus_n,
            "n_species_any_raw_climate_axis":int(table.status.eq("available_seven_abiotic_axis_decomposition").sum()),
            "n_species_at_least_2_coarse_photo_labels_with_5_each":int(
                (table.status.eq("available_seven_abiotic_axis_decomposition")&
                 table.frequency_separation_opportunity.fillna(False)).sum()),
            "per_axis_results":summary,
            "genus_bootstrap":genus_ci,
            "species_denominator_diagnostics":table.status.value_counts().to_dict(),
        }
    result={
        "schema":"fcp_original_photo447500_nine_palette_7climate_exact_cov_decomposition_v1",
        "status":"source_run_completed_posthoc_compositional_covariance_not_genetic_morphs",
        "date_jst":"2026-10-24",
        "source_original_photos":TOTAL_ROWS,
        "source_original_species":TOTAL_TAXA,
        "source_classifiable_photos":CLASSIFIABLE_ROWS,
        "photo_source_sha256":SOURCE_SHA,
        "bioclim_archive_sha256":BIO_ARCHIVE_SHA,
        "solar_archive_sha256":SRAD_ZIP_SHA,
        "axes":list(AXES),
        "n_month_matched_permutations":PERMS,
        "genus_bootstrap_iterations":BOOTS,
        "photo_sampling":"original classified flower ROI nine-palette; all 4475 taxa retained; no missing hue imputation",
        "raw_covariance_identity":"total=between_original_four_photo_state_labels+within_original_four_photo_state_labels exactly for each species and environment",
        "signed_white_score":"covariance of continuous white-palette fraction with species-standardized observed climate; not genetic white morph frequency",
        "hue_vector_norm":"norm of nine palette/climate covariance after subtracting month-permutation norm expectation; cannot interpret within/between magnitudes as additive causal variance shares",
        "geo_month_residual":"calendar month and 3D Earth-centred coordinate linear projection on standardized climate; retain only >=5% original climate variance after projection; not a full nonlinear space/spatiotemporal correction",
        "climate_cells":"at least 5 occupied 10-arcminute WorldClim cells; 40 eligible photos; at least 20 observers",
        "multiplicity":"Holm 7 climate axes within each post hoc endpoint family, 199 permutation worlds; many endpoint families remain exploratory",
        "temperate_vs_hot_slope_unproved":True,
        "continuous_vs_discrete_genetic_morphs_not_identified":True,
        "camera_organ_background_month_bias_possible":True,
        "no_actual_fitness_genomics_or_heritable_morph_measurement":True,
        "results":ledger,
    }
    (args.out/"summary.json").write_text(
        json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("FCP4475_9COLOUR_7CLIMATE_COVARIANCE_COMPLETE",json.dumps({
        kind:{
            "eligible":v["n_species_any_raw_climate_axis"],
            "multi_supported":v["n_species_at_least_2_coarse_photo_labels_with_5_each"],
            "summary":{reg:{ax:{"n":a["n_species"],"total":a.get("components",{}).get("total")}
                         for ax,a in val.items()}
                       for reg,val in v["per_axis_results"].items()}}
        for kind,v in ledger.items()},allow_nan=False))
if __name__=="__main__":
    main()
