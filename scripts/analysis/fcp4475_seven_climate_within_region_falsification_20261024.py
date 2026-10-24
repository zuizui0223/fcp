#!/usr/bin/env python3
"""Within-region falsification of 447,500-source FCP photographed hue–climate.

The previously authenticated seven-climate temperature/white association used
linear XYZ geographic controls. This independent analysis replaces that with
FIXED 5-degree and 1-degree geographic blocks and a species calendar-month
term, then compares original photo hues with 199 spatial-block-and-month
permutations. One-photo-per-observer sensitivity is independent of photo hues.

This is PHOTO covariance only, not adaptation or genetic morph frequencies.
The 4 photographic colour labels enter only an exact algebraic covariance
partition, never a genotype or inherited flower colour assertion.
"""
from __future__ import annotations

import argparse
import json
import hashlib
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"/"analysis"))
from fcp_full4475_hypothesis_census_20261021 import (
    load_original, sample_climate, all_coords_ok, TOTAL_TAXA,
    TOTAL_ROWS, CLASSIFIABLE_ROWS, SOURCE_SHA, BIO_ARCHIVE_SHA, PHOTO_LABELS
)
from fcp_all4475_seven_abiotic_axes_20261021 import add_srad, SRAD_ZIP_SHA, AXES
from fcp4475_nine_colour_seven_climate_covariance_20261024 import (
    fixed_colour_component_basis
)

SCALES=(5,1)
PERMS=199
BOOTS=499
SEED=20261024
MIN_CLASS=40
MIN_OBSERVERS=20
MIN_RASTER_CELLS=5
MIN_INFORMATIVE_PHOTOS=12
MIN_INFORMATIVE_STRATA=2
MIN_CLIMATE_RESID_FRAC=.05
COMP=("total","between_original_photo_labels","within_original_photo_labels")
NAX=len(AXES)


def check(yes,why):
    if not yes:raise RuntimeError("FCP_PHOTO_CLIMATE_REGIONAL_FALSIFICATION: "+why)


def climate_block(g,width):
    lat=g.latitude.to_numpy(float)
    lon=g.longitude.to_numpy(float)
    blocks=np.column_stack((
        np.floor((lat+90)/width),np.floor((lon+180)/width)
    )).astype(int)
    _, bidx=np.unique(blocks,axis=0,return_inverse=True)
    return bidx


def residualize_fixed_blocks_month(climate,block,month):
    n=len(block)
    nb=int(block.max())+1
    design=np.concatenate((
        np.eye(nb,dtype=float)[block],
        np.eye(13,dtype=float)[month,1:]
    ),axis=1)
    coeff=np.linalg.lstsq(design,climate,rcond=1e-9)[0]
    remaining=climate-design@coeff
    check(remaining.shape==climate.shape,
          "geography/month fixed effect wrong shape")
    return remaining


def permutable_groups(g,block,month):
    """Null: exchange measured photos only within block and calendar month."""
    strata=block*13+month
    groups=[np.flatnonzero(strata==k)
            for k in np.unique(strata)]
    return groups


def informative_support(env,rows,groups):
    """How many photo observations actually have within-stratum climate contrasts?"""
    info=np.zeros((env.shape[1],),int)
    groups_count=np.zeros_like(info)
    cell_counts=[]
    for idx in groups:
        if len(idx)<2:continue
        distinct=np.unique(rows[idx],axis=0).shape[0]
        if distinct<2:continue
        for a in range(env.shape[1]):
            if np.ptp(env[idx,a])>1e-7:
                info[a]+=len(idx)
                groups_count[a]+=1
        cell_counts.append(distinct)
    return info,groups_count


def one_taxon(original,one_observer,width,seed=SEED,perms=PERMS):
    name=str(original.species.iloc[0])
    row={
        "inat_taxon_id":int(original.inat_taxon_id.iloc[0]),
        "species":name,"genus":name.split()[0],
        "source_cohort":str(original.fcp_source_cohort.iloc[0]),
        "n_source_classifiable":int(len(original)),
        "one_photo_per_observer":bool(one_observer),
        "block_width_degree":int(width),
        "status":"not_estimable",
        "n_photos_with_complete_climate":0,
        "n_source_observers":0,
        "n_worldclim_cells":0,
        "n_regional_blocks":0,
        "n_mixed_photo_label_blocks":0,
    }
    if len(original)<MIN_CLASS:
        row["status"]="fewer_than_40_classified"
        return row,None,None,None

    g=original
    if one_observer:
        g=(g.dropna(subset=["observer_id"])
             .sort_values(["observer_id","photo_id"],kind="stable")
             .drop_duplicates("observer_id")
             .sort_values("photo_id",kind="stable"))
    climate=g[list(AXES)].to_numpy(float)
    valid=all_coords_ok(g)&np.isfinite(climate).all(axis=1)&g.date.notna().to_numpy()
    g=g.loc[valid].reset_index(drop=True)
    n=len(g)
    row["n_photos_with_complete_climate"]=int(n)
    row["n_source_observers"]=int(g.observer_id.nunique())
    if n<MIN_CLASS or g.observer_id.nunique()<MIN_OBSERVERS:
        row["status"]="fewer_than_40_complete_or_20_observers"
        return row,None,None,None
    coords=np.column_stack((
        np.floor((90-g.latitude.to_numpy(float))*6),
        np.floor((g.longitude.to_numpy(float)+180)*6)
    )).astype(int)
    nc=np.unique(coords,axis=0).shape[0]
    row["n_worldclim_cells"]=int(nc)
    if nc<MIN_RASTER_CELLS:
        row["status"]="fewer_than_5_original_climate_grid_cells"
        return row,None,None,None
    block=climate_block(g,width)
    row["n_regional_blocks"]=int(block.max()+1)
    mon=g.date.dt.month.to_numpy(int)
    Z=g[list(AXES)].to_numpy(float)
    means=Z.mean(axis=0)
    sd=Z.std(axis=0)
    climatic=sd>np.maximum(1e-7,2e-10*np.maximum(np.abs(means),1))
    Z=(Z-means)/np.where(climatic,sd,1)
    R=residualize_fixed_blocks_month(Z,block,mon)
    ratio=np.mean(R**2,axis=0)/np.maximum(np.mean(Z**2,axis=0),1e-14)
    groups=permutable_groups(g,block,mon)
    info,groupcounts=informative_support(Z,coords,groups)
    row["n_source_photo_strata"]=int(len(groups))
    row["n_swappable_photos"]=int(sum(len(a) for a in groups if len(a)>1))
    row["n_informative_BIO5_photos"]=int(info[2])
    row["n_informative_BIO5_strata"]=int(groupcounts[2])
    row["BIO5_within_region_climate_residual_variance_fraction"]=float(ratio[2]) if climatic[2] else None
    available=(climatic
               &(ratio>=MIN_CLIMATE_RESID_FRAC)
               &(info>=MIN_INFORMATIVE_PHOTOS)
               &(groupcounts>=MIN_INFORMATIVE_STRATA))
    photo=pd.Categorical(g.morph,categories=PHOTO_LABELS).codes
    check(np.all(photo>=0),"unrecognized original original-photo category")
    counts=np.bincount(photo,minlength=4)
    well_supported=int((counts>=5).sum())>=2
    # For group-frequency interpretation, actually observe >=2 distinct
    # original photo labels in a region with >=2 source photos per label.
    mixed=0
    for b in np.unique(block):
        c=np.bincount(photo[block==b],minlength=4)
        mixed+=int((c>=2).sum()>=2)
    row["n_mixed_photo_label_blocks"]=int(mixed)
    between_allowed=well_supported and mixed>=1
    y=np.stack(g.palette_norm9.to_numpy()).astype(float)
    check(y.shape==(n,9) and np.allclose(y.sum(axis=1),1,rtol=0,atol=1e-10),
          "original photo nine-channel source normalization changed")
    basis=fixed_colour_component_basis(y,photo)
    cov=(R.T@basis/n).reshape(NAX,3,9)
    observed=np.stack([cov[:,:,0],np.linalg.norm(cov,axis=-1)],axis=0)
    null=np.empty((perms,2,NAX,3),float)
    rng=np.random.default_rng(seed+int(row["inat_taxon_id"])*31
                               +int(one_observer)*19+width)
    for i in range(perms):
        order=np.arange(n)
        for idx in groups:
            if len(idx)>1:
                order[idx]=rng.permutation(idx)
        C=(R[order].T@basis/n).reshape(NAX,3,9)
        null[i,0]=C[:,:,0]
        null[i,1]=np.linalg.norm(C,axis=-1)
    null_mean=null.mean(axis=0)
    for a,axis in enumerate(AXES):
        row[f"{axis}_climate_resid_variance_fraction"]=float(ratio[a]) if climatic[a] else None
        row[f"{axis}_n_informative_photos"]=int(info[a])
        row[f"{axis}_n_informative_block_month_groups"]=int(groupcounts[a])
        row[f"{axis}_available"]=bool(available[a])
        if not available[a]:continue
        for j,cname in enumerate(COMP):
            if j==1 and not between_allowed:continue
            for z,metric in enumerate(("signed_white","hue_vector_norm")):
                val=float(observed[z,a,j])
                row[f"{axis}_{cname}_{metric}"]=val
                row[f"{axis}_{cname}_{metric}_excess"]=float(val-null_mean[z,a,j])
    row["status"]="regional_within_block_test_run"
    # between component per taxon needs regional mixed original labels
    allowed=np.stack([available,available&between_allowed,available],axis=1)
    return row,observed,null,allowed


def summarize(rows,observers,nulls,allowed,scale,thin):
    counts=np.zeros((NAX,3),int)
    obs=np.zeros((2,NAX,3),float)
    nul=np.zeros((PERMS,2,NAX,3),float)
    for row,x,y,mask in zip(rows,observers,nulls,allowed):
        if x is None:continue
        for a in range(NAX):
            for j in range(3):
                if not mask[a,j]:continue
                counts[a,j]+=1
                obs[:,a,j]+=x[:,a,j]
                nul[:,:,a,j]+=y[:,:,a,j]
    result={}
    for a,axis in enumerate(AXES):
        comps={}
        for j,cname in enumerate(COMP):
            n=int(counts[a,j])
            if n==0:
                comps[cname]={"n_species":0,"status":"not_identifiable"}
                continue
            o=obs[:,a,j]/n
            p=nul[:,:,a,j]/n
            signed0=float(o[0]-p[:,0].mean())
            signedp=float((1+np.count_nonzero(
                abs(p[:,0]-p[:,0].mean())>=abs(signed0)))/(PERMS+1))
            norm0=float(o[1]-p[:,1].mean())
            normp=float((1+np.count_nonzero(p[:,1]>=o[1]))/(PERMS+1))
            comps[cname]={
                "n_species":n,
                "signed_white_covariance_excess_vs_block_month_null":signed0,
                "signed_two_sided_p_199":signedp,
                "hue_vector_covariance_norm_excess_vs_block_month_null":norm0,
                "norm_upper_p_199":normp,
                "observed_white_covariance":float(o[0]),
                "null_white_covariance":float(p[:,0].mean()),
            }
        result[axis]={
            "n_total_support":int(counts[a,0]),
            "n_between_label_mixed_within_blocks":int(counts[a,1]),
            "components":comps,
        }
    return result


def holm(ps):
    vals=np.array(ps,float)
    ix=np.argsort(vals)
    out=np.zeros_like(vals)
    running=0.
    for i,j in enumerate(ix):
        running=max(running,min(1.,(len(vals)-i)*vals[j]))
        out[j]=running
    return out


def add_holm(panel):
    for j,cname in enumerate(COMP):
        for key in ("signed_two_sided_p_199","norm_upper_p_199"):
            ps=[panel[axis]["components"][cname].get(key,1)
                for axis in AXES]
            corrected=holm(ps)
            for axis,p in zip(AXES,corrected):
                if panel[axis]["components"][cname]["n_species"]:
                    panel[axis]["components"][cname][key+"_holm7"]=float(p)


def genus_ci(frame,panel,reps=BOOTS,seed=SEED):
    # Data have zero or more per-axis effect columns; never impute absent
    # species as having a zero environmental association.
    d=pd.DataFrame(frame)
    genus,names=pd.factorize(d.genus,sort=True)
    G=len(names)
    rng=np.random.default_rng(seed)
    indices=rng.integers(G,size=(reps,G))
    out={}
    for cname in COMP:
        for metric in ("signed_white","hue_vector_norm"):
            cols=[f"{axis}_{cname}_{metric}_excess" for axis in AXES]
            val=d.reindex(columns=cols).to_numpy(float)
            valid=np.isfinite(val)
            sumG=np.zeros((G,NAX),float)
            countG=np.zeros((G,NAX),float)
            np.add.at(sumG,genus,np.where(valid,val,0))
            np.add.at(countG,genus,valid.astype(float))
            bootsum=sumG[indices].sum(axis=1)
            bootn=countG[indices].sum(axis=1)
            boot=np.divide(bootsum,bootn,out=np.full_like(bootsum,np.nan),where=bootn>0)
            for i,axis in enumerate(AXES):
                if not panel[axis]["components"][cname]["n_species"]:continue
                samples=boot[:,i]
                samples=samples[np.isfinite(samples)]
                if len(samples)<.95*reps:continue
                out.setdefault(axis,{})[cname+"_"+metric]={
                    "marginal_genus95":[float(v) for v in np.quantile(samples,[.025,.975])],
                    "genus_bootstrap_draws":int(len(samples))}
    return out,G


def self_test():
    y=np.array([
        [.9,.1]+[0]*7,[.9,.1]+[0]*7,[.9,.1]+[0]*7,[.9,.1]+[0]*7,
        [.1,.9]+[0]*7,[.1,.9]+[0]*7,[.1,.9]+[0]*7,[.1,.9]+[0]*7
    ],float)
    z=np.array([-2,-1,1,2,-2,-1,1,2],float)
    block=np.array([0,0,0,0,1,1,1,1],int)
    month=np.repeat([6],8)
    resid=residualize_fixed_blocks_month(z[:,None],block,month).ravel()
    a=fixed_colour_component_basis(y,np.array([0]*4+[1]*4))
    check(abs(float(resid@a[:,0]/8))<1e-12,
          "confounded colour versus region did not vanish within regions")
    y2=y.copy()
    y2[:,0]=.5+z*.1
    y2[:,1]=1-y2[:,0]
    b=fixed_colour_component_basis(y2,np.array([0]*4+[1]*4))
    check(float(resid@b[:,0]/8)>.005,
          "true within-region color/climate shift was not detected")
    check(not np.any(holm([.01,.04,.2])<0),"Holm correction invalid")
    print("FCP_CLIMATE_WITHIN_REGION_SYNTHETIC_POSITIVE_AND_CONFOUND_NEGATIVE_PASS")


def main():
    p=argparse.ArgumentParser()
    for field in ("photos","worldclim-bio","worldclim-srad","out"):
        p.add_argument("--"+field,type=Path)
    p.add_argument("--self-test",action="store_true")
    args=p.parse_args()
    self_test()
    if args.self_test:return
    check(args.photos and args.worldclim_bio and args.worldclim_srad and args.out,
          "all three SHA-locked sources and output directory required")
    source, classified=load_original(args.photos)
    classified=sample_climate(classified,args.worldclim_bio)
    classified=add_srad(classified,args.worldclim_srad)
    check(len(source)==TOTAL_ROWS and len(classified)==CLASSIFIABLE_ROWS,
          "original 447500-photo source denominator changed")
    args.out.mkdir(parents=True,exist_ok=True)
    products={}
    for scale in SCALES:
        for thin in (False,True):
            name=("observer_thinned" if thin else "original")
            rows=[]; obs=[]; nul=[]; allowed=[]
            for j,(tid,g) in enumerate(classified.groupby("inat_taxon_id",sort=True)):
                row,x,y,mask=one_taxon(g,thin,scale)
                rows.append(row);obs.append(x);nul.append(y);allowed.append(mask)
                if (j+1)%1000==0:
                    print("FCP_WITHIN_REGION_TAXA_SCANNED",scale,name,j+1,flush=True)
            check(len(rows)==TOTAL_TAXA and
                  sum(row["n_source_classifiable"] for row in rows)==CLASSIFIABLE_ROWS,
                  "missing original taxa or classified denominator")
            panel=summarize(rows,obs,nul,allowed,scale,thin)
            add_holm(panel)
            intervals,genera=genus_ci(rows,panel,seed=SEED+scale+int(thin)*13)
            for axis,entry in panel.items():
                entry["genus_intervals"]=intervals.get(axis,{})
            d=pd.DataFrame(rows)
            dest=args.out/f"original_all4475_{name}_{scale}degree_regional_photo9_climate7.csv"
            d.to_csv(dest,index=False,float_format="%.12g")
            products[f"{scale}degree_{name}"]={
                "n_original_species_in_ledger":len(d),
                "n_taxa_having_any_regional_test":int(d.status.eq("regional_within_block_test_run").sum()),
                "n_BIO5_informative":panel["BIO5"]["n_total_support"],
                "n_BIO5_mixed_photo_label":panel["BIO5"]["n_between_label_mixed_within_blocks"],
                "n_genera_bootstrap":genera,
                "failure_status_counts":d.status.value_counts().to_dict(),
                "seven_abiotic_axes":panel,
            }
    result={
        "schema":"fcp4475_original_source_photo9_climate7_region_month_fixed_falsification_v1",
        "source_original_447500_sha256":SOURCE_SHA,
        "worldclim_bio_archive_sha256":BIO_ARCHIVE_SHA,
        "worldclim_srad_archive_sha256":SRAD_ZIP_SHA,
        "source_n_species":TOTAL_TAXA,
        "source_n_photos":TOTAL_ROWS,
        "source_n_4state_classified":CLASSIFIABLE_ROWS,
        "source_photo_colour_bins_are_not_genetic_morphs":True,
        "region_width_degrees":list(SCALES),
        "region_coordinates":"fixed latitude/longitude degree grids, not equal-area; geolocation precision <=5 km",
        "controls":"separate within-species region fixed effects and month fixed effects; permute original photo climate residual within fixed region × calendar month; preserve category labels; 1 photo per observer separate",
        "eligible_min_original_photos":MIN_CLASS,
        "eligible_min_observers":MIN_OBSERVERS,
        "eligible_min_10arcmin_raster_cells":MIN_RASTER_CELLS,
        "minimum_informative_photos":MIN_INFORMATIVE_PHOTOS,
        "minimum_informative_region_month_groups":MIN_INFORMATIVE_STRATA,
        "minimum_within_region_climate_residual_variance_fraction":MIN_CLIMATE_RESID_FRAC,
        "original_photo_labels_need_two_labels_each_five_for_between_and_two_each_locally":True,
        "photo_months_pooled_across_years":True,
        "photographic_negative_background_not_independently_controlled":True,
        "conditional_on_observed_uneven_sampling_not_randomized_climate_or_photography":True,
        "ineligible_are_not_zeros":True,
        "taxa_weights_equal_not_photo_pair_pseudoreplication":True,
        "bootstrap_genera_draws":BOOTS,
        "month_region_null_draws":PERMS,
        "seven_axis_holm_within_each_exploratory_endpoint_only":True,
        "earlier_strict_5km_same_year_month_photographic_null_remains_unsupported":True,
        "no_causal_temperature_or_genetic_FCP_adaptation_inference":True,
        "results":products
    }
    (args.out/"region_month_climate_falsification_summary.json").write_text(
        json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("FCP4475_REGIONAL_CLIMATE_FALSIFICATION_COMPLETE",
          json.dumps({k:{"n_taxa":v["n_taxa_having_any_regional_test"],
                         "bio5":v["seven_abiotic_axes"]["BIO5"]}
                      for k,v in products.items()},sort_keys=True,allow_nan=False))


if __name__=="__main__":
    main()
