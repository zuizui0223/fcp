#!/usr/bin/env python3
"""FCP authentic 447,500 photographs: climate *domains* and nine-hue directions.

The correlated seven WorldClim predictors are represented by five fixed,
PHOTO-OUTCOME-BLIND coordinates:
 thermal_level = mean(species-z(BIO1), species-z(BIO5))
 precip_amount = mean(species-z(BIO12), species-z(BIO14))
 BIO4 temperature seasonality, BIO15 precipitation seasonality, SRADmean.
These five source domains are mutually adjusted after broad XYZ geographic
and calendar-month controls. This intentionally tests a GROUP-level shared
thermal/moisture gradient, NOT seven individually independent causal drivers.

No biological petal colour, heritable morph, local adaptation, biochemical
pigment investment or flower-specific background correction is inferred.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"/"analysis"))
from fcp_full4475_hypothesis_census_20261021 import (
    load_original, sample_climate, all_coords_ok,
    SOURCE_SHA, BIO_ARCHIVE_SHA, TOTAL_ROWS, TOTAL_TAXA, CLASSIFIABLE_ROWS,
    PHOTO_LABELS, WCOLOURS
)
from fcp_all4475_seven_abiotic_axes_20261021 import (
    add_srad, SRAD_ZIP_SHA, AXES
)
from fcp4475_nine_colour_seven_climate_covariance_20261024 import (
    design_geo_month, fixed_colour_component_basis,
    MIN_CLASS, MIN_CELLS, MIN_DISTINCT_OBSERVERS, MIN_GROUP
)

DOMAINS=("thermal_level","precip_amount","BIO4_seasonality",
         "BIO15_seasonality","SRADmean")
COMP=("total","between_original_photo_labels","within_original_photo_labels")
K=len(DOMAINS)
H=len(WCOLOURS)
B=999
SEED=20261024
CLIMATE_MIN_REMAIN=.05
Q=np.array([8.]+[-1.]*8)/np.sqrt(72.)


def must(ok,why):
    if not ok:
        raise RuntimeError("FCP4475_SHARED_CLIMATE_DOMAIN: "+why)


def compose_macros(z):
    """Outcome-blind signed 5D macro-climate basis for original seven normals."""
    z=np.asarray(z,float)
    must(z.ndim==2 and z.shape[1]==7,"source predictor axis shape")
    m=np.column_stack((
        (z[:,0]+z[:,2])/2,
        (z[:,3]+z[:,4])/2,
        z[:,1],z[:,5],z[:,6],
    ))
    return m


def one_month_null(x,basis,months):
    """Exact expectation of a one-variable within-month source-photo shuffle."""
    n=len(x)
    must(basis.shape==(n,3*H) and len(months)==n,
         "photo hue/month null expectation shape mismatch")
    v=np.zeros((3*H,),float)
    for m in np.unique(months):
        ids=np.flatnonzero(months==m)
        v+=(x[ids].sum()/n)*basis[ids].mean(axis=0)
    return v.reshape((3,H))


def one_species(g,observer_thinned=False):
    species=str(g.species.iloc[0])
    row={
       "inat_taxon_id":int(g.inat_taxon_id.iloc[0]),
       "species":species,
       "genus":species.split()[0],
       "source_cohort":str(g.fcp_source_cohort.iloc[0]),
       "n_original_classifiable":int(len(g)),
       "n_photo_climate_complete":0,
       "n_original_observers":0,
       "n_original_raster_cells":0,
       "n_photo_bins_ge5":0,
       "one_photo_per_observer":bool(observer_thinned),
       "status":"not_estimable",
    }
    empty=np.full((K,3,H),np.nan)
    if len(g)<MIN_CLASS:
        row["status"]="less_than_40_classifiable"
        return row,empty,empty.copy()
    if observer_thinned:
        g=(g.dropna(subset=["observer_id"])
           .sort_values(["observer_id","photo_id"],kind="stable")
           .drop_duplicates("observer_id")
           .sort_values("photo_id",kind="stable"))
    valid=(all_coords_ok(g)&g.date.notna().to_numpy()
           &np.isfinite(g[list(AXES)].to_numpy(float)).all(axis=1))
    g=g.loc[valid].reset_index(drop=True)
    n=len(g)
    row["n_photo_climate_complete"]=int(n)
    row["n_original_observers"]=int(g.observer_id.nunique())
    if n<MIN_CLASS or row["n_original_observers"]<MIN_DISTINCT_OBSERVERS:
        row["status"]="insufficient_geoclim_photo_or_independent_observers"
        return row,empty,empty.copy()
    rr=np.column_stack((
        np.floor((90-g.latitude.to_numpy(float))*6),
        np.floor((g.longitude.to_numpy(float)+180)*6),
    )).astype(int)
    ncell=np.unique(rr,axis=0).shape[0]
    row["n_original_raster_cells"]=int(ncell)
    if ncell<MIN_CELLS:
        row["status"]="fewer_than_5_original_climate_raster_cells"
        return row,empty,empty.copy()
    Y=np.stack(g.palette_norm9.to_numpy()).astype(float)
    must(Y.shape==(n,H) and np.isfinite(Y).all()
         and np.allclose(Y.sum(axis=1),1,rtol=0,atol=1e-10),
         "original continuous nine-palette photo measurement invalid")
    labels=pd.Categorical(g.morph,categories=PHOTO_LABELS).codes
    must(np.all(labels>=0),"unknown frozen photo class")
    counts=np.bincount(labels,minlength=4)
    row["n_photo_bins_ge5"]=int(np.sum(counts>=MIN_GROUP))

    env=g[list(AXES)].to_numpy(float)
    means=env.mean(axis=0)
    sd=env.std(axis=0)
    ok=sd>np.maximum(1e-7,2e-10*np.maximum(1,np.abs(means)))
    Z=(env-means)/np.where(ok,sd,1)
    M=compose_macros(Z)
    D,months=design_geo_month(g)
    R=M-D@np.linalg.lstsq(D,M,rcond=1e-8)[0]
    before=(M*M).mean(axis=0)
    target=np.full_like(R,np.nan)
    retained=np.zeros(K)
    for j in range(K):
        other=[k for k in range(K) if k!=j]
        X=R[:,other]
        v=R[:,j]-X@np.linalg.lstsq(X,R[:,j],rcond=1e-9)[0]
        retained[j]=float(np.mean(v*v))
        target[:,j]=v
    retainfrac=retained/np.maximum(before,1e-12)
    eligible=(before>1e-10)&(retainfrac>=CLIMATE_MIN_REMAIN)
    for j,domain in enumerate(DOMAINS):
        row[domain+"_climate_original_variance"]=float(before[j])
        row[domain+"_climate_unique_variance"]=float(retained[j])
        row[domain+"_unique_variance_fraction"]=float(retainfrac[j])
        row[domain+"_identifiable"]=bool(eligible[j])
    row["n_identifiable_climate_domains"]=int(np.sum(eligible))
    row["thermal_original_BIO1_BIO5_corr"]=(
        float(np.corrcoef(Z[:,0],Z[:,2])[0,1])
        if ok[0] and ok[2] else None)
    row["precip_original_BIO12_BIO14_corr"]=(
        float(np.corrcoef(Z[:,3],Z[:,4])[0,1])
        if ok[3] and ok[4] else None)
    if not eligible.any():
        row["status"]="climate_domains_not_identifiable_after_covariate_controls"
        return row,empty,empty.copy()

    basis=fixed_colour_component_basis(Y,labels)
    beta=np.full((K,3,H),np.nan)
    cov=np.full_like(beta,np.nan)
    for j in range(K):
        if not eligible[j]:continue
        x=target[:,j]
        observed=(x@basis/n).reshape(3,H)
        perm=one_month_null(x,basis,months)
        result=observed-perm
        must(np.max(np.abs(result[0]-result[1]-result[2]))<1e-11,
             "source photo label covariance additive identity broken")
        cov[j]=result
        beta[j]=result/retained[j]
    if row["n_photo_bins_ge5"]<2:
        cov[:,1]=np.nan
        beta[:,1]=np.nan
    row["status"]="source_photo_macroclimate_estimable"
    return row,beta,cov


def bootstrap(rows,arrays,seed,rep=B):
    d=pd.DataFrame(rows)
    arr=np.stack(arrays)
    must(arr.shape==(TOTAL_TAXA,K,3,H),"original taxonomic 9-hue tensor shape drift")
    label,genus=pd.factorize(d.genus,sort=True)
    G=len(genus)
    a=np.zeros((G,K,3,H),float)
    n=np.zeros_like(a)
    finite=np.isfinite(arr)
    np.add.at(a,label,np.where(finite,arr,0))
    np.add.at(n,label,finite.astype(float))
    count=n.sum(axis=0)
    mean=np.divide(a.sum(axis=0),count,
                   out=np.full_like(count,np.nan),where=count>0)
    draws=np.random.default_rng(seed).integers(G,size=(rep,G))
    # Sparse enough for deterministic direct genus-bootstrap sufficient stats.
    w=np.zeros((rep,G),float)
    np.add.at(w,(np.repeat(np.arange(rep),G),draws.ravel()),1)
    ba=w@a.reshape(G,-1)
    bn=w@n.reshape(G,-1)
    estimates=np.divide(ba,bn,out=np.full_like(ba,np.nan),where=bn>0)
    return d,arr,mean,count,estimates.reshape((rep,K,3,H)),G


def joint_boot_intervals(point,rep,alpha=.05):
    point=np.asarray(point,float)
    rep=np.asarray(rep,float)
    must(rep.shape==(B,*point.shape),"familywise genus intervals wrong draw shape")
    flat=point.ravel()
    BS=rep.reshape(B,-1)
    finite=np.isfinite(flat)&np.isfinite(BS).all(axis=0)
    idx=np.flatnonzero(finite)
    c=np.full((len(flat),2),np.nan)
    m=np.full_like(c,np.nan)
    threshold=None
    if len(idx):
        sd=BS[:,idx].std(axis=0,ddof=1)
        yes=sd>1e-15
        idx=idx[yes]
        sd=sd[yes]
        if len(idx):
            threshold=float(np.quantile(np.abs((BS[:,idx]-flat[idx])/sd).max(axis=1),1-alpha))
            c[idx,0]=flat[idx]-threshold*sd
            c[idx,1]=flat[idx]+threshold*sd
            m[idx]=np.quantile(BS[:,idx],[alpha/2,1-alpha/2],axis=0).T
    return c.reshape(*point.shape,2),m.reshape(*point.shape,2),{
        "family_size_evaluable":len(idx),
        "maxT_genus_critical":threshold,
        "genus_draws":B,
    }


def alignment(v):
    if not np.isfinite(v).all():return None
    den=float(np.dot(v,v))
    return float((np.dot(v,Q)**2)/den) if den>1e-15 else None


def match_two_domains(df,betas,seed):
    sel=[DOMAINS.index("thermal_level"),DOMAINS.index("precip_amount")]
    x=betas[:,sel,0,:]
    ok=np.isfinite(x).all(axis=(1,2))
    sub=x[ok]
    if len(sub)<40:
        return {"status":"not_enough_2domain_original_source_photo_taxa","n_taxa":len(sub)}
    code,genera=pd.factorize(df.loc[ok,"genus"],sort=True)
    G=len(genera)
    sm=np.zeros((G,2,H),float)
    n=np.bincount(code,minlength=G)
    np.add.at(sm,code,sub)
    mean=sm.sum(axis=0)/n.sum()
    draw=np.random.default_rng(seed).integers(G,size=(B,G))
    w=np.zeros((B,G),float)
    np.add.at(w,(np.repeat(np.arange(B),G),draw.ravel()),1)
    bm=(w@sm.reshape(G,-1)).reshape((B,2,H))
    denom=w@n
    bm/=denom[:,None,None]
    en=np.zeros((B,2))
    for j in range(2):
        en[:,j]=np.einsum("nh,h->n",bm[:,j,:],Q)**2/np.maximum(
                   np.sum(bm[:,j,:]**2,axis=1),1e-16)
    delta=en[:,0]-en[:,1]
    ci=np.quantile(delta,[.025,.975])
    return {
        "status":"same_original_species_temperature_vs_precipitation_macro_domain",
        "n_species":len(sub),
        "n_genera":G,
        "thermal_white_axis_alignment":alignment(mean[0]),
        "precip_white_axis_alignment":alignment(mean[1]),
        "difference_white_axis_alignment":alignment(mean[0])-alignment(mean[1]),
        "paired_genus_bootstrap95":[float(v) for v in ci],
        "each_alignment_genus95":[[float(v) for v in np.quantile(en[:,j],[.025,.975])]
                                   for j in (0,1)],
        "genetic_morph_explained_variance":False,
        "not_a_preregistered_confirmation":True,
    }


def common_photo_bin(rows,arrays,domain,seed):
    d=pd.DataFrame(rows)
    j=DOMAINS.index(domain)
    v=np.stack(arrays)[:,j,:,0]
    good=np.isfinite(v).all(axis=1)
    sub=v[good]
    if len(sub)<50:return {"status":"not_identifiable","n_taxa":len(sub)}
    must(np.max(np.abs(sub[:,0]-sub[:,1]-sub[:,2]))<1e-9,
         "photo original-label components do not add to same taxon total")
    g,gn=pd.factorize(d.loc[good,"genus"],sort=True)
    G=len(gn)
    sums=np.zeros((G,3))
    n=np.bincount(g,minlength=G)
    np.add.at(sums,g,sub)
    draw=np.random.default_rng(seed).integers(G,size=(B,G))
    boots=sums[draw].sum(axis=1)/n[draw].sum(axis=1)[:,None]
    mu=sub.mean(axis=0)
    ratio=float(mu[1]/mu[0]) if abs(mu[0])>1e-9 else None
    shares=boots[:,1]/boots[:,0]
    shares=shares[np.isfinite(shares)]
    return {
        "status":"within_and_between_original_photo_bins_common_taxa",
        "n_taxa":len(sub),
        "n_genera":G,
        "mean_photo_white_partial_slope":{
            COMP[i]:{
                "mean":float(mu[i]),
                "genus_marginal95":[float(v) for v in np.quantile(boots[:,i],[.025,.975])]
            } for i in range(3)},
        "between_photo_bin_fraction_of_signed_mean_if_defined":ratio,
        "genus_bootstrap_ratio95":[float(v) for v in np.quantile(shares,[.025,.975])],
        "not_genetic_morph_frequency":True,
    }


def summarize(rows,arrays,thin):
    d,mat,mean,count,draw,G=bootstrap(rows,arrays,SEED+int(thin)*297)
    projections=np.einsum("kjh,h->kj",mean,Q)
    bproj=np.einsum("bkjh,h->bkj",draw,Q)
    chrom=mean-projections[...,None]*Q
    bchrom=draw-bproj[...,None]*Q
    output={}
    for j,comp in enumerate(COMP):
        wci,wm,wa=joint_boot_intervals(projections[:,j],bproj[:,:,j])
        hci,hm,ha=joint_boot_intervals(mean[:,j,:],draw[:,:,j,:])
        cci,cm,ca=joint_boot_intervals(chrom[:,j,1:],bchrom[:,:,j,1:])
        ds={}
        for k,domain in enumerate(DOMAINS):
            n=int(count[k,j,0])
            if n==0:
                ds[domain]={"status":"not_estimable","n_species":0}
                continue
            x=mat[:,k,j,0]
            x=x[np.isfinite(x)]
            must(len(x)==n,"source nine-hue white fraction denominator mismatch")
            ds[domain]={
                "status":"observational_photo_hue_macroclimate_conditional",
                "n_species":n,
                "n_genera":int(d.loc[np.isfinite(mat[:,k,j,0]),"genus"].nunique()),
                "white_photo_fraction_mean_partial_slope":float(mean[k,j,0]),
                "white_photo_fraction_median_partial_slope":float(np.median(x)),
                "fraction_species_white_slope_positive":float(np.mean(x>0)),
                "q_achromatic_mean_photo_colour_slope":float(projections[k,j]),
                "q_white_genus_simultaneous_five_domains95":[float(v) for v in wci[k]],
                "q_white_family5_excludes_zero":bool(wci[k,0]>0 or wci[k,1]<0),
                "q_white_genus_marginal95":[float(v) for v in wm[k]],
                "white_axis_alignment_of_mean_photo_colour_vector":alignment(mean[k,j]),
                "hue_vector_norm":float(np.linalg.norm(mean[k,j])),
                "eight_chromatic_hue_components":{
                    WCOLOURS[i+1]:{
                        "mean":float(chrom[k,j,i+1]),
                        "genus_simultaneous_five_by_eight95":[float(v) for v in cci[k,i]],
                        "family40_excludes_zero":bool(cci[k,i,0]>0 or cci[k,i,1]<0),
                    } for i in range(8)},
                "full_nine_source_photo_palette":{
                    WCOLOURS[i]:{
                        "mean":float(mean[k,j,i]),
                        "genus_simultaneous_five_by_nine95":[float(v) for v in hci[k,i]],
                        "family45_excludes_zero":bool(hci[k,i,0]>0 or hci[k,i,1]<0),
                    } for i in range(H)},
            }
        output[comp]={
            "white_projection_five_macroclimate_simultaneous_family":wa,
            "chromatic_five_by_eight_macroclimate_simultaneous_family":ca,
            "photo_palette_five_by_nine_simultaneous_family":ha,
            "per_domain":ds,
        }
    cohorts={}
    for cohort,subset in d.groupby("source_cohort",sort=True):
        ids=subset.index.to_numpy(int)
        focused={}
        for macro in ("thermal_level","precip_amount"):
            k=DOMAINS.index(macro)
            x=mat[ids,k,0,0]
            eligible=np.isfinite(x)
            focused[macro]={
                "n_species":int(eligible.sum()),
                "mean_white_photo_fraction_partial_slope":float(x[eligible].mean()) if eligible.any() else None,
                "fraction_species_positive":float(np.mean(x[eligible]>0)) if eligible.any() else None,
            }
        cohorts[str(cohort)]={"n_source_species":len(subset),"domains":focused}
    return {
        "n_original_taxa_in_source":len(d),
        "n_taxa_with_any_domain":int(d.status.eq("source_photo_macroclimate_estimable").sum()),
        "n_taxa_2_photo_bins_ge5":int((d.status.eq("source_photo_macroclimate_estimable")&
                                        d.n_photo_bins_ge5.ge(2)).sum()),
        "n_genus_bootstrap_clusters":G,
        "status_counts":d.status.value_counts().to_dict(),
        "components":output,
        "five_source_photo_cohorts":cohorts,
        "matched_climatic_photo_colour_geometry":match_two_domains(d,mat,SEED+int(thin)*13),
        "same_taxon_photo_coarse_bin_partition":{
            name:common_photo_bin(rows,arrays,name,SEED+int(thin)*17+DOMAINS.index(name))
            for name in ("thermal_level","precip_amount")},
    },d,mat


def synthetic():
    n=40
    rng=np.random.default_rng(SEED)
    s=rng.normal(size=(n,7))
    s[:,2]=s[:,0]+.5*rng.normal(size=n)
    s[:,4]=s[:,3]+.6*rng.normal(size=n)
    m=compose_macros(s)
    must(m.shape==(n,5)
         and np.allclose(m[:,0],(s[:,0]+s[:,2])/2)
         and np.allclose(m[:,1],(s[:,3]+s[:,4])/2),
         "original shared temperature/moisture construction wrong")
    y=rng.dirichlet(np.ones(H),size=n)
    group=np.repeat([0,1],20)
    basis=fixed_colour_component_basis(y,group)
    month=np.repeat([6,7],20)
    x=rng.normal(size=n)
    expected=one_month_null(x,basis,month)
    xp=x.copy()
    # Enumerate exact permutations of four constructed source photos in
    # a small fixture rather than 20! permutations.
    z=np.arange(4,dtype=float)
    yy=np.array([[.8,.2]+[0.]*7,[.7,.3]+[0.]*7,
                 [.2,.8]+[0.]*7,[.1,.9]+[0.]*7])
    bb=fixed_colour_component_basis(yy,np.array([0,0,1,1]))
    mo=np.array([6,6,7,7])
    exact=one_month_null(z,bb,mo)
    import itertools
    perm=[]
    for a in itertools.product([0,1],repeat=2):
        order=[0,1,2,3]
        for i,swap in enumerate(a):
            if swap:order[2*i:2*i+2]=order[2*i:2*i+2][::-1]
        perm.append((z[order]@bb/4).reshape(3,H))
    must(np.allclose(np.mean(perm,axis=0),exact,atol=1e-12),
         "monthly original climate permutation analytic mean not exact")
    v=np.array([8]+[-1]*8,float)
    v-=v.mean()
    must(abs(np.dot(v,Q)**2/(v@v)-1)<1e-12,
         "white-versus-chromatic compositional projection is not unit aligned")
    # Known thermal signal independent of moisture nuisance survives:
    t=rng.normal(size=n)
    p=rng.normal(size=n)
    clim=np.column_stack((t,p,rng.normal(size=(n,3))))
    R=clim-np.mean(clim,axis=0)
    tR=R[:,0]-R[:,1:]@np.linalg.lstsq(R[:,1:],R[:,0],rcond=1e-9)[0]
    must(abs(float(tR@R[:,1]/n))<1e-12,"temperature still collinear with precipitation")
    pigment=.2+.1*t
    must(float(tR@pigment/n)>.001,"true thermal photo source covariance discarded")
    chk=np.stack([rng.dirichlet(np.ones(H)) for _ in range(10)])
    must(np.isfinite(chk).all(),"original source palette data fixture malformed")
    print("FCP_SHARED_CLIMATE_9HUE_THERMAL_MOISTURE_SELFTEST_PASS")


def main():
    p=argparse.ArgumentParser()
    for n in ("photos","worldclim-bio","worldclim-srad","out"):
        p.add_argument("--"+n,type=Path)
    p.add_argument("--self-test",action="store_true")
    args=p.parse_args()
    synthetic()
    if args.self_test:return
    must(args.photos and args.worldclim_bio and args.worldclim_srad and args.out,
         "original photo and both WorldClim source archive inputs required")
    full,photo=load_original(args.photos)
    photo=sample_climate(photo,args.worldclim_bio)
    photo=add_srad(photo,args.worldclim_srad)
    must(len(full)==TOTAL_ROWS and len(photo)==CLASSIFIABLE_ROWS
         and photo.inat_taxon_id.nunique()==TOTAL_TAXA,
         "not exactly original 447500 photos / 227107 classifiable / 4475 taxa")
    args.out.mkdir(parents=True,exist_ok=True)
    panels={}
    for thin in (False,True):
        name="observer_thinned" if thin else "original"
        rows=[];BETA=[];COV=[]
        for ii,(taxid,g) in enumerate(photo.groupby("inat_taxon_id",sort=True)):
            row,beta,cov=one_species(g,thin)
            rows.append(row);BETA.append(beta);COV.append(cov)
            if (ii+1)%1000==0:
                print("FCP_SHARED_THERMAL_MOISTURE_TAXA_PROGRESS",name,ii+1,flush=True)
        panel,df,mat=summarize(rows,BETA,thin)
        must(len(df)==TOTAL_TAXA and df.inat_taxon_id.is_unique
             and int(df.n_original_classifiable.sum())==CLASSIFIABLE_ROWS,
             "incomplete original 4475 taxonomic ledger")
        df.to_csv(args.out/f"fcp4475_{name}_thermal_moisture_support.csv",
                  index=False,float_format="%.12g")
        np.savez_compressed(
            args.out/f"fcp4475_{name}_thermal_moisture_nine_hue_coefficients.npz",
            inat_taxon_id=df.inat_taxon_id.to_numpy(np.int64),
            conditional_slope=mat,
            signed_photo_climate_covariance=np.stack(COV),
            domains=np.array(DOMAINS),
            photo_hues=np.array(WCOLOURS),
            photographed_components=np.array(COMP)
        )
        panels[name]=panel
    obj={
       "schema":"fcp4475_shared_thermal_moisture_source_photo_nine_hue_macroclimate_v1",
       "date_jst":"2026-10-24",
       "source_447500_photo_sha256":SOURCE_SHA,
       "source_WorldClim_BIO_sha256":BIO_ARCHIVE_SHA,
       "source_WorldClim_SRAD_sha256":SRAD_ZIP_SHA,
       "original_photos":TOTAL_ROWS,
       "original_species":TOTAL_TAXA,
       "source_classifiable_photos":CLASSIFIABLE_ROWS,
       "climate_original_seven_axes":list(AXES),
       "fixed_five_climate_domains":list(DOMAINS),
       "thermal_definition":"(original within-species zBIO1 + zBIO5)/2",
       "precip_definition":"(original within-species zBIO12 + zBIO14)/2",
       "climate_domain_adjustment":"thermal-level and precipitation-amount and temperature seasonality BIO4 and precipitation seasonality BIO15 and mean solar, each conditional on the other four; geography linear XYZ and photo calendar month projected first",
       "calendar_month_permutation_exact_expectation":True,
       "temperature_alone_unique_cause_not_claimed":True,
       "source_photo_bin_classification_not_heritable_morph":True,
       "photo_based_colour_not_independent_petal_spectrometry":True,
       "no_verified_botanical_flower_or_camera_background_adaptation_fitness":True,
       "n_whole_genus_resamples":B,
       "species_equal_genus_cluster_inference":True,
       "climate_unique_variance_fraction_minimum":CLIMATE_MIN_REMAIN,
       "source_photo_observer_one_per_taxon_colour_blind":True,
       "not_identifiable_species_retained_as_missing_not_zero":True,
       "bootstrap_families":"independent 5-dimensional achromatic, 5x8 chromatic, 5x9 original source hues for each source component/photo panel; exploratory across families",
       "results":panels,
    }
    (args.out/"shared_thermal_moisture_nine_hue_source_result.json").write_text(
        json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("FCP_SHARED_THERMAL_MOISTURE_FULL_SOURCE_RESULT",json.dumps({
        mode:{
          "n_any":v["n_taxa_with_any_domain"],
          "n_two_bins":v["n_taxa_2_photo_bins_ge5"],
          "thermal":v["components"]["total"]["per_domain"]["thermal_level"],
          "precip":v["components"]["total"]["per_domain"]["precip_amount"],
          "matched":v["matched_climatic_photo_colour_geometry"],
          "between":v["same_taxon_photo_coarse_bin_partition"]["thermal_level"],
          "five_cohorts":v["five_source_photo_cohorts"]
        } for mode,v in panels.items()
    },sort_keys=True,allow_nan=False))

if __name__=="__main__":
    main()
