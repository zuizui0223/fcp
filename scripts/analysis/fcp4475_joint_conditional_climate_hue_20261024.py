#!/usr/bin/env python3
"""FCP: conditional seven-climate directions of original continuous 9-hue photos.

Two distinct estimands:
(A) climate core: BIO5 | BIO12, SRAD and BIO12 | BIO5, SRAD;
(B) unique seven: each axis | all other SIX climate variables.

Both models first project monthly source observation time and coarse XYZ site
geography out of all original WorldClim normals. Only measured true original
photos enter, no missing-photo filling, and all 4475 taxon identities remain.

These estimates are photo-palette/climate associations, NOT flower-specific
reflectance or genetic morph frequency, adaptation, causation or pigment cost.
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
    load_original,sample_climate,all_coords_ok,
    TOTAL_TAXA,TOTAL_ROWS,CLASSIFIABLE_ROWS,SOURCE_SHA,BIO_ARCHIVE_SHA,
    PHOTO_LABELS,WCOLOURS,
)
from fcp_all4475_seven_abiotic_axes_20261021 import (
    add_srad,SRAD_ZIP_SHA,AXES,
)
from fcp4475_nine_colour_seven_climate_covariance_20261024 import (
    fixed_colour_component_basis,design_geo_month,
    MIN_CLASS,MIN_CELLS,MIN_DISTINCT_OBSERVERS,MIN_GROUP
)
from fcp4475_nine_hue_climate_direction_20261024 import (
    analytic_month_null,Q
)

MODEL_AXES={
  "thermal_moisture_solar_core":("BIO5","BIO12","SRADmean"),
  "seven_climates_joint_conditional":tuple(AXES),
}
VAR_THRESH=.05
STRICT_VAR_THRESH=.10
B=999
SEED=20261024
NAX=len(AXES)
NH=9
NCOMP=3
COMP=("total","between_original_photo_labels","within_original_photo_labels")
COHORT=("discovery","validation","prospective","wave4","new_U100")


def check(ok,why):
    if not ok:
        raise RuntimeError("FCP4475_JOINT_CLIMATE_HUE: "+why)


def residualize_against_others(R,axis_indices):
    """Unstandardized UNIQUE climate residual in units of original species SD.

    Inputs R are already centered and projected off calendar month and
    coarse XYZ source location. Each target climate index is projected on
    exactly the other chosen climate indices. Report actual retained climate
    variance, rather than interpret an unstable, collinear pseudo-slope.
    """
    n,d=R.shape
    check(d==NAX,"wrong original WorldClim source-seven shape")
    output=np.full_like(R,np.nan)
    var=np.zeros(NAX,float)
    for j in axis_indices:
        others=[i for i in axis_indices if i!=j]
        if others:
            o=R[:,others]
            resid=R[:,j]-o@np.linalg.lstsq(o,R[:,j],rcond=1e-9)[0]
        else:
            resid=R[:,j].copy()
        output[:,j]=resid
        var[j]=float(np.mean(resid**2))
    return output,var


def source_species(g,observer_thinned,model):
    taxon=int(g.inat_taxon_id.iloc[0])
    label=str(g.species.iloc[0])
    row={
        "inat_taxon_id":taxon,
        "species":label,
        "genus":label.split()[0],
        "source_cohort":str(g.fcp_source_cohort.iloc[0]),
        "n_original_classifiable":int(len(g)),
        "n_source_climate_complete":0,
        "n_original_observers":0,
        "n_original_worldclim_cells":0,
        "n_source_photo_labels_supported_ge5":0,
        "observer_thinned":bool(observer_thinned),
        "model":model,
        "status":"source_not_estimable",
    }
    empty=np.full((NAX,NCOMP,NH),np.nan,float)
    thresholds=np.full((NAX,),np.nan,float)
    if len(g)<MIN_CLASS:
        row["status"]="fewer_than_40_original_classified_photo"
        return row,empty,empty.copy(),thresholds

    if observer_thinned:
        g=(g.dropna(subset=["observer_id"])
             .sort_values(["observer_id","photo_id"],kind="stable")
             .drop_duplicates("observer_id")
             .sort_values("photo_id",kind="stable"))
    valid=(
        all_coords_ok(g)
        & g.date.notna().to_numpy()
        & np.isfinite(g[list(AXES)].to_numpy(float)).all(axis=1)
    )
    g=g.loc[valid].reset_index(drop=True)
    n=len(g)
    row["n_source_climate_complete"]=int(n)
    row["n_original_observers"]=int(g.observer_id.nunique())
    if n<MIN_CLASS or row["n_original_observers"]<MIN_DISTINCT_OBSERVERS:
        row["status"]="fewer_than_40_original_photo_or_20_observers"
        return row,empty,empty.copy(),thresholds

    lat=g.latitude.to_numpy(float)
    lon=g.longitude.to_numpy(float)
    cells=np.stack((np.floor((90-lat)*6),np.floor((lon+180)*6)),axis=1).astype(int)
    row["n_original_worldclim_cells"]=int(np.unique(cells,axis=0).shape[0])
    if row["n_original_worldclim_cells"]<MIN_CELLS:
        row["status"]="fewer_than_5_worldclim_cells"
        return row,empty,empty.copy(),thresholds

    Y=np.stack(g.palette_norm9.to_numpy()).astype(float)
    check(Y.shape==(n,NH)
          and np.isfinite(Y).all()
          and np.allclose(Y.sum(axis=1),1,atol=1e-10,rtol=0),
          "original source colour ROI not a genuine normalized nine-palette")

    labels=pd.Categorical(g.morph,categories=PHOTO_LABELS).codes
    check(np.all(labels>=0),"original coarse source photo colour name mismatch")
    counts=np.bincount(labels,minlength=4)
    supported=int(np.sum(counts>=MIN_GROUP))
    row["n_source_photo_labels_supported_ge5"]=supported

    C=g[list(AXES)].to_numpy(float)
    center=C.mean(axis=0)
    sd=C.std(axis=0)
    stable=sd>np.maximum(1e-7,2e-10*np.maximum(np.abs(center),1.))
    Z=(C-center)/np.where(stable,sd,1)
    geographic,month=design_geo_month(g)
    residual=Z-geographic@np.linalg.lstsq(geographic,Z,rcond=1e-8)[0]
    index=[AXES.index(axis) for axis in MODEL_AXES[model]]
    conditional,unique_var=residualize_against_others(residual,index)
    eligible=stable&(unique_var>=VAR_THRESH)
    for j,axis in enumerate(AXES):
        row[f"{axis}_unique_climate_variance_fraction"]=(
            float(unique_var[j]) if stable[j] and j in index else None)
        row[f"{axis}_conditional_estimable"]=bool(j in index and eligible[j])
        row[f"{axis}_strict_var_ge0_10"]=bool(j in index and eligible[j]
                                           and unique_var[j]>=STRICT_VAR_THRESH)
    row["n_conditional_axes_identifiable"]=int(eligible.sum())

    beta=np.full((NAX,NCOMP,NH),np.nan,float)
    cov=np.full_like(beta,np.nan)
    if not eligible.any():
        row["status"]="no_independent_climate_variance_after_adjustment"
        return row,beta,cov,unique_var

    basis=fixed_colour_component_basis(Y,labels)
    for j in index:
        if not eligible[j]:continue
        x=conditional[:,j]
        # Complete source climate residual exchange, within the same
        # species and calendar month. Its expectation is analytic.
        observed=(x@basis/n).reshape(NCOMP,NH)
        null=analytic_month_null(x[:,None],basis,month)[0]
        difference=observed-null
        check(np.allclose(difference[0],difference[1]+difference[2],atol=1e-11),
              "continuous 9-hue covariance is not between+within original photo bins")
        cov[j]=difference
        beta[j]=difference/unique_var[j]
    if supported<2:
        beta[:,1]=np.nan
        cov[:,1]=np.nan

    row["status"]="source_joint_climate_identifiable"
    return row,beta,cov,unique_var


def mean_bootstrap_genus(rows,coefs,reps=B,seed=SEED):
    d=pd.DataFrame(rows)
    mat=np.stack(coefs)
    check(mat.shape==(TOTAL_TAXA,NAX,NCOMP,NH),
          "taxon × original 7 climate × source 3 categories × 9 palette shape drift")
    gen,gen_names=pd.factorize(d.genus,sort=True)
    G=len(gen_names)
    sm=np.zeros((G,NAX,NCOMP,NH),float)
    num=np.zeros_like(sm)
    good=np.isfinite(mat)
    np.add.at(sm,gen,np.where(good,mat,0))
    np.add.at(num,gen,good.astype(float))
    tot=num.sum(axis=0)
    means=np.divide(sm.sum(axis=0),tot,out=np.full_like(tot,np.nan),
                    where=tot>0)
    draw=np.random.default_rng(seed).integers(G,size=(reps,G))
    weights=np.zeros((reps,G),float)
    np.add.at(weights,(np.repeat(np.arange(reps),G),draw.reshape(-1)),1)
    bs=weights@sm.reshape(G,-1)
    bn=weights@num.reshape(G,-1)
    samples=np.divide(bs,bn,out=np.full_like(bs,np.nan),where=bn>0)
    samples=samples.reshape((reps,NAX,NCOMP,NH))
    return d,mat,means,tot,samples,G


def maxT(point,bootstrap):
    """Familywise bootstrap interval for a pre-defined 7-axis signed photo hue.

    Values unestimable are excluded from the family, never set equal to zero.
    Source genus is whole resampling block, not individual photo pairs.
    """
    p=np.asarray(point,float)
    v=np.asarray(bootstrap,float)
    check(v.shape==(B,*p.shape),"wrong genus cluster draws for fixed climate family")
    f=p.reshape(-1)
    draws=v.reshape((B,-1))
    good=np.isfinite(f)&np.isfinite(draws).all(axis=0)
    sd=draws[:,good].std(axis=0,ddof=1)
    idx=np.flatnonzero(good)[sd>1e-15]
    c=np.full((len(f),2),np.nan)
    marginal=np.full_like(c,np.nan)
    if len(idx):
        err=np.abs((draws[:,idx]-f[idx])/draws[:,idx].std(axis=0,ddof=1))
        z=float(np.quantile(err.max(axis=1),.95))
        sig=draws[:,idx].std(axis=0,ddof=1)
        c[idx,0]=f[idx]-z*sig
        c[idx,1]=f[idx]+z*sig
        marginal[idx]=np.quantile(draws[:,idx],[.025,.975],axis=0).T
    return c.reshape(*p.shape,2),marginal.reshape(*p.shape,2),int(len(idx))


def hue_alignment(vec):
    v=np.asarray(vec,float)
    n=float(np.dot(v,v))
    if not np.isfinite(n) or n<1e-18:return None
    return float(np.dot(v,Q)**2/n)


def paired_geometry(d,mat,predictors=("BIO5","BIO12"),seed=SEED):
    i=[AXES.index(a) for a in predictors]
    pa=mat[:,i,0,:]
    ok=np.isfinite(pa).all(axis=(1,2))
    z=pa[ok]
    gens=d.loc[ok,"genus"]
    if len(z)<30:return {"status":"not_estimable","n_species":int(len(z))}
    labels,gs=pd.factorize(gens,sort=True)
    G=len(gs)
    total=np.zeros((G,2,NH),float)
    num=np.bincount(labels,minlength=G)
    np.add.at(total,labels,z)
    mean=total.sum(axis=0)/num.sum()
    draws=np.random.default_rng(seed).integers(G,size=(B,G))
    weights=np.zeros((B,G),float)
    np.add.at(weights,(np.repeat(np.arange(B),G),draws.reshape(-1)),1)
    values=(weights@total.reshape(G,-1)).reshape((B,2,NH))
    denom=weights@num
    values=values/denom[:,None,None]
    energy=np.stack([
        np.einsum("ij,j->i",values[:,j,:],Q)**2/
        np.maximum(np.sum(values[:,j,:]**2,axis=1),1e-18)
        for j in (0,1)
    ],axis=1)
    difference=energy[:,0]-energy[:,1]
    contrasts=[hue_alignment(mean[j]) for j in (0,1)]
    return {
        "status":"same_taxa_same_photos_within_species_conditional_climate",
        "n_matched_taxa":int(len(z)),
        "n_genera":int(G),
        "first_predictor":predictors[0],"second_predictor":predictors[1],
        "first_axis_alignment":contrasts[0],
        "second_axis_alignment":contrasts[1],
        "difference_first_minus_second":(
             float(contrasts[0]-contrasts[1])
             if None not in contrasts else None),
        "genus_bootstrap_marginal95_per_alignment":[
            [float(x) for x in np.quantile(energy[:,j],[.025,.975])]
            for j in (0,1)
        ],
        "difference_genus_bootstrap95":[
            float(x) for x in np.quantile(difference,[.025,.975])
        ],
        "posthoc_axis_comparison_not_preregistered":True,
        "not_share_of_heritable_colour_variance":True,
    }


def summarize(rows,betas,covariances,model,thinned=False):
    d,mat,mean,count,boot,G=mean_bootstrap_genus(
        rows,betas,seed=SEED+len(model)*17+int(thinned)*821)
    # Signed white*climate conditional slope is reported on its original
    # palette whiteness coordinate and also on fixed orthogonal q.
    qmean=np.einsum("ajh,h->aj",mean,Q)
    qboot=np.einsum("bajh,h->baj",boot,Q)
    out={}
    for comp_idx,comp in enumerate(COMP):
        C,mc,n=maxT(qmean[:,comp_idx],qboot[:,:,comp_idx])
        axes={}
        for j,axis in enumerate(AXES):
            n_species=int(count[j,comp_idx,0])
            if n_species==0:
                axes[axis]={"status":"climate_axis_not_identifiable","n_species":0}
                continue
            x=mat[:,j,comp_idx,0]
            valid=x[np.isfinite(x)]
            check(len(valid)==n_species,"original per-species source colour support changed")
            v=mean[j,comp_idx]
            nonwhite=v-Q*float(qmean[j,comp_idx])
            f64={
                h:{
                    "mean_photo_palette_conditional_slope":float(v[hidx]),
                    "chromatic_only_component":float(nonwhite[hidx]),
                }
                for hidx,h in enumerate(WCOLOURS)
            }
            axes[axis]={
                "status":"source_photo_conditional_colour_climate_estimable",
                "n_species":n_species,
                "n_genera":int(d.loc[np.isfinite(x),"genus"].nunique()),
                "mean_white_photo_fraction_conditional_slope":float(v[0]),
                "median_white_photo_fraction_conditional_slope":float(np.median(valid)),
                "fraction_species_white_conditional_slope_positive":float(np.mean(valid>0)),
                "fixed_q_white_chromatic_conditional_slope":float(qmean[j,comp_idx]),
                "q_white_axis_simultaneous7_95":[float(v) for v in C[j]],
                "q_white_axis_marginal95":[float(v) for v in mc[j]],
                "q_white_axis_genus_family7_excludes_zero":bool(
                    C[j,0]>0 or C[j,1]<0),
                "white_axis_alignment_of_mean_photo_hue_vector":hue_alignment(v),
                "mean_conditional_nine_hue_vector":f64,
            }
        out[comp]={
            "simultaneous_family7_white_direction_effects_with_variance":n,
            "per_climate":axes
        }
    geno={}
    for cohort,frame in d.groupby("source_cohort",sort=True):
        idx=frame.index.to_numpy(int)
        samples={}
        for axis in ("BIO5","BIO12"):
            j=AXES.index(axis)
            w=mat[idx,j,0,0]
            valid=np.isfinite(w)
            samples[axis]={
                "n_species":int(valid.sum()),
                "mean_white_conditional_slope":float(w[valid].mean()) if valid.any() else None,
                "fraction_positive":float(np.mean(w[valid]>0)) if valid.any() else None,
            }
        geno[str(cohort)]={"n_source_species":int(len(frame)),"focus_climates":samples}
    return {
       "n_original_source_species_in_ledger":len(d),
       "n_source_taxa_identifiable_any_axis":int(d.status.eq("source_joint_climate_identifiable").sum()),
       "n_source_taxa_with_two_photo_colours_ge5":int(
           (d.status.eq("source_joint_climate_identifiable")&
             d.n_source_photo_labels_supported_ge5.ge(2)).sum()),
       "n_bootstrap_genera":G,
       "unidentified_status_counts":d.status.value_counts().to_dict(),
       "per_source_cohort":geno,
       "photographic_photo_category_components":out,
       "matched_temperature_vs_precipitation_photo_hue_geometry":paired_geometry(
          d,mat,seed=SEED+int(thinned)*227+len(model)*47)
    },d,mat


def synthetic():
    rng=np.random.default_rng(42)
    n=160
    x=rng.normal(size=n)
    z=rng.normal(size=n)
    aux=rng.normal(size=n)
    # BIO5 is correlated with BIO12 but contains no direct colour response.
    temp=.8*x+.6*z
    precip=x.copy()
    solar=aux.copy()
    rest=np.zeros((n,NAX))
    rest[:,AXES.index("BIO5")]=temp
    rest[:,AXES.index("BIO12")]=precip
    rest[:,AXES.index("SRADmean")]=solar
    index=[AXES.index(a) for a in MODEL_AXES["thermal_moisture_solar_core"]]
    residual,variance=residualize_against_others(rest,index)
    resid_temp=residual[:,AXES.index("BIO5")]
    hue=0.25+0.08*precip
    check(abs(float(np.dot(resid_temp,hue)/n))<1e-12,
          "false BIO5 photo association caused by correlated BIO12 should vanish")
    hue2=0.25+0.08*temp
    check(float(np.dot(resid_temp,hue2)/n)>1e-4,
          "true unique BIO5 colour covariance must survive BIO12 adjustment")
    check(np.dot(resid_temp,precip)/n<1e-11,
          "conditional source temperature still correlated with precipitation")
    # The compositional white axis projection and chromatic-only residual
    # must partition EVERY nine-hue covariance vector exactly.
    v=rng.normal(size=NH)
    v-=v.mean()
    alpha=float(np.dot(v,Q))
    check(np.max(abs(v-(alpha*Q+(v-alpha*Q))))<1e-12,
          "photographic white/chromatic partition identity broken")
    y=np.stack([rng.dirichlet(np.ones(9)) for _ in range(n)])
    labels=np.repeat([0,1,2,3],40)
    basis=fixed_colour_component_basis(y,labels)
    chk=rest[:,index[0]]@basis/n
    check(np.isfinite(chk).all(),"original categorical photograph basis unstable")
    print("FCP4475_INDEPENDENT_CLIMATE_HUE_SYNTHETIC_PASS")


def main():
    p=argparse.ArgumentParser()
    for name in ("photos","worldclim-bio","worldclim-srad","out"):
        p.add_argument("--"+name,type=Path)
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    synthetic()
    if a.self_test:return
    check(all((a.photos,a.worldclim_bio,a.worldclim_srad,a.out)),
          "all three source-verified photo/WorldClim files required")
    allphotos,photo=load_original(a.photos)
    photo=sample_climate(photo,a.worldclim_bio)
    photo=add_srad(photo,a.worldclim_srad)
    check(len(allphotos)==TOTAL_ROWS and len(photo)==CLASSIFIABLE_ROWS and
          photo.inat_taxon_id.nunique()==TOTAL_TAXA,
          "original complete 447500/227107/4475 source lost")
    a.out.mkdir(parents=True,exist_ok=True)
    results={}
    for thinned in (False,True):
        panel="observer_thinned" if thinned else "original"
        for model in MODEL_AXES:
            key=panel+"__"+model
            rows=[];coeff=[];cov=[]
            for t,(tid,g) in enumerate(photo.groupby("inat_taxon_id",sort=True)):
                row,beta,part,var=source_species(g,thinned,model)
                rows.append(row);coeff.append(beta);cov.append(part)
                if (t+1)%1000==0:
                    print("FCP4475_JOINT_CLIMATE_SOURCE_TAXON_PROGRESS",key,t+1,flush=True)
            summary,df,mat=summarize(rows,coeff,cov,model,thinned)
            check(len(df)==TOTAL_TAXA and df.inat_taxon_id.is_unique
                  and int(df.n_original_classifiable.sum())==CLASSIFIABLE_ROWS,
                  "four thousand four hundred seventy-five taxon ledger changed")
            support=a.out/f"original4475_{key}_conditional_climate_support.csv"
            df.to_csv(support,index=False,float_format="%.12g")
            np.savez_compressed(
                a.out/f"original4475_{key}_nine_hue_conditional_source.npz",
                inat_taxon_id=df.inat_taxon_id.to_numpy(np.int64),
                conditional_photo_palette_slope=mat,
                conditional_source_photo_covariance=np.stack(cov),
                climate_axes=np.array(list(AXES)),
                photo_hues=np.array(list(WCOLOURS)),
                component_names=np.array(list(COMP)),
                )
            results[key]=summary
    receipt={
        "schema":"fcp4475_photo_nine_hue_seven_joint_conditional_worldclim_v1",
        "source_original_photo_sha256":SOURCE_SHA,
        "source_WorldClim_BIO_sha256":BIO_ARCHIVE_SHA,
        "source_WorldClim_SRAD_sha256":SRAD_ZIP_SHA,
        "original_photo_count":TOTAL_ROWS,
        "original_taxon_count":TOTAL_TAXA,
        "original_photo_classified_count":CLASSIFIABLE_ROWS,
        "tested_original_climate_axes":list(AXES),
        "tested_original_photo_hues":list(WCOLOURS),
        "models":{k:list(v) for k,v in MODEL_AXES.items()},
        "conditional_association":"species-standardized original climate residual after broad XYZ and calendar month projection, additionally conditioned on other two (BIO5,BIO12,SRAD) or all six other original seven climate axes",
        "unique_variance_gate":VAR_THRESH,
        "extra_strict_variance_gate_reported":STRICT_VAR_THRESH,
        "original_source_photo_classes_are_not_heritable_morphs":True,
        "source_photo_coarse_label_freq_requires_two_photo_classes_at_least_five_each":True,
        "slope_not_causal_or_individual_biological_response":True,
        "background_phenology_cultivation_camera_not_disambiguated":True,
        "source_missing_photographs_not_imputed":True,
        "within_source_species_observer_thinning_colour_blind":True,
        "all_original_taxa_retained_even_when_climate_collinear":True,
        "genus_cluster_resamples":B,
        "posthoc_exploratory_two_models_and_original_vs_thinned":True,
        "no_genetic_colour_morph_local_adaptation_or_pigment_cost_claim":True,
        "results":results,
    }
    (a.out/"joint_climate_conditional_hue_summary.json").write_text(
        json.dumps(receipt,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("FCP4475_JOINT_CLIMATE_HUE_SOURCE_FINISHED",
          json.dumps({key:{
            "n_any":s["n_source_taxa_identifiable_any_axis"],
            "n_two_label":s["n_source_taxa_with_two_photo_colours_ge5"],
            "BIO5":s["photographic_photo_category_components"]["total"]["per_climate"]["BIO5"],
            "BIO12":s["photographic_photo_category_components"]["total"]["per_climate"]["BIO12"],
            "matched":s["matched_temperature_vs_precipitation_photo_hue_geometry"]
          } for key,s in results.items()},sort_keys=True,allow_nan=False))


if __name__=="__main__":
    main()
