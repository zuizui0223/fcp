#!/usr/bin/env python3
"""FCP 447500 source photographs: nine hue direction across seven climates.

All source-photo climate associations are observational. No source genetic
morph, flower-organ, pigment-cost, local adaptation, or fitness claims.

This follow-up tests whether an observed association of photographed white
fraction with BIO5 reflects a reproducible WHITE vs CHROMATIC direction,
or redistribution among different chromatic hues. This exact orthogonal
geometry is distinct from original categorical colour-bin frequencies.

Per-species multivariate climate-colour COVARIANCE, not a causal slope:
 C_i = Cov_i(Z_resid, colour vector) -
       E[Cov_i(Z_resid permuted within calendar month, colour vector)].
The expectation is computed ANALYTICALLY, avoiding Monte Carlo error in a
signed linear covariance estimand. The term "colour vector" is a source
photograph ROI nine-palette mixture, not spectrally measured floral petals.

Bootstrap whole genera; 7-axis simultaneous family for the white fraction and
7x8-family (56) for colour redistribution INDEPENDENT of the white fraction.
The families are exploratory and not jointly protected across 3 components
and both original/observer-thinned runs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"/"analysis"))
from fcp_full4475_hypothesis_census_20261021 import (
    TOTAL_TAXA,TOTAL_ROWS,CLASSIFIABLE_ROWS,SOURCE_SHA,BIO_ARCHIVE_SHA,
    load_original,sample_climate,all_coords_ok,PHOTO_LABELS,WCOLOURS
)
from fcp_all4475_seven_abiotic_axes_20261021 import (
    add_srad,SRAD_ZIP_SHA,AXES
)
from fcp4475_nine_colour_seven_climate_covariance_20261024 import (
    design_geo_month,fixed_colour_component_basis,
    MIN_CLASS,MIN_CELLS,MIN_DISTINCT_OBSERVERS,MIN_GROUP,
    RESIDUAL_CLIMATE_VARIANCE_MIN
)

SEED=20261024
B=1999
NAX=len(AXES)
NC=len(WCOLOURS)
COMP=("total","between_photo_labels","within_photo_labels")
Q=np.array([8.]+[-1.]*8)/np.sqrt(72.)  # fixed white-versus-other palette direction


def check(v,message):
    if not v:
        raise RuntimeError("FCP4475_COLOUR_CLIMATE_DIRECTION: "+message)


def original_sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for data in iter(lambda:f.read(1<<20),b""):h.update(data)
    return h.hexdigest()


def analytic_month_null(climate,source_basis,month):
    """Exact conditional permutation EXPECTATION under within-calendar-month shuffle.

    climate: n x 7 source photo-associated climate residual; basis: n x 27
    original 9 hue total/between/within colour basis; month: n calendar label.
    Shuffling the entire seven-climate vector within its calendar month group
    preserves cross-climate associations in the permutation opportunity.
    """
    climate=np.asarray(climate,float)
    source_basis=np.asarray(source_basis,float)
    n=len(month)
    check(climate.shape==(n,NAX) and source_basis.shape==(n,3*NC),
          "analytic original-photo month-null source shape mismatch")
    null=np.zeros((NAX,3*NC),float)
    for m in np.unique(month):
        idx=np.flatnonzero(month==m)
        null+=np.outer(climate[idx].sum(axis=0),source_basis[idx].mean(axis=0))/n
    return null.reshape(NAX,3,NC)


def one_species(source,observer_thinned=False):
    s=str(source.species.iloc[0])
    row={
        "inat_taxon_id":int(source.inat_taxon_id.iloc[0]),
        "species":s,
        "genus":s.split()[0],
        "source_cohort":str(source.fcp_source_cohort.iloc[0]),
        "n_original_classified":int(len(source)),
        "one_photo_per_observer":bool(observer_thinned),
        "n_complete_photo":0,
        "n_original_observers":0,
        "n_10arcminute_worldclim_cells":0,
        "n_photo_classes_supported_5each":0,
        "status":"source_not_estimable",
    }
    empty=np.full((NAX,3,NC),np.nan,float)
    if len(source)<MIN_CLASS:
        row["status"]="less_than_40_classified_source_photos"
        return row,empty

    g=source
    if observer_thinned:
        g=(g.dropna(subset=["observer_id"])
             .sort_values(["observer_id","photo_id"],kind="stable")
             .drop_duplicates("observer_id")
             .sort_values("photo_id",kind="stable"))
    valid=(all_coords_ok(g)&g.date.notna().to_numpy()
           &np.isfinite(g[list(AXES)].to_numpy(float)).all(axis=1))
    g=g.loc[valid].reset_index(drop=True)
    n=len(g)
    row["n_complete_photo"]=int(n)
    row["n_original_observers"]=int(g.observer_id.nunique())
    if n<MIN_CLASS or row["n_original_observers"]<MIN_DISTINCT_OBSERVERS:
        row["status"]="insufficient_geovalid_classified_photo_or_observer"
        return row,empty

    coords=np.column_stack((
        np.floor((90-g.latitude.to_numpy(float))*6),
        np.floor((g.longitude.to_numpy(float)+180)*6)
    )).astype(np.int64)
    cells=int(np.unique(coords,axis=0).shape[0])
    row["n_10arcminute_worldclim_cells"]=cells
    if cells<MIN_CELLS:
        row["status"]="insufficient_distinct_climate_raster_cells"
        return row,empty

    Y=np.stack(g.palette_norm9.to_numpy()).astype(float)
    check(Y.shape==(n,NC) and np.isfinite(Y).all() and
          np.allclose(Y.sum(axis=1),1,atol=1e-10,rtol=0),
          "original nine-colour palette differs from original source")
    labels=pd.Categorical(g.morph,categories=PHOTO_LABELS).codes
    check(np.all(labels>=0),"unrecognized pre-existing source photo class")
    colour_counts=np.bincount(labels,minlength=4)
    row["n_photo_classes_supported_5each"]=int(np.count_nonzero(colour_counts>=MIN_GROUP))
    class_between_eligible=row["n_photo_classes_supported_5each"]>=2

    Z=g[list(AXES)].to_numpy(float)
    average=Z.mean(axis=0)
    scale=Z.std(axis=0)
    varying=(scale>np.maximum(1e-7,2e-10*np.maximum(np.abs(average),1)))
    Z=(Z-average)/np.where(varying,scale,1)
    D,month=design_geo_month(g)
    R=Z-D@np.linalg.lstsq(D,Z,rcond=1e-8)[0]
    residual_frac=np.mean(R**2,axis=0)/np.maximum(np.mean(Z**2,axis=0),1e-12)
    eligible=varying&(residual_frac>=RESIDUAL_CLIMATE_VARIANCE_MIN)
    row["n_climate_axes_identifiable"]=int(np.count_nonzero(eligible))
    row["BIO5_climate_residual_variance_fraction"]=(
       float(residual_frac[2]) if varying[2] else None
    )
    if not eligible.any():
        row["status"]="geography_month_absorbs_climate_or_climate_constant"
        return row,empty

    basis=fixed_colour_component_basis(Y,labels)
    observed=(R.T@basis/n).reshape(NAX,3,NC)
    null=analytic_month_null(R,basis,month)
    delta=observed-null
    check(np.max(np.abs(delta[:,0]-delta[:,1]-delta[:,2]))<1e-11,
          "total white/palette covariance does not equal components")
    delta[~eligible]=np.nan
    if not class_between_eligible:
        delta[:,1]=np.nan
    row["status"]="photo_climate_direction_estimable"
    row["between_coarse_photo_bins_estimable"]=bool(class_between_eligible)

    return row,delta


def mean_cluster_bootstrap(rows,arrays,seed=SEED,reps=B):
    """Species-equal mean; independent genus cluster re-draw for uncertainty."""
    d=pd.DataFrame(rows)
    matrix=np.stack(arrays).astype(float)
    check(matrix.shape==(TOTAL_TAXA,NAX,3,NC),
          "original 4475 taxon x seven climate x three components x nine hue shape invalid")
    geno,unique=pd.factorize(d.genus,sort=True)
    G=len(unique)
    sums=np.zeros((G,NAX,3,NC),float)
    counts=np.zeros_like(sums)
    present=np.isfinite(matrix)
    np.add.at(sums,geno,np.where(present,matrix,0))
    np.add.at(counts,geno,present.astype(float))
    denominator=counts.sum(axis=0)
    obs=np.divide(sums.sum(axis=0),denominator,
                  out=np.full((NAX,3,NC),np.nan),where=denominator>0)
    draw=np.random.default_rng(seed).integers(G,size=(reps,G))
    # Count-resampling matrix uses much less memory than the huge broadcast
    # of (reps x G x 7 x 3 x 9) source-data tensors.
    weights=np.zeros((reps,G),float)
    np.add.at(weights,(np.repeat(np.arange(reps),G),draw.reshape(-1)),1)
    numer=(weights@sums.reshape(G,-1))
    denom=(weights@counts.reshape(G,-1))
    samples=np.divide(numer,denom,out=np.full_like(numer,np.nan),where=denom>0)
    samples=samples.reshape(reps,NAX,3,NC)
    return d,matrix,obs,denominator,samples,G


def simultaneous_family(point,replicates):
    """point shape (...), replicates shape (B,...) for one chosen family."""
    p=np.asarray(point,float)
    boot=np.asarray(replicates,float)
    check(boot.shape==(B,*p.shape),"source bootstrap family shape mismatch")
    flat=p.reshape(-1)
    each=boot.reshape(B,-1)
    finite=(np.isfinite(flat)&np.isfinite(each).all(axis=0))
    check(finite.any(),"no effect size available for genus family")
    sd=each[:,finite].std(axis=0,ddof=1)
    good=sd>1e-14
    keep=np.flatnonzero(finite)[good]
    ci=np.full((len(flat),2),np.nan,float)
    marginal=np.full_like(ci,np.nan)
    family_n=0
    cutoff=None
    if len(keep):
        margin=np.abs((each[:,keep]-flat[keep])/sd[good]).max(axis=1)
        cutoff=float(np.quantile(margin,.95))
        ci[keep,0]=flat[keep]-cutoff*sd[good]
        ci[keep,1]=flat[keep]+cutoff*sd[good]
        marginal[keep]=np.quantile(each[:,keep],[.025,.975],axis=0).T
        family_n=len(keep)
    return ci.reshape(*p.shape,2),marginal.reshape(*p.shape,2),{
        "family_effects_with_variance":family_n,"max_t_quantile_95":cutoff,
        "bootstrap_cluster_resamples":B,
    }


def project_achromatic_chromatic(matrix):
    """For compositional vectors whose nine hues sum to zero, exact split."""
    w=np.sum(matrix*Q,axis=-1)  # +/- 1.06066 times WHITE component
    chrom=matrix-w[...,None]*Q
    return w,chrom


def summarise_panel(rows,arrays,observer_thinned=False):
    d,mat,mean,count,boots,genus_n=mean_cluster_bootstrap(
        rows,arrays,seed=SEED+int(observer_thinned)*517)
    check(len(d)==TOTAL_TAXA and int(d.n_original_classified.sum())==CLASSIFIABLE_ROWS,
          "source taxonomy or photo records no longer count to original")
    w,chrom=project_achromatic_chromatic(mean)
    w_boot,ch_boot=project_achromatic_chromatic(boots)
    # Main family: signed white component across 7 climatic axes.
    # Separately compare eight chromatic hue directions with no white component
    # under a 7-axis x 8-colour family.
    results={}
    for j,cname in enumerate(COMP):
        wCI,wMarg,wInfo=simultaneous_family(w[:,j],w_boot[:,:,j])
        hCI,hMarg,hInfo=simultaneous_family(mean[:,j,:],boots[:,:,j,:])
        cCI,cMarg,cInfo=simultaneous_family(chrom[:,j,1:],ch_boot[:,:,j,1:])
        climates={}
        for a,axis in enumerate(AXES):
            taxa=int(count[a,j,0])
            vals=mat[:,a,j,0]
            valid=vals[np.isfinite(vals)]
            if taxa==0:
                climates[axis]={"status":"not_estimable"}
                continue
            check(len(valid)==taxa,"missing whiteness denominator mismatch")
            whiteMean=float(mean[a,j,0])
            qv=float(w[a,j])
            chromVector=chrom[a,j,1:]
            hueNorm=float(np.linalg.norm(mean[a,j]))
            white_fraction=(float(qv*qv/(hueNorm*hueNorm))
                            if hueNorm>1e-12 else None)
            reps=list()
            for h,colour in enumerate(WCOLOURS):
                ch={
                    "estimated_signed_photo_palette_covariance":float(mean[a,j,h]),
                    "genus_bootstrap_marginal95":[float(v) for v in hMarg[a,h]],
                    "genus_bootstrap_simultaneous_63_95":[float(v) for v in hCI[a,h]],
                    "simultaneous63_excludes_zero":bool(
                        hCI[a,h,0]>0 or hCI[a,h,1]<0),
                }
                if h>0:
                    ch["chromatic_only_component"]=float(chrom[a,j,h])
                    ch["chromatic_only_genus_simultaneous_56_95"]=[
                        float(v) for v in cCI[a,h-1]]
                    ch["chromatic_56_excludes_zero"]=bool(
                        cCI[a,h-1,0]>0 or cCI[a,h-1,1]<0)
                reps.append((colour,ch))
            climates[axis]={
                "status":"photographic_original_source_climate_geometry",
                "n_species":taxa,
                "n_genera":int(d.loc[np.isfinite(mat[:,a,j,0]),"genus"].nunique()),
                "white_positive_fraction_species":float(np.mean(valid>0)),
                "white_negative_fraction_species":float(np.mean(valid<0)),
                "white_median_species_covariance":float(np.median(valid)),
                "white_mean_species_covariance":whiteMean,
                "white_genus_marginal95":[float(v) for v in hMarg[a,0]],
                "white_genus_simultaneous_63_95":[float(v) for v in hCI[a,0]],
                "white_family_63_excludes_zero":bool(hCI[a,0,0]>0 or hCI[a,0,1]<0),
                "fixed_achromatic_white_projection":qv,
                "white_projection_genus_simultaneous_7_95":[float(v) for v in wCI[a]],
                "white_axis_7_excludes_zero":bool(wCI[a,0]>0 or wCI[a,1]<0),
                "photo_colour_mean_vector_norm":hueNorm,
                "white_axis_energy_fraction_of_mean_vector":white_fraction,
                "chromatic_only_mean_vector_norm":float(np.linalg.norm(chromVector)),
                "chromatic_only_56_axis_any_excludes_zero":bool(np.any(
                    (cCI[a,:,0]>0)|(cCI[a,:,1]<0))),
                "photographed_nine_hue_components":dict(reps),
            }
        results[cname]={
            "white_projected_simultaneous_seven_axis_family":wInfo,
            "full_photo_palette_simultaneous_63_family":hInfo,
            "chromatic_only_seven_by_eight_simultaneous_56_family":cInfo,
            "axes":climates,
        }
    by_cohort={}
    for cohort,sub in d.groupby("source_cohort",sort=True):
        idx=sub.index.to_numpy(int)
        out={}
        for a,axis in enumerate(AXES):
            vals=mat[idx,a,0,0]
            vals=vals[np.isfinite(vals)]
            out[axis]={
                "n_eligible":len(vals),
                "mean_white_signed_covariance":float(vals.mean()) if len(vals) else None,
                "white_fraction_positive":float(np.mean(vals>0)) if len(vals) else None
            }
        by_cohort[str(cohort)]={
            "n_taxa_in_source":len(sub),
            "by_climate":out,
        }
    return {
        "n_taxa_in_original_roster":len(d),
        "n_taxa_with_any_identifiable_seven_climate":int(
            d.status.eq("photo_climate_direction_estimable").sum()),
        "n_taxa_with_two_coarse_photo_categories_5_each":int(
            (d.status.eq("photo_climate_direction_estimable")&
             d.n_photo_classes_supported_5each.ge(2)).sum()),
        "n_genera_bootstrap":genus_n,
        "failures_and_estimation_status_counts":d.status.value_counts().to_dict(),
        "source_cohort_sign_replication":by_cohort,
        "component_results":results,
    },d,mat


def main():
    p=argparse.ArgumentParser()
    for k in ("photos","worldclim-bio","worldclim-srad","out"):
        p.add_argument("--"+k,type=Path)
    p.add_argument("--self-test",action="store_true")
    opts=p.parse_args()
    if opts.self_test:
        rng=np.random.default_rng(2)
        n=6
        Y=np.stack([np.full(9,0.) for i in range(n)])
        Y[:,0]=np.array([.3,.4,.2,.5,.1,.3])
        Y[:,1]=1-Y[:,0]
        C=np.array([0,0,1,1,2,2])
        basis=fixed_colour_component_basis(Y,C)
        climate=rng.normal(size=(n,NAX))
        month=np.array([4,4,5,5,6,6])
        original=analytic_month_null(climate,basis,month)
        # Exact permutation expectation over 2x2x2=8 possible assignments.
        import itertools
        all_cov=[]
        for sel in itertools.product([0,1],repeat=3):
            idx=np.arange(n)
            for j,t in enumerate(sel):
                if t:idx[2*j:2*j+2]=idx[2*j:2*j+2][::-1]
            all_cov.append((climate[idx].T@basis/n).reshape(NAX,3,NC))
        check(np.allclose(np.mean(all_cov,axis=0),original,atol=1e-13),
              "exact source month-permutation expectation differs from exhaustive worlds")
        z=rng.normal(size=(NAX,3,NC))
        z-=z.mean(axis=-1,keepdims=True)
        white,chroma=project_achromatic_chromatic(z)
        check(np.max(np.abs(white[...,None]*Q+chroma-z))<1e-12,
              "white vs chromatic orthogonal source geometry identity broken")
        check(np.max(np.abs(chroma[...,0]))<1e-12,
              "chromatic hue contrast incorrectly retains white coordinate")
        print("FCP4475_EXACT_PHOTO_MONTH_NULL_HUE_DIRECTION_SELFTEST_PASS")
        return

    check(opts.photos and opts.worldclim_bio and opts.worldclim_srad and opts.out,
          "source original photo, BIO and monthly solar archives required")
    full,photo=load_original(opts.photos)
    photo=sample_climate(photo,opts.worldclim_bio)
    photo=add_srad(photo,opts.worldclim_srad)
    check(len(full)==TOTAL_ROWS and len(photo)==CLASSIFIABLE_ROWS and
          photo.inat_taxon_id.nunique()==TOTAL_TAXA,
          "original 447500/227107/4475 original-photo source missing")
    opts.out.mkdir(parents=True,exist_ok=True)
    results={}
    for thin in (False,True):
        panel="observer_thinned" if thin else "original"
        rows=[];arrays=[]
        for j,(tid,g) in enumerate(photo.groupby("inat_taxon_id",sort=True)):
            row,matrix=one_species(g,thin)
            rows.append(row);arrays.append(matrix)
            if (j+1)%1000==0:
                print("FCP4475_CONTINUOUS_HUE_CLIMATE_SOURCE_TAXA",panel,j+1,flush=True)
        summary,df,mat=summarise_panel(rows,arrays,thin)
        results[panel]=summary
        df.to_csv(opts.out/f"original4475_{panel}_photo_climate_colour_direction_denominators.csv",index=False)
        # Exact nine-photo-hue × seven-climate per taxon coefficients, taxonomy-only;
        # no original photo IDs, observer IDs or exact coordinates published.
        coefficients=pd.DataFrame({
            "inat_taxon_id":df.inat_taxon_id,
            "species":df.species,
            "source_cohort":df.source_cohort,
        })
        for a,axis in enumerate(AXES):
            for j,comp in enumerate(COMP):
                for k,h in enumerate(WCOLOURS):
                    coefficients[f"{axis}_{comp}_{h}_signed_covariance_month_null_excess"]=mat[:,a,j,k]
        coefficients.to_csv(opts.out/f"original4475_{panel}_per_taxon_seven_climate_nine_hue.csv",index=False,float_format="%.12g")
    data={
        "schema":"fcp4475_original_photo_nine_hue_seven_climate_white_vs_chromatic_covariance_v1",
        "date_jst":"2026-10-24",
        "original_photo_sha256":SOURCE_SHA,
        "worldclim_BIO_archive_sha256":BIO_ARCHIVE_SHA,
        "worldclim_SRAD_archive_sha256":SRAD_ZIP_SHA,
        "source_n_original_photos":TOTAL_ROWS,
        "source_n_original_taxa":TOTAL_TAXA,
        "source_n_classified_photo":CLASSIFIABLE_ROWS,
        "source_nine_hues":list(WCOLOURS),
        "climate_axes":list(AXES),
        "covariance_nomenclature":"species-standardised observed climate minus global XYZ-and-photo-month projected trend; covariance with nine-hue original photo ROI; exact analytic same-month permutation expected value removed",
        "analytic_null_rule":"source climate/photo vectors permuted independently within each source species x calendar month; expectation evaluated exactly, not p value",
        "white_achromatic_axis":"fixed q=(8,-1,-1,-1,-1,-1,-1,-1,-1)/sqrt(72); relative photographed whiteness, not pigment synthesis pathway",
        "chromatic_only_subspace":"remove orthogonal q projection; nine photographic fractions still compositional",
        "biological_scope":"photographic source measurements of within-named-species variation; no independent botanist perianth validation, genetic inheritance, fitness or adaptation",
        "other_climates_correlated":True,
        "five_cohorts_same_photo_provider_and_assay":True,
        "bootstrap_genus_clusters":B,
        "species_equal_inference":True,
        "hypothesis_testing":"genus-bootstrap simultaneous max-t seven climate axes for white direction; separate family of seven by eight chromatic-only nonwhite hues and 7x9 original hue components; multiple component/source panels exploratory",
        "no_causal_climate_or_pigment_cost_claim":True,
        "all_4475_species_retained_even_if_not_estimable":True,
        "results":results,
    }
    (opts.out/"seven_climate_nine_hue_source_signatures.json").write_text(
        json.dumps(data,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("FCP4475_NINE_HUE_CLIMATE_DIRECTION_COMPLETE",json.dumps({
        p:{
          "eligible":v["n_taxa_with_any_identifiable_seven_climate"],
          "n_genus":v["n_genera_bootstrap"],
          "n_2classes":v["n_taxa_with_two_coarse_photo_categories_5_each"],
          "BIO5":v["component_results"]["total"]["axes"]["BIO5"],
          "five_cohorts":{k:co["by_climate"]["BIO5"] for k,co in v["source_cohort_sign_replication"].items()},
        } for p,v in results.items()},allow_nan=False,sort_keys=True,default=str))

if __name__ == "__main__":
    main()
