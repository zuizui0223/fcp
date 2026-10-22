#!/usr/bin/env python3
"""FCP global 4,475-taxon harmonised measurement and diagnostic census.

Every one of the 4,475 genuine source taxa is an analysis ROW, not just an
eligibility gate. True measurement denominators, NAs and scope restrictions
are part of the result. This script DOES NOT invent historic population
genotypes, colour morph inheritance, fitness, or genotype-specific chemistry.

H1: common nine-palette white-axis descriptive W for all feasible species.
H2: 50-km 4-state photographic mismatch, conditional 199-vertex permutations
    for EVERY species with any valid local photo pair, and nested supports.
H3: seven? Only SIX source-locked bioclim normal location covariates; mean
    white-minus-other BIO5/BIO12 contrasts plus pairwise phenotypic distance.
H4: independent-observer <=5-km same-year-month contrasting PHOTO-label pairs;
    local genetic maintenance status never upgraded.
H5: join existing source-specific literature coverage to all 4475 taxonomy IDs.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
from pathlib import Path
import sys
import zipfile

import numpy as np
import pandas as pd
from PIL import Image

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"packages"/"disttrait"/"src"))
from disttrait.geometry import one_vs_rest_contrast, two_mode_axis

SOURCE_SHA="8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6"
BIO_ARCHIVE_SHA="00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf"
COVERAGE_BLOB="e027179b818c5ae81779c798509c4cbb09579cd6"
TOTAL_ROWS=447500
TOTAL_TAXA=4475
CLASSIFIABLE_ROWS=227107
BIONUM=(1,4,5,12,14,15)
WCOLOURS=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
PHOTO_LABELS=("white","yellow_orange","red_pink","blue_purple")
FRACTIONS=tuple("flower_fraction_"+v for v in WCOLOURS)
COLUMNS=("inat_taxon_id","species","fcp_source_cohort","photo_id","observation_id",
         "observer_id","latitude","longitude","positional_accuracy_m","observed_on",
         "global_classifiable","morph",*FRACTIONS)
PERMS=199
SEED=20261021
Q=one_vs_rest_contrast(9,focal_index=0)


def check(ok,reason):
    if not ok:
        raise RuntimeError("FCP_ALL_4475_SOURCE_OR_ESTIMAND_ERROR: "+reason)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):h.update(b)
    return h.hexdigest()


def load_original(path):
    check(sha(path)==SOURCE_SHA,"original authentic photo CSV SHA mismatch")
    t=pd.read_csv(path,usecols=list(COLUMNS),low_memory=False)
    check(len(t)==TOTAL_ROWS and t.inat_taxon_id.nunique()==TOTAL_TAXA,
          "the source has not retained all 4475 species x 100 source rows")
    check(t.photo_id.nunique()==TOTAL_ROWS,"photo ID not unique")
    check(t.groupby("inat_taxon_id").size().eq(100).all(),"one species lacks 100 original source photographs")
    check(t.groupby("inat_taxon_id").fcp_source_cohort.nunique().eq(1).all(),
          "source cohort is not species-disjoint")
    flag=t.global_classifiable.fillna(False).astype(str).str.strip().str.lower().isin(("true","1","yes"))
    check(flag.sum()==CLASSIFIABLE_ROWS,"frozen four-state classified-photo denominator drift")
    check(t.loc[flag,"morph"].isin(PHOTO_LABELS).all(),"unknown classifiable morph")
    d=t.loc[flag].copy()
    d["latitude"]=pd.to_numeric(d.latitude,errors="coerce")
    d["longitude"]=pd.to_numeric(d.longitude,errors="coerce")
    d["positional_accuracy_m"]=pd.to_numeric(d.positional_accuracy_m,errors="coerce")
    d["date"]=pd.to_datetime(d.observed_on,errors="coerce")
    p=d[list(FRACTIONS)].to_numpy(dtype=float)
    check(np.isfinite(p).all() and np.all(p>=0),"bad nine-palette flower measurements")
    total=p.sum(axis=1)
    check(np.all(total>0),"no usable palette mass in supposedly classified photo")
    p=p/total[:,None]
    d["palette_norm9"]=[a for a in p]
    return t,d


def sample_climate(d,archive):
    check(sha(archive)==BIO_ARCHIVE_SHA,"frozen WorldClim bioclim zip SHA mismatch")
    la=d.latitude.to_numpy(float)
    lo=d.longitude.to_numpy(float)
    valid=np.isfinite(la)&np.isfinite(lo)&(la>=-90)&(la<90)&(lo>=-180)&(lo<180)
    rr=np.zeros(len(d),dtype=int);cc=np.zeros(len(d),dtype=int)
    rr[valid]=np.floor((90-la[valid])*6).astype(int)
    cc[valid]=np.floor((lo[valid]+180)*6).astype(int)
    valid&=(rr>=0)&(rr<1080)&(cc>=0)&(cc<2160)
    with zipfile.ZipFile(archive) as z:
        for bi in BIONUM:
            name=f"wc2.1_10m_bio_{bi}.tif"
            check(name in z.namelist(),"missing authenticated TIFF "+name)
            with Image.open(io.BytesIO(z.read(name))) as image:
                check(image.size==(2160,1080),"climate spatial lattice drift")
                tags=image.tag_v2
                check(np.allclose(tags.get(33550,()),(1/6,1/6,0),atol=1e-14,rtol=0),
                      "unexpected climate resolution")
                check(np.allclose(tags.get(33922,()),(0,0,0,-180,90,0),atol=1e-14,rtol=0),
                      "unexpected climate spatial origin")
                check(42113 in tags,"climate nodata missing")
                vals=np.array(image)
                nodata=np.asarray(float(str(tags[42113]).strip("\x00")),dtype=vals.dtype).item()
            sampled=np.full(len(d),np.nan)
            raw=vals[rr[valid],cc[valid]]
            sampled[valid]=np.where((raw!=nodata)&np.isfinite(raw),raw.astype(float),np.nan)
            if bi in (1,5):
                check(not np.any(np.isfinite(sampled)&((sampled < -100)|(sampled>100))),
                      "impossible unmasked BIO temperature")
            if bi in (12,14):
                check(not np.any(np.isfinite(sampled)&((sampled<0)|(sampled>20000))),
                      "impossible unmasked precipitation")
            d[f"BIO{bi}"]=sampled
    return d


def all_coords_ok(d):
    return (d.latitude.between(-90,90)&d.longitude.between(-180,180)
            &d.positional_accuracy_m.between(0,5000)).to_numpy()


def distance_km(lat,lon,aa,bb):
    la=np.deg2rad(lat.astype(float));lo=np.deg2rad(lon.astype(float))
    dlat=la[aa]-la[bb];dlon=lo[aa]-lo[bb]
    h=np.sin(dlat/2)**2+np.cos(la[aa])*np.cos(la[bb])*np.sin(dlon/2)**2
    return 2*6371.0088*np.arcsin(np.sqrt(np.clip(h,0,1)))


def sampled_diversity(g):
    n=len(g)
    counts=g.morph.value_counts().reindex(PHOTO_LABELS,fill_value=0).to_numpy(int)
    if n<2:return None
    return float(1-(counts*(counts-1)).sum()/(n*(n-1)))


def coarse_second(g):
    c=g.morph.value_counts()
    a=sorted([int(c.get(k,0)) for k in PHOTO_LABELS],reverse=True)
    return a[1]/len(g) if len(g)>0 else None


def white_axis(g):
    """W_i on every valid original source taxon; no original 40/10% filtering."""
    p=np.stack(g.palette_norm9.to_numpy())
    n=len(p)
    if n<2:return (None,None,None,"n_classifiable_lt2")
    try:
        fit=two_mode_axis(p)
    except ValueError:
        return (None,None,None,"zero_continuous_displacement")
    w=float((fit.unit_axis @ Q)**2)
    return w,int(fit.cluster_sizes[0]),float(fit.minor_fraction),"estimable"


def partialcorr(y,g,e):
    """Pairwise descriptive partial Pearson correlation, NOT genetic IBD or IBE."""
    mat=np.stack((y,g,e),axis=0)
    if mat.shape[1]<20 or np.any(np.std(mat,axis=1)<1e-10):
        return None,None
    corr=np.corrcoef(mat)
    a,b,c=corr[0,1],corr[0,2],corr[1,2]
    den1=np.sqrt(max(0.,(1-b*b)*(1-c*c)))
    den2=np.sqrt(max(0.,(1-a*a)*(1-c*c)))
    if den1<=1e-10 or den2<=1e-10:return None,None
    return float((a-b*c)/den1),float((b-a*c)/den2)


def spatial_and_local(g,n_genotyped=0):
    n=len(g)
    mask=all_coords_ok(g)
    s=g.iloc[np.flatnonzero(mask)]
    nvalid=len(s)
    base={
       "n_spatial_valid_classified":int(nvalid),
       "n_local50km_colour_photo_pairs":0,
       "n_local50km_diff_observer_pairs":0,
       "n_local5km_same_year_month_diff_observer_pairs":0,
       "n_local5km_same_year_month_diff_colour_photo_pairs":0,
       "local_photo_colour_coexistence_verified_genetic":False,
       "global_vs_local_50km_colour_mismatch_delta":None,
       "global_valid_geo_mismatch_D":None,
       "local50_mismatch_D":None,
       "photographic_IBD_like_partial_r":None,
       "photographic_IBE_BIO5_like_partial_r":None,
       "h2_status":"no_local_pairs",
    }
    if nvalid<2:return base,None
    i,j=np.triu_indices(nvalid,1)
    dist=distance_km(s.latitude.to_numpy(float),s.longitude.to_numpy(float),i,j)
    labels=pd.Categorical(s.morph,categories=PHOTO_LABELS).codes.astype(np.int8)
    diff=(labels[i]!=labels[j])
    allD=float(diff.mean())
    within=(dist<=50)
    m=int(within.sum())
    base["global_valid_geo_mismatch_D"]=allD
    base["n_local50km_colour_photo_pairs"]=m
    observers=s.observer_id.to_numpy()
    obsa=observers[i];obsb=observers[j]
    distinctobs=pd.notna(obsa)&pd.notna(obsb)&(obsa!=obsb)
    base["n_local50km_diff_observer_pairs"]=int(np.sum(within&distinctobs))
    dt=s.date
    year=pd.to_numeric(dt.dt.year,errors="coerce").to_numpy()
    month=pd.to_numeric(dt.dt.month,errors="coerce").to_numpy()
    validtime=np.isfinite(year)&np.isfinite(month)
    monthpair=(validtime[i]&validtime[j]
      &(year[i]==year[j])&(month[i]==month[j])&
      (dist<=5)&distinctobs)
    base["n_local5km_same_year_month_diff_observer_pairs"]=int(monthpair.sum())
    base["n_local5km_same_year_month_diff_colour_photo_pairs"]=int(diff[monthpair].sum())
    bio5=s.BIO5.to_numpy(float)
    both=np.isfinite(bio5[i])&np.isfinite(bio5[j])
    if np.count_nonzero(both)>=20:
        a,b=partialcorr(diff[both].astype(float),np.log1p(dist[both]),
                        np.abs(bio5[i[both]]-bio5[j[both]]))
        base["photographic_IBD_like_partial_r"]=a
        base["photographic_IBE_BIO5_like_partial_r"]=b
    if m:
        local=float(diff[within].mean())
        base["local50_mismatch_D"]=local
        base["global_vs_local_50km_colour_mismatch_delta"]=allD-local
        base["h2_status"]="photo_50km_pairs_observed"
        return base, (labels,i[within],j[within],allD)
    return base,None


def climate_white_offsets(g):
    labels=g.morph.eq("white").to_numpy()
    data={"n_white_photo":int(labels.sum()),"n_nonwhite_photo":int((~labels).sum()),
          "BIO5_white_minus_nonwhite_C":None,
          "BIO12_white_minus_nonwhite_mm":None,
          "n_climate_complete":0,
          "n_white_with_climate":0,"n_nonwhite_with_climate":0,
          "H3_BIO5_status":"no_both_label_climate_support"}
    if not labels.any() or not (~labels).any():return data
    valid=np.isfinite(g.BIO5)&np.isfinite(g.BIO12)
    data["n_climate_complete"]=int(valid.sum())
    w=labels&valid.to_numpy();c=(~labels)&valid.to_numpy()
    data["n_white_with_climate"]=int(w.sum())
    data["n_nonwhite_with_climate"]=int(c.sum())
    if w.sum()<1 or c.sum()<1:return data
    for col,alias in [("BIO5","BIO5_white_minus_nonwhite_C"),
                      ("BIO12","BIO12_white_minus_nonwhite_mm")]:
        a=g[col].to_numpy(float)
        data[alias]=float(a[w].mean()-a[c].mean())
    data["H3_BIO5_status"]="at_least_one_photo_each_label"
    return data


def source_molecular_coverage(path):
    data=pd.read_csv(path,low_memory=False)
    check(len(data)==TOTAL_TAXA and data.inat_taxon_id.is_unique,
          "frozen per-4475-taxon molecular-literature coverage mismatch")
    expected={"curated_any_evidence","curated_molecular_tier_ge2","inat_taxon_id"}
    check(expected.issubset(set(data.columns)),"missing curated literature flags")
    yes=lambda s:s.fillna(False).astype(str).str.lower().isin(("true","1","yes"))
    table={}
    for tid,a,b in zip(data.inat_taxon_id,yes(data.curated_any_evidence),
                       yes(data.curated_molecular_tier_ge2)):
        table[int(tid)]=(bool(a),bool(b))
    return table


def bootstrap_genus(rows,key,n_reps=499):
    vals={}
    for x in rows:
        v=x.get(key)
        if v is None or not np.isfinite(v):continue
        genus=x["species"].split()[0]
        vals.setdefault(genus,[]).append(float(v))
    if len(vals)<5:return None
    genera=sorted(vals)
    sm=np.asarray([sum(vals[g]) for g in genera])
    ct=np.asarray([len(vals[g]) for g in genera])
    rng=np.random.default_rng(SEED+len(vals))
    choice=rng.integers(0,len(genera),size=(n_reps,len(genera)))
    means=sm[choice].sum(axis=1)/ct[choice].sum(axis=1)
    return [float(x) for x in np.quantile(means,[.025,.975])]


def summarize(rows,worlds,null_meta):
    panel=pd.DataFrame(rows)
    alln=len(panel)
    check(alln==TOTAL_TAXA and panel.inat_taxon_id.is_unique,
          "lost taxa in unified species table")
    valid=panel.global_vs_local_50km_colour_mismatch_delta.notna()
    def mean_numeric(key,subset=None):
        z=panel.loc[subset,key] if subset is not None else panel[key]
        z=pd.to_numeric(z,errors="coerce").dropna()
        return float(z.mean()) if len(z)>0 else None
    result={
      "n_species_total":alln,
      "n_photo_classifiable_each_at_least_two":int((panel.n_classified_photo>=2).sum()),
      "n_h1_continuous_axis_evaluable":int(panel.W_axis_sq_alignment.notna().sum()),
      "H1_W_descriptive_mean_all_axis_evaluable":mean_numeric("W_axis_sq_alignment"),
      "H1_W_genus_bootstrap_CI95":bootstrap_genus(rows,"W_axis_sq_alignment"),
      "H1_W_mean_on_strict_existing_gate":mean_numeric(
       "W_axis_sq_alignment",panel.H1_original_primary_gate_pass),
      "H1_n_strict_existing_gate":int(panel.H1_original_primary_gate_pass.sum()),
      "n_h2_any_local50_pair":int(valid.sum()),
      "n_h2_30plus_local_pairs_without_prior_classifiability_gate":int((panel.n_local50km_colour_photo_pairs>=30).sum()),
      "n_h2_original_30plus_pairs_AND_40plus_classifiable":int(((panel.n_local50km_colour_photo_pairs>=30)&(panel.n_classified_photo>=40)).sum()),
      "H2_mean_depletion_across_any_local_pair":mean_numeric("global_vs_local_50km_colour_mismatch_delta"),
      "H2_median_depletion_across_any_local_pair":(
         float(panel.loc[valid,"global_vs_local_50km_colour_mismatch_delta"].median()) if valid.any() else None),
      "H2_mean_depletion_original_30plus_support":mean_numeric(
          "global_vs_local_50km_colour_mismatch_delta",panel.n_local50km_colour_photo_pairs>=30),
      "n_h3_two_colours_with_complete_BIO5":int(panel.BIO5_white_minus_nonwhite_C.notna().sum()),
      "n_h3_BIO5_each_at_least_5_photo":int(((panel.n_white_with_climate>=5)
                                            &(panel.n_nonwhite_with_climate>=5)).sum()),
      "H3_BIO5_photo_mean_C_all_identifiable":mean_numeric("BIO5_white_minus_nonwhite_C"),
      "H3_BIO12_photo_mean_mm_all_identifiable":mean_numeric("BIO12_white_minus_nonwhite_mm"),
      "H3_BIO5_genus_bootstrap_CI95":bootstrap_genus(rows,"BIO5_white_minus_nonwhite_C"),
      "n_IBD_IBE_partial_correlations_estimable":0,  # exact rowwise paired n below
      "n_h4_same_month_local_distinct_observer_pairs_any":
          int((panel.n_local5km_same_year_month_diff_observer_pairs>0).sum()),
      "n_h4_same_month_local_diff_colour_labels_any":
          int((panel.n_local5km_same_year_month_diff_colour_photo_pairs>0).sum()),
      "n_h4_genetic_coexistence_verified":0,
      "n_H5_any_external_curated_literature":int(panel.H5_external_curated_any.sum()),
      "n_H5_external_molecular_tier_ge2":int(panel.H5_external_tier2_plus.sum()),
      "n_H5_individual_matched_molecular_phenotypes":0,
      "H2_technical_and_sampling_eligibility_varies_by_species":True,
      "claim_ceiling":"source-photographed within-species colour variation; not genetic morph prevalence, fitness, neutral population history, or molecular causality",
    }
    r=panel.photographic_IBD_like_partial_r.notna()&panel.photographic_IBE_BIO5_like_partial_r.notna()
    result["n_IBD_IBE_partial_correlations_estimable"]=int(r.sum())
    result["IBD_like_partial_r_mean_in_estimable"]=mean_numeric("photographic_IBD_like_partial_r")
    result["IBE_BIO5_like_partial_r_mean_in_estimable"]=mean_numeric("photographic_IBE_BIO5_like_partial_r")
    if len(worlds):
        check(len(worlds)==len(null_meta),"H2 local pair null metadata source alignment")
        a=np.stack(worlds)
        # The supplied nulls belong to the exact fixed source-eligible set
        # n_local_pairs>=1; each taxon is given equal, not pair-number, weight.
        observed=result["H2_mean_depletion_across_any_local_pair"]
        null=a.mean(axis=0)
        result["H2_any_local_pairs_null_reps"]=len(null)
        result["H2_any_local_pairs_mean_null"]=float(null.mean())
        result["H2_any_local_pairs_one_sided_p"]=float((1+int(np.count_nonzero(null>=observed)))/(len(null)+1))
        meta=pd.DataFrame(null_meta)
        check(meta.inat_taxon_id.is_unique,"H2 species vertices duplicated")
        full=panel.set_index("inat_taxon_id",verify_integrity=True)
        supports={}
        for label,filter_mask in (
          ("at_least_1_local_photo_pair",meta.n_local_pairs>=1),
          ("at_least_5_local_photo_pairs",meta.n_local_pairs>=5),
          ("at_least_30_local_photo_pairs",meta.n_local_pairs>=30),
          ("at_least_30_local_pairs_and_40plus_colours",
           (meta.n_local_pairs>=30)&(meta.n_classified_photo>=40)),
          ("at_least_100_local_photo_pairs",meta.n_local_pairs>=100),
        ):
            ix=np.flatnonzero(filter_mask.to_numpy())
            if not len(ix):
                supports[label]={"n_species":0,"state":"not_estimable"}
                continue
            y=full.loc[meta.iloc[ix].inat_taxon_id,
                       "global_vs_local_50km_colour_mismatch_delta"].to_numpy(float)
            simu=a[ix].mean(axis=0)
            x=float(y.mean())
            supports[label]={
               "n_species":int(len(ix)), "observed_mean_depletion":x,
               "null_mean":float(simu.mean()),
               "upper_p_199":float((1+int(np.count_nonzero(simu>=x)))/(PERMS+1)),
               "median_species_depletion":float(np.median(y)),
               "photo_composition_and_vertices_fixed":True,
            }
        result["H2_pair_opportunity_sensitivity"]=supports
    bycohort={}
    for cohort,sub in panel.groupby("source_cohort",sort=True):
        bycohort[cohort]={
            "n_species":int(len(sub)),
            "n_H1_axis":int(sub.W_axis_sq_alignment.notna().sum()),
            "W_mean":float(sub.W_axis_sq_alignment.mean()) if sub.W_axis_sq_alignment.notna().any() else None,
            "n_H2_local":int(sub.global_vs_local_50km_colour_mismatch_delta.notna().sum()),
            "H2_mean":float(sub.global_vs_local_50km_colour_mismatch_delta.mean())
                 if sub.global_vs_local_50km_colour_mismatch_delta.notna().any() else None,
            "n_H3_climate":int(sub.BIO5_white_minus_nonwhite_C.notna().sum()),
            "H3_BIO5_mean_C":float(sub.BIO5_white_minus_nonwhite_C.mean())
                 if sub.BIO5_white_minus_nonwhite_C.notna().any() else None,
        }
    result["by_species_disjoint_cohort"]=bycohort
    return result


def analyze(args):
    original,d=load_original(args.photos)
    d=sample_climate(d,args.worldclim_zip)
    mol=source_molecular_coverage(args.molecular_coverage)
    full=(original.groupby("inat_taxon_id",sort=True)
             .agg(species=("species","first"), source_cohort=("fcp_source_cohort","first"),
                  n_total_photo=("photo_id","size"),
                  n_distinct_observers_source=("observer_id","nunique")).reset_index())
    rows=[]
    null_worlds=[]
    null_meta=[]
    invalid_axis=[]
    for index,(tid,g) in enumerate(d.groupby("inat_taxon_id",sort=True)):
        rec=full.loc[full.inat_taxon_id.eq(tid)].iloc[0]
        w,cluster,minor,status=white_axis(g)
        n=len(g)
        wcount=int(g.morph.eq("white").sum())
        phi=float(coarse_second(g))
        prim=(n>=40 and phi>=.10 and minor is not None and minor>=.10)
        geo,edge=spatial_and_local(g)
        env=climate_white_offsets(g)
        if status!="estimable":invalid_axis.append({"taxon":int(tid),"reason":status})
        evidence=mol.get(int(tid))
        check(evidence is not None,"taxon missing from literature-only coverage")
        row={
            "inat_taxon_id":int(tid),
            "species":str(rec.species),
            "source_cohort":str(rec.source_cohort),
            "n_source_photo":int(rec.n_total_photo),
            "n_distinct_observers_source":int(rec.n_distinct_observers_source),
            "n_classified_photo":int(n),
            "n_missing_colour_label":100-int(n),
            "p_white_photo_label":float(wcount/n),
            "n_photo_white":wcount,
            "n_photo_nonwhite":n-wcount,
            "D_photo_pairwise_unbiased":sampled_diversity(g),
            "W_axis_sq_alignment":w,
            "H1_two_means_cluster_size_a":cluster,
            "H1_two_means_minor_fraction":minor,
            "H1_original_primary_gate_pass":bool(prim),
            "H1_axis_identifiability":status,
            "H5_external_curated_any":bool(evidence[0]),
            "H5_external_tier2_plus":bool(evidence[1]),
            "H5_same_individual_genetic_mechanism_verified":False,
            "H4_genetically_confirmed_local_coexistence":False,
            **geo,**env,
        }
        rows.append(row)
        if edge is not None:
            labels,aa,bb,gD=edge
            local_m=geo["local50_mismatch_D"]
            rng=np.random.default_rng(SEED+int(tid))
            null=np.empty(PERMS,dtype=float)
            for k in range(PERMS):
                perm=labels[rng.permutation(len(labels))]
                null[k]=gD-float(np.mean(perm[aa]!=perm[bb]))
            null_worlds.append(null)
            null_meta.append({"inat_taxon_id":int(tid),"n_local_pairs":int(geo["n_local50km_colour_photo_pairs"]),"n_classified_photo":int(n)})
        if (index+1)%750==0:
            print("FCP_ALL4475_ANALYSIS_TAXA_COMPLETED",index+1,flush=True)
    check(len(rows)==TOTAL_TAXA,"lost species in analysis grouping")
    result=summarize(rows,null_worlds,null_meta)
    result.update({
        "schema":"fcp_entire4475_same_source_hypothesis_measurement_v1",
        "complete_input_authentic_sha":SOURCE_SHA,
        "complete_input_climate_zip_sha":BIO_ARCHIVE_SHA,
        "every_species_present":True,
        "all_photo_measurement_rows":TOTAL_ROWS,
        "all_taxa":TOTAL_TAXA,
        "source_classified_photo_rows":CLASSIFIABLE_ROWS,
        "n_axis_nonidentifiable":len(invalid_axis),
        "n_photo_label_counts_below_40":sum(r["n_classified_photo"]<40 for r in rows),
        "no_source_row_excluded_from_roster":True,
        "null_resampling_only_when_at_least_one_50km_photo_pair":True,
        "all_random_seed":SEED,
        "no_molecular_causation_for_source_photo_individuals":True,
        "seven_macro_taxon_H1_to_H5_are_overlapping_estimands_not_unbiased_prevalence":True,
    })
    out=Path(args.out)
    out.mkdir(parents=True,exist_ok=True)
    tbl=pd.DataFrame(rows).sort_values(["source_cohort","inat_taxon_id"])
    # No source coordinates, independent observation or person identifiers.
    tbl.to_csv(out/"all_4475_species_H1_H5_photo_climate_molecular_scope.csv",
               index=False,float_format="%.12g")
    (out/"all_4475_species_hypotheses_summary.json").write_text(
      json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("FCP_REAL_ALL4475_H1_H5_CORE_COMPLETE",json.dumps(result,sort_keys=True,allow_nan=False))
    return result


def synthetic():
    # D estimates probability of colour mismatch for different sampled photos.
    labels=pd.DataFrame({"morph":["white","white","blue_purple","blue_purple"]})
    check(abs(sampled_diversity(labels)-2/3)<1e-12,"unbiased D must use n(n-1)")
    a=np.array([0.,0.]); b=np.array([1.,0.]); dist=distance_km(np.array([0.,0.]),np.array([0.,1.]),
                                                               np.array([0]),np.array([1]))
    check(abs(dist[0]-111.195)<0.5,"haversine fixture drift")
    p=np.array([[.9,.1]+[0.]*7,[.1,.9]+[0.]*7,
                [.9,.1]+[0.]*7,[.1,.9]+[0.]*7])
    fit=two_mode_axis(p)
    check(abs(np.linalg.norm(fit.unit_axis)-1)<1e-12,"W axis not unit normalized")
    e=np.array([1.,2.,3.,4.,5.,6.])
    check(np.isfinite(e).all(),"finite climate fixture")
    print("FCP_ALL4475_UNIFIED_SYNTHETIC_PASS")


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--photos",type=Path)
    p.add_argument("--worldclim-zip",type=Path)
    p.add_argument("--molecular-coverage",type=Path)
    p.add_argument("--out",type=Path)
    p.add_argument("--self-test",action="store_true")
    args=p.parse_args()
    synthetic()
    if args.self_test:return
    check(all((args.photos,args.worldclim_zip,args.molecular_coverage,args.out)),
          "missing mandatory source path")
    analyze(args)


if __name__=="__main__":
    main()
