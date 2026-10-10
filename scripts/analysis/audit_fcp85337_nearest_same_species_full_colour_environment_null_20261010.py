#!/usr/bin/env python3
"""Exact source same-species nearest-site FOUR-COLOUR abiotic permutation.

Freeze ONE closest pair of public-source photographed sites per species using
only original photo IDs/coordinates. For 50/100/250/500km capped pairs,
orthogonalize signed environment differences to signed geographic/elevation
differences using covariates only, not source flower-colour labels.

Within those same source-species pairs, swap original four-state photo colours
between the two locations under fixed 999 signs; test all environmental
variables with max-statistic FWER, never report favourable single features
alone. No genetic/fitness or causal adaptation inference.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd

from compare_fcp_original85337_species_fixed_multiabiotic_20261010 import (
    match_original,populations,
)
from compare_fcp_global42111_full_climate_sun_soil_blocks_20261010 import (
    BLOCKS,CLASSES,
)
from audit_fcp85337_same_species_microgeographic_exchange_20261010 import (
    geometry,EARTH_KM,
)

SCHEMA="fcp_original85337_nearest_same_species_four_colour_geo_adjusted_abiotic_sign_null_v1"
RADII_KM=(50,100,250,500)
N_PERM=999
SEED=20261010
MIN_ORIGINAL_PAIRS=60
MIN_ORIGINAL_COLOUR_MISMATCH_PAIRS=30
MIN_ORIGINAL_GENERA=20
GEO_DIFF=("d_lat","d_abslat","d_lon_sin","d_lon_cos","d_elev")
ENV_GROUPS={k:tuple(v) for k,v in BLOCKS.items() if k!="elevation"}
FEATURES_CLIMATE=tuple(x for k,v in ENV_GROUPS.items() if k!="soil" for x in v)
FEATURES_SOIL=FEATURES_CLIMATE+ENV_GROUPS["soil"]
ALL_FIELDS=tuple(set(FEATURES_SOIL)|{"wc_elevation_m"})


def deterministic_closest_original_pairs(source:pd.DataFrame)->pd.DataFrame:
    required={"inat_taxon_id","photo_id","observation_id","cell_id","species",
              "morph","latitude","longitude",*ALL_FIELDS}
    if not required.issubset(source):
        raise ValueError("Original source original-photo local-climate fields unavailable")
    if source.photo_id.duplicated().any() or source[["inat_taxon_id","cell_id"]].duplicated().any():
        raise ValueError("Original one-photo per taxon×geographical cell identity changed")
    rows=[]
    for taxon,p in source.groupby("inat_taxon_id",sort=True):
        if len(p)<2:continue
        p=p.sort_values("photo_id",kind="stable").reset_index(drop=True)
        distance=geometry(p.latitude.to_numpy(float),p.longitude.to_numpy(float))
        ii,jj=np.triu_indices(len(p),1)
        k=int(np.argmin(distance[ii,jj]))
        a=p.iloc[int(ii[k])];b=p.iloc[int(jj[k])]
        if a.cell_id==b.cell_id or a.photo_id==b.photo_id:
            raise RuntimeError("Same original photographed species-cell erroneously paired")
        la,lb=float(a.latitude),float(b.latitude)
        loa,lob=float(a.longitude),float(b.longitude)
        result={
            "inat_taxon_id":int(taxon),
            "genus":str(a.species).split()[0],
            "photo_id_a":int(a.photo_id),
            "photo_id_b":int(b.photo_id),
            "cell_id_a":int(a.cell_id),
            "cell_id_b":int(b.cell_id),
            "morph_a":str(a.morph),
            "morph_b":str(b.morph),
            "mismatch":int(str(a.morph)!=str(b.morph)),
            "distance_km":float(distance[ii[k],jj[k]]),
            "d_lat":lb-la,
            "d_abslat":abs(lb)-abs(la),
            "d_lon_sin":float(np.sin(np.deg2rad(lob))-np.sin(np.deg2rad(loa))),
            "d_lon_cos":float(np.cos(np.deg2rad(lob))-np.cos(np.deg2rad(loa))),
            "d_elev":float(b.wc_elevation_m-a.wc_elevation_m),
        }
        for name in FEATURES_SOIL:
            av=pd.to_numeric(a[name],errors="coerce")
            bv=pd.to_numeric(b[name],errors="coerce")
            result["d_"+name]=float(bv-av) if np.isfinite(av) and np.isfinite(bv) else np.nan
        rows.append(result)
    d=pd.DataFrame(rows)
    if len(d) and (d.inat_taxon_id.duplicated().any() or
                   d[["photo_id_a","photo_id_b"]].stack().duplicated().any()):
        raise RuntimeError("Original source photographic species was used in multiple fixed pairs")
    return d


def geography_residualized_gradients(subset:pd.DataFrame,fields:tuple[str,...])->np.ndarray:
    """Remove only identifiable linear signed geographical displacement.

    Regress signed environmental contrasts against signed lat, |lat|,
    longitude sine/cos and real photographed-site elevation contrasts across
    source pairs, without consulting original flower-colour labels.
    """
    X=subset[["d_"+k for k in fields]].to_numpy(float)
    Z=subset[list(GEO_DIFF)].to_numpy(float)
    if not np.isfinite(X).all() or not np.isfinite(Z).all():
        raise ValueError("Source pair model missing true-site climate/geographic differences")
    mu=Z.mean(axis=0);sd=Z.std(axis=0)
    sd[sd<1e-9]=1.
    z=(Z-mu)/sd
    z=np.column_stack((np.ones(len(z)),z))
    coef=np.linalg.solve(z.T@z+np.diag([0.,1e-5,1e-5,1e-5,1e-5,1e-5]),z.T@X)
    residual=X-z@coef
    rms=np.sqrt(np.mean(np.square(residual),axis=0))
    good=rms>=1e-8
    residual[:,good]=residual[:,good]/rms[good]
    residual[:,~good]=0.
    return residual


def sign_null(subset:pd.DataFrame,fields:tuple[str,...],*,nperm:int=N_PERM,seed:int=SEED)->dict:
    observed_color=np.asarray([CLASSES.index(x) for x in subset.morph_a],int)
    second_color=np.asarray([CLASSES.index(x) for x in subset.morph_b],int)
    yy=np.eye(4)[second_color]-np.eye(4)[observed_color]
    X=geography_residualized_gradients(subset,fields)
    weighted=X[:,:,None]*yy[:,None,:]
    norm=np.sqrt(np.square(weighted).sum(axis=(0,2)))
    good=norm>1e-8
    observed=np.zeros(len(fields))
    observed[good]=np.linalg.norm(weighted.sum(axis=0)[good],axis=1)/norm[good]
    rng=np.random.default_rng(seed)
    signs=rng.choice(np.array([-1.,1.]),size=(nperm,len(subset)),replace=True)
    perm_sums=np.einsum("rn,nfc->rfc",signs,weighted,optimize=True)
    null=np.zeros((nperm,len(fields)),float)
    null[:,good]=np.linalg.norm(perm_sums[:,good,:],axis=2)/norm[None,good]
    allmax=null.max(axis=1)
    effects={}
    for j,feature in enumerate(fields):
        t=float(observed[j])
        p_unadjusted=float((1+np.sum(null[:,j]>=t-1e-12))/(nperm+1))
        p_max=float((1+np.sum(allmax>=t-1e-12))/(nperm+1))
        effects[feature]={
            "source_signed_gradient_multiclass_photo_statistic":t,
            "one_sided_label_swap_p_unadjusted":p_unadjusted,
            "within_radius_maxT_across_all_individual_features_FWER_p":p_max,
            "conservative_8_nested_cohort_radius_tests_Bonferroni_maxT_p":min(1.,p_max*8),
            "signed_photo_colour_components":{c:float(v) for c,v in zip(CLASSES,weighted.sum(axis=0)[j])},
            "zero_residualized_environmental_gradient":bool(not good[j]),
        }
    groups={}
    for g,features in ENV_GROUPS.items():
        members=[fields.index(f) for f in features if f in fields]
        if not members:continue
        stat=float(max(observed[members]))
        null_max=np.max(null[:,members],axis=1)
        p=float((1+np.sum(null_max>=stat-1e-12))/(nperm+1))
        groups[g]={
            "max_multiclass_gradient_statistic_within_physical_block":stat,
            "within_block_label_swap_p":p,
            "maxT_over_ALL_environment_features_FWER_p":
                float((1+np.sum(allmax>=stat-1e-12))/(nperm+1)),
            "min_distance_and_cohort_adjusted_p":min(
                1.,8*float((1+np.sum(allmax>=stat-1e-12))/(nperm+1))),
            "n_environmental_features":len(members),
        }
    return {
        "n_fixed_original_source_photo_pairs":len(subset),
        "n_four_colour_photograph_mismatched_pairs":int((observed_color!=second_color).sum()),
        "n_different_original_source_genera":int(subset.genus.nunique()),
        "n_source_photo_species":int(subset.inat_taxon_id.nunique()),
        "geo_difference_residualization_without_photo_colour_outcomes":True,
        "one_photo_pair_per_species_fixed_by_nearest_actual_photo_site_distance":True,
        "n_within_pair_photo_colour_label_swap_nulls":nperm,
        "preserves_exact_source_species_pair_and_two_photo_colour_counts":True,
        "environmental_feature_results":effects,
        "all_named_environmental_block_results":groups,
        "multiple_testing":"Within radius maxT over ALL individual variables, and conservative x8 over nested radii and climate/soil cohorts; posthoc exploratory",
    }


def run(original:pd.DataFrame,site:pd.DataFrame,*,strict=True,nperm:int=N_PERM)->dict:
    linked,original_receipt=match_original(original,site,strict=strict)
    pools,coverage=populations(linked)
    if strict and (original_receipt["source_original_taxon_cell_photo_rows"]!=85337 or
                   original_receipt["source_four_state_colour_classifiable"]!=39075):
        raise RuntimeError("Original source full photo labels changed")
    climate=pools["CLIMATE_SOURCE_ONLY"]
    paired=deterministic_closest_original_pairs(climate)
    if paired.empty:
        raise RuntimeError("No original repeated photos can be paired")
    studies={}
    for radius in RADII_KM:
        within=paired.loc[paired.distance_km.le(radius+1e-8)]
        categories={}
        for name,features in (("climate_without_soil",FEATURES_CLIMATE),
                              ("climate_and_soil",FEATURES_SOIL)):
            d=within
            if name=="climate_and_soil":
                d=d.loc[d[["d_"+v for v in BLOCKS["soil"]]].notna().all(axis=1)]
            count=len(d)
            mismatch=int(d.mismatch.sum())
            out={
                "n_all_original_photo_climate_pairs_before_radius":len(paired),
                "n_original_source_species_pairs_at_radius":count,
                "n_original_source_photo_colour_mismatched_pairs":mismatch,
                "n_distinct_nominal_source_genera":int(d.genus.nunique()),
                "n_original_source_photo_pairs_with_same_four_state":count-mismatch,
                "status":"HOLD_INSUFFICIENT_INTRASPECIFIC_NEARBY_FOUR_COLOUR_PAIR_SUPPORT",
            }
            if count>=MIN_ORIGINAL_PAIRS and mismatch>=MIN_ORIGINAL_COLOUR_MISMATCH_PAIRS and d.genus.nunique()>=MIN_ORIGINAL_GENERA:
                out["status"]="ORIGINAL_SOURCE_GEO_RESIDUALIZED_FOUR_COLOUR_PAIR_PERMUTATION"
                out["null"]=sign_null(d,features,nperm=nperm,seed=SEED+int(radius)+(1 if name=="climate_and_soil" else 0))
            categories[name]=out
        studies[str(int(radius))]=categories
    return {
        "schema":SCHEMA,
        "date_jst":"2026-10-10",
        "status":"SOURCE_ORIGINAL_SAME_SPECIES_NEAREST_SITE_FOUR_COLOUR_GEOGRAPHIC_PAIR_EXPLORATION",
        "original_photo_ledger":original_receipt,
        "original_abiotic_coverage":coverage,
        "original_climate_complete_nearest_pair_species":len(paired),
        "distance_caps_km_predefined_all_reported":list(RADII_KM),
        "source_nearest_pair_selection_only_by_photo_ID_and_real_coordinates":True,
        "source_photo_colour_outcomes_unchanged":True,
        "source_4colour_phenotype_orientation_conditional_label_swap":True,
        "no_independent_1499_or_2000_730_taxa_opened":True,
        "all_environment_variables_reported_without_selecting_a_winner":True,
        "studies":studies,
        "hard_nonclaims":[
            "Nearest across-source-cell photographed species pairs are not replicate observations of the same plant genotype or breeding population",
            "Linear signed lat/lon/elevation difference residualization and <=50/100/250/500km cap are not a full spatial random field",
            "Conditional source label swap assumes exchangeability of photographed floral-colour orientation at two old sites and is not a causal randomization trial",
            "Mismatched original four-state photograph categories do not prove an inherited pigment gradient or adaptive fitness effect",
            "Any test is retrospective; nested distances and soil missingness cohorts are overlapping not independent confirmations",
            "Multiple 35 environmental features are maxT-corrected within each radius and conservatively Bonferroni corrected across 8 nested comparisons",
        ],
    }


def main():
    a=argparse.ArgumentParser()
    a.add_argument("--original-taxon-cell",required=True,type=Path)
    a.add_argument("--expanded-all-original-photo-sites",required=True,type=Path)
    a.add_argument("--outdir",required=True,type=Path)
    x=a.parse_args()
    result=run(pd.read_csv(x.original_taxon_cell,low_memory=False),
               pd.read_csv(x.expanded_all_original_photo_sites,low_memory=False))
    x.outdir.mkdir(parents=True,exist_ok=True)
    (x.outdir/"result.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "schema":SCHEMA,
        "original_source":result["original_photo_ledger"],
        "studies":{k:{name:{
            "status":v["status"],"n":v["n_original_source_species_pairs_at_radius"],
            "n_mismatch":v["n_original_source_photo_colour_mismatched_pairs"],
            "environment_block_null":v.get("null",{}).get("all_named_environmental_block_results")
        } for name,v in cohort.items()} for k,cohort in result["studies"].items()}
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
