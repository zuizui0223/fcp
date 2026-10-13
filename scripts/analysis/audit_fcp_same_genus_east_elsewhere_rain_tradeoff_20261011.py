#!/usr/bin/env python3
"""Source-locked ecological falsification: does rain/flower-colour transport among the SAME genera in east and other longitudes?

Not a generic global colour predictor. The comparison deliberately retains only
botanical genera with independent sampled species in both geographic sectors.
Three further controls: common observed rainfall support, photo-quality-only
samples, and original SAME named species re-photographed in both sectors.

This is post-hoc original iNaturalist photo colour, not organ-verified pigment,
genotypic morphs, fitness or an identified climatic selective trade-off.
"""
from __future__ import annotations
import argparse, hashlib, json, itertools
from pathlib import Path
import numpy as np
import pandas as pd
from audit_fcp_full42111_rain_white_bioregion_20261011 import (
    get_source, eligible, species_means, COLS, HUES, SHA
)

RANDOM_SEED = 2026101103
BOOTSTRAPS = 1999
MIN_REGION_SPECIES = (3, 5, 10)
SCHEMA = "fcp_global42111_identical_genus_east_rest_abiotic_flower_photo_v1"


def shared_genus_subset(c: pd.DataFrame, min_sp: int) -> tuple[pd.DataFrame, int]:
    """Each region must contain >=min_sp DISTINCT original source taxa in SAME genus."""
    d = c.copy()
    d["east"] = d.longitude.ge(90.) & d.longitude.le(180.)
    counts = d.groupby(["genus", "east"]).inat_taxon_id.nunique().unstack(fill_value=0)
    if not {True, False}.issubset(counts.columns):
        raise ValueError("missing one geographic side of source photographs")
    names = counts.index[(counts[True] >= min_sp) & (counts[False] >= min_sp)]
    out = d.loc[d.genus.isin(names)].copy()
    if out.empty: return out, 0
    if out.groupby(["genus", "east"]).inat_taxon_id.nunique().min() < min_sp:
        raise RuntimeError("a genus is not supported by distinct species on both sides")
    return out, len(names)


def sufficient_statistics(source: pd.DataFrame, mu: np.ndarray, sd: np.ndarray,
                          minimum_species: int, equal_genus: bool) -> dict:
    """Within-genus species-equal photographic hue and abiotic slopes."""
    sp = species_means(source, mu, sd)
    groups: dict[str, tuple[np.ndarray,np.ndarray,int]] = {}
    if len(sp) == 0: return groups
    for genus, g in sp.groupby("genus", sort=True):
        if len(g) < minimum_species: continue
        xx = g[list(COLS)].to_numpy(float)
        yy = g[list(HUES)].to_numpy(float)
        if not np.isfinite(xx).all() or not np.isfinite(yy).all():
            raise RuntimeError("unexpected original source photo climate or hue missing")
        xx = xx - xx.mean(axis=0)
        yy = yy - yy.mean(axis=0)
        weight = 1.0 / len(g) if equal_genus else 1.0
        groups[str(genus)] = (weight*(xx.T@xx), weight*(xx.T@yy), len(g))
    return groups


def paired_genus_bootstrap(a: dict, b: dict, *, repetitions: int,
                           seed: int, include_white_blue_contrast: bool=False) -> dict:
    """Resample IDENTICAL genera on both sides to estimate response difference."""
    common = sorted(set(a) & set(b))
    n = len(common)
    if n < 30: return {"status":"HOLD_FEWER_THAN_30_SHARED_GENERA","n_shared_genera":n}
    mat_a = np.stack([a[g][0] for g in common])
    mat_b = np.stack([b[g][0] for g in common])
    rhs_a = np.stack([a[g][1] for g in common])
    rhs_b = np.stack([b[g][1] for g in common])
    m1, m2 = mat_a.sum(axis=0),mat_b.sum(axis=0)
    if np.linalg.cond(m1)>1e10 or np.linalg.cond(m2)>1e10:
        return {"status":"HOLD_ILL_CONDITIONED","n_shared_genera":n}
    base_a = np.linalg.solve(m1, rhs_a.sum(axis=0))[0,:]
    base_b = np.linalg.solve(m2, rhs_b.sum(axis=0))[0,:]
    rng = np.random.default_rng(seed)
    draws = rng.multinomial(n, np.repeat(1/n,n), size=repetitions)
    matrices_a = np.einsum("rn,nij->rij",draws,mat_a,optimize=True)
    matrices_b = np.einsum("rn,nij->rij",draws,mat_b,optimize=True)
    ys_a = np.einsum("rn,nij->rij",draws,rhs_a,optimize=True)
    ys_b = np.einsum("rn,nij->rij",draws,rhs_b,optimize=True)
    try:
        v_a = np.linalg.solve(matrices_a,ys_a)[:,0,:]
        v_b = np.linalg.solve(matrices_b,ys_b)[:,0,:]
    except np.linalg.LinAlgError:
        return {"status":"HOLD_SINGULAR_BOOTSTRAP","n_shared_genera":n}
    result = {}
    for j,hue in enumerate(HUES):
        diff = v_a[:,j] - v_b[:,j]
        result[hue] = {
          "east_rainfall_slope":float(base_a[j]),
          "elsewhere_rainfall_slope":float(base_b[j]),
          "east_minus_elsewhere_rainfall_slope":float(base_a[j]-base_b[j]),
          "east_genus_resample_CI95":np.quantile(v_a[:,j],[.025,.975]).tolist(),
          "elsewhere_genus_resample_CI95":np.quantile(v_b[:,j],[.025,.975]).tolist(),
          "paired_same_genus_difference_CI95":np.quantile(diff,[.025,.975]).tolist(),
          "one_sided_p_difference_not_positive":float((1+(diff<=0).sum())/(repetitions+1))
        }
    # The default return is IDENTICAL to the 2026-10-11 parent result.
    # A separately requested, post-outcome colour contrast uses the SAME
    # matched genus draws (thus retaining cross-hue bootstrap covariance).
    out={"status":"OK","n_shared_genera":n,
      "n_original_distinct_species_E":int(sum(a[g][2] for g in common)),
      "n_original_distinct_species_elsewhere":int(sum(b[g][2] for g in common)),
      "bootstrap_same_genus_replicates":repetitions,
      "response_four_original_photo_hues":result}
    if include_white_blue_contrast:
        wi=HUES.index("white"); bi=HUES.index("blue_purple")
        paired_region=v_a-v_b
        sampled=(paired_region[:,wi]-paired_region[:,bi])
        observed=float((base_a[wi]-base_b[wi])-(base_a[bi]-base_b[bi]))
        names=list(itertools.combinations(range(len(HUES)),2))
        pair_point=np.array([(base_a[i]-base_b[i])-(base_a[j]-base_b[j])
                             for i,j in names])
        pair_samples=np.column_stack([paired_region[:,i]-paired_region[:,j]
                                      for i,j in names])
        pair_sd=pair_samples.std(axis=0,ddof=1)
        keep=pair_sd>1e-12
        max_t=0.0 if not keep.any() else float(np.quantile(
            np.max(np.abs((pair_samples[:,keep]-pair_point[keep])/pair_sd[keep]),axis=1),.95))
        ix=names.index((wi,bi))
        white_ci=result["white"]["paired_same_genus_difference_CI95"]
        blue_ci=result["blue_purple"]["paired_same_genus_difference_CI95"]
        out["posthoc_white_blue_region_rain_contrast"]={
          "definition":"(white-blue)_east minus (white-blue)_elsewhere",
          "unit":"photo_colour_probability_per_global_BIO12_SD",
          "point":observed,
          "east_white_minus_blue":float(base_a[wi]-base_a[bi]),
          "elsewhere_white_minus_blue":float(base_b[wi]-base_b[bi]),
          "paired_genus_bootstrap_percentile_CI95":np.quantile(sampled,[.025,.975]).tolist(),
          "bootstrap_positive_fraction":float((sampled>0).mean()),
          "bootstrap_nonpositive_fraction":float((sampled<=0).mean()),
          "bonferroni_conservative_one_sided_95pct_lower":float(white_ci[0]-blue_ci[1]),
          "bonferroni_conservative_one_sided_95pct_upper":float(white_ci[1]-blue_ci[0]),
          "six_hue_pair_maxT_critical_95":max_t,
          "six_pair_simultaneous_bootstrap_approx_CI95":[
            float(observed-max_t*pair_sd[ix]),float(observed+max_t*pair_sd[ix])],
          "six_pair_family_only_region_not_multiplicity_adjusted":True,
          "retrospective_selected_east_sector_and_colour_pair":True,
          "not_a_null_hypothesis_p_value":True}
    return out


def paired_same_species(source: pd.DataFrame, mu:np.ndarray,sd:np.ndarray,
                        min_photo_each:int) -> dict:
    """Focal source-species intercept eliminated by east-minus-other differences."""
    f=source.copy()
    f["east"] = f.longitude.ge(90.) & f.longitude.le(180.)
    grouping = ["inat_taxon_id", "east"]
    a={v:"mean" for v in (*COLS,"white")}
    a.update({"photo_id":"size","genus":"first"})
    grouped=f.groupby(grouping,sort=True,as_index=False).agg(a).rename(columns={"photo_id":"n_photos"})
    supp=grouped.groupby("inat_taxon_id").agg(n_regions=("east","size"),
                                               minimum_photos=("n_photos","min"))
    ids=supp.index[(supp.n_regions==2)&(supp.minimum_photos>=min_photo_each)]
    matched=grouped[grouped.inat_taxon_id.isin(ids)]
    east=matched[matched.east].set_index("inat_taxon_id").sort_index()
    rest=matched[~matched.east].set_index("inat_taxon_id").sort_index()
    if len(east)<100 or not east.index.equals(rest.index):
        return {"status":"HOLD_INSUFFICIENT_MATCHED_SOURCE_SPECIES","n_species":len(east)}
    X=np.column_stack((np.ones(len(east)),
           (east[list(COLS)].to_numpy(float)-rest[list(COLS)].to_numpy(float))/sd[None,:]))
    yy=(east.white-rest.white).to_numpy(float)
    beta=np.linalg.lstsq(X,yy,rcond=None)[0]
    residual=yy-X@beta
    genera,ind=np.unique(east.genus.to_numpy(str),return_inverse=True)
    k=len(genera);n=len(yy);p=X.shape[1]
    inv=np.linalg.inv(X.T@X)
    scores=np.zeros((k,p))
    np.add.at(scores,ind,X*residual[:,None])
    vcov=inv@(scores.T@scores)@inv*(k/(k-1))*((n-1)/(n-p))
    se=np.sqrt(np.maximum(0,np.diag(vcov)))
    j=1
    return {"status":"OK","n_same_taxa_both_geographies":n,
      "n_genera":k,"n_original_photo_rows":int(matched.n_photos.sum()),
      "source_species_fixed_comparison":True,
      "delta_white_E_minus_elsewhere_mean":float(yy.mean()),
      "annual_rainfall_BIO12_species_paired_slope":float(beta[j]),
      "annual_rainfall_BIO12_genus_cluster_CI95":[float(beta[j]-1.96*se[j]),float(beta[j]+1.96*se[j])],
      "no_true_species_population_genotype_adaptation_test":True}


def run(path:Path,repetitions:int=BOOTSTRAPS,include_white_blue_contrast:bool=False)->dict:
    source=get_source(path)
    opportunities,c,roi,mu,sd=eligible(source)
    if len(source)!=100543 or len(c)!=45836 or c.inat_taxon_id.nunique()!=25028:
        raise RuntimeError("Original complete FCP photographed population drifted")
    if not c.morph.isin(HUES).all():
        raise RuntimeError("focal source has invalid four-state photo hue")
    c=c.copy()
    c["east"]=c.longitude.ge(90.) & c.longitude.le(180.)
    # Rainfall common support is fixed independently of colour from actual photo sites.
    qa=c[c.east].wc_bio12.quantile([.1,.9])
    qb=c[~c.east].wc_bio12.quantile([.1,.9])
    support=[float(max(qa.loc[.1],qb.loc[.1])),float(min(qa.loc[.9],qb.loc[.9]))]
    if not 0<support[0]<support[1]:raise RuntimeError("No shared source photo rainfall support")
    qarea=float(c.flower_effective_pixels.median())
    qnuis=float(c.nuisance_pixel_fraction.median())
    cases={
       "all_classified_original_source_photos":c,
       "common_rainfall_10_90_both_regions":c.loc[c.wc_bio12.between(*support)].copy(),
       "original_photo_quality_strict":c.loc[(c.flower_effective_pixels>=qarea)&
                                               (c.nuisance_pixel_fraction<=qnuis)].copy()
    }
    comparisons={}
    for label,frame in cases.items():
        comparisons[label]={}
        for k in MIN_REGION_SPECIES:
            chosen,ng=shared_genus_subset(frame,k)
            if ng<30:
                comparisons[label][str(k)]={"status":"HOLD_INSUFFICIENT_GENERA",
                  "n_shared_genera":ng};continue
            n_species=int(chosen.inat_taxon_id.nunique())
            n_photos=len(chosen)
            both_taxa=int(chosen.groupby("inat_taxon_id").east.nunique().eq(2).sum())
            sides={}
            for tag,equal_genus in (("equal_species",False),("equal_genera",True)):
                left=sufficient_statistics(chosen[chosen.east],mu,sd,k,equal_genus)
                right=sufficient_statistics(chosen[~chosen.east],mu,sd,k,equal_genus)
                v=paired_genus_bootstrap(left,right,repetitions=repetitions,
                                          seed=RANDOM_SEED + k + (0 if not equal_genus else 100)
                                          + {"all_classified_original_source_photos":0,
                                             "common_rainfall_10_90_both_regions":1000,
                                             "original_photo_quality_strict":2000}[label],
                                          include_white_blue_contrast=include_white_blue_contrast)
                sides[tag]=v
            comparisons[label][str(k)]={
                "source_qualified_shared_genera":ng,
                "n_unique_source_species":n_species,
                "n_original_photo_records":n_photos,
                "n_same_named_species_observed_both_regions":both_taxa,
                "models":sides,
                "species_composition_changed_inside_shared_genus":True
            }
    paired={str(k):paired_same_species(c,mu,sd,k) for k in (1,2,3,5)}
    return {"schema":SCHEMA,"date_jst":"2026-10-11",
      "source_archive_sha256":SHA,
      "original_unique_photo_opportunities":len(source),
      "original_4state_classified_photo_ids":int(source.classifiable.sum()),
      "original_classifiable_geoclimate_photo_ids":len(c),
      "original_classifiable_geoclimate_species":int(c.inat_taxon_id.nunique()),
      "source_raw_rainfall_SD_mm":float(sd[0]),
      "geography":{"east_definition":"photo site longitude >=90 and <=180 E",
        "elsewhere_definition":"all other original geographic longitude sectors",
        "east_photo_rainfall_quantile_10_90_mm":qa.tolist(),
        "elsewhere_photo_rainfall_quantile_10_90_mm":qb.tolist(),
        "common_rainfall_interval_mm":support},
      "quality_photo_medians":{"flower_effective_pixels":qarea,"nuisance_fraction":qnuis},
      "environment_adjustments":list(COLS),
      "n_genus_bootstrap_replicates":repetitions,
      "same_genera_across_geographic_sectors":comparisons,
      "same_species_across_geographic_sectors":paired,
      "hard_nonclaims":[
       "Same genus in both sectors is not same species, dated phylogenetic control or shared pigment biochemical mechanism.",
       "Source flora photographic labels are not individual genetic pigmentation/fitness payoffs.",
       "Sector 90-180E singled out after seeing prior exploratory signals: difference is posthoc, not prospective transport.",
       "Retrospective source photographs and classification opportunity vary by geography and background quality.",
       "Same species photographed in both sectors is not replicated reciprocal-transplant/genetic environment reaction.",
       "Even a same-genus geographically varying slope cannot establish climatic colour adaptation, herbivore protection, pollinator payoff or pigment synthesis costs.",
       "No raw original FCP photo colour or independent H1/H2 prospective data changed."
      ],"frozen_original_main_H1_H2_unchanged":True}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-archive",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    p.add_argument("--bootstrap",type=int,default=BOOTSTRAPS)
    p.add_argument("--posthoc-white-blue",action="store_true",help="Opt-in derived paired genus posthoc contrast; default parent output unchanged")
    a=p.parse_args()
    z=run(a.source_archive,repetitions=a.bootstrap,
          include_white_blue_contrast=a.posthoc_white_blue)
    a.outdir.mkdir(parents=True,exist_ok=True)
    (a.outdir/"result.json").write_text(json.dumps(z,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({
      "global_classified_source_photos":z["original_classifiable_geoclimate_photo_ids"],
      "common_rainfall_interval_mm":z["geography"]["common_rainfall_interval_mm"],
      "shared_genus_cases":z["same_genera_across_geographic_sectors"],
      "same_species_paired":z["same_species_across_geographic_sectors"]
    },indent=2))
if __name__=="__main__":main()
