#!/usr/bin/env python3
"""Independent downstream audit of authenticated original-photo mode/CV outcomes.

Reads the artifact produced by GitHub Actions run 39923312111. No photographs,
observer identifiers, or exact coordinate fields are expected in input/output.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd

EXPECTED_N = 1761
EXPECTED_MODE = {"not_resolved_as_multimodal":1691,
                 "unstable_photo_mode_candidate":43,
                 "recurrent_photographic_mode_candidate":27}
EXPECTED_GEO = {"geography_underidentified":1312,
                "no_clear_predictive_geometry":440,
                "transferable_broad_cline_candidate":8,
                "local_nonlinearity_only_candidate":1}


def audit(d):
    assert len(d) == EXPECTED_N
    assert d.inat_taxon_id.is_unique
    assert d.mode_status.value_counts().to_dict() == EXPECTED_MODE
    assert d.spatial_signature.value_counts().to_dict() == EXPECTED_GEO
    assert int(d.spatial_block_estimable.sum()) == 449
    assert int(d.observer_estimable.sum()) == 1566
    assert d.n_original_classifiable.ge(40).all()
    assert d.species.astype(str).str.len().min() > 0
    assert all(not any(x in c for x in ("photo_id", "observer_id",
                                         "latitude", "longitude"))
               for c in d.columns)
    def num(x):
        x=pd.to_numeric(x, errors="coerce").dropna()
        if len(x)==0:return None
        return {"n":int(len(x)), "mean":float(x.mean()),
                "median":float(x.median()),
                "q025":float(x.quantile(.025)),
                "q975":float(x.quantile(.975)),
                "fraction_gt_0":float((x>0).mean()),
                "fraction_gt_0p08":float((x>.08).mean())}
    transfer=d.loc[d.spatial_block_estimable].copy()
    observer=d.loc[d.observer_estimable].copy()
    mode=d.loc[d.mode_status=="recurrent_photographic_mode_candidate"].copy()
    candidate_columns=["inat_taxon_id","species","genus","n_original_classifiable",
                       "n_spatial_classifiable","photographic_coarse_label_count",
                       "photographic_pc1_explained_variance","bic_selected_components",
                       "density_peaks","delta_bic_vs_unimodal",
                       "stable_observer_subsamples","spatial_signature",
                       "coarse_50km_depletion"]
    ranked=(mode.sort_values(["stable_observer_subsamples",
                              "delta_bic_vs_unimodal"],
                             ascending=False,kind="stable")
            [candidate_columns])
    assert len(ranked)==27
    mode_groups={k:{"n":int(len(g)),
                     "coarse_P2_mean":float(g.coarse_50km_depletion.mean()),
                     "pc1_median_variance":float(g.photographic_pc1_explained_variance.median()),
                     "observed_colour_classes":g.photographic_coarse_label_count.value_counts().sort_index().to_dict()}
                 for k,g in d.groupby("mode_status")}
    shape_crosstab=pd.crosstab(d.mode_status,d.spatial_signature)
    result={
      "schema":"fcp_authenticated_modes_photo_spatial_diagnostic_v1",
      "source_run_id":39923312111,
      "n_original_source_photo":447500,
      "n_original_classifiable":227107,
      "n_original_p2_strict_taxa":1761,
      "mode_counts":EXPECTED_MODE,"shape_counts":EXPECTED_GEO,
      "mode_by_shape":{str(c):{str(r):int(v) for r,v in shape_crosstab[c].items()}
                       for c in shape_crosstab.columns},
      "mode_groups":mode_groups,
      "mode_component_number_in_27":{str(k):int(v) for k,v in
                                    mode.bic_selected_components.value_counts().sort_index().items()},
      "mode_detected_peak_number_in_27":{str(k):int(v) for k,v in
                                    mode.density_peaks.value_counts().sort_index().items()},
      "geo_block_cline_R2_relative_month":num(transfer.spatial_block_cline_gain_vs_month),
      "observer_cline_R2_relative_month":num(observer.observer_cline_gain_vs_month),
      "observer_nonlin_increment":num(observer.observer_nonlin_increment_vs_cline),
      "observer_nonlin_increment_in_geo_block_eligible":num(transfer.observer_nonlin_increment_vs_cline),
      "shared_three_shape_candidates":d.loc[d.spatial_signature.isin(
          ["transferable_broad_cline_candidate",
           "local_nonlinearity_only_candidate",
           "broad_trend_plus_local_nonlinearity_candidate"]),
          ["species","mode_status","spatial_signature"]].to_dict(orient="records"),
      "top_mode_taxa":"see top_candidates.csv; requires blind fresh-petal review",
      "note":"statistical model signatures only: mode K is not genetic morph K, and CV success/failure is not identified geographic population history",
    }
    return result,ranked




def adjust_mode_vs_spatial(mode_table, continuous_table, reps=499):
    """Post hoc association: photo-PC1 peaks vs independent continuous JSD P2."""
    needed = {"taxon_id", "global_jsd", "pairs_50km",
              "excess_jsd_50km_month", "n_classified_original",
              "coarse_50km_delta"}
    assert needed.issubset(continuous_table.columns), sorted(needed-set(continuous_table))
    assert continuous_table.taxon_id.is_unique and len(continuous_table)==1761
    keep=continuous_table[list(needed)].copy().rename(columns={"taxon_id":"inat_taxon_id"})
    d=mode_table.merge(keep,on="inat_taxon_id",how="left",validate="one_to_one",
                       indicator=True)
    assert d._merge.eq("both").all() and len(d)==1761
    assert np.max(np.abs(d.coarse_50km_depletion-d.coarse_50km_delta))<1e-10
    assert d.global_jsd.notna().all() and d.excess_jsd_50km_month.notna().all()
    admitted=d.loc[d.mode_status.isin(
        ["recurrent_photographic_mode_candidate",
         "not_resolved_as_multimodal"])].copy()
    assert len(admitted)==1718
    admitted["stable_mode"]=admitted.mode_status.eq("recurrent_photographic_mode_candidate").astype(int)
    def vector(sensitive=False):
        M=[np.ones(len(admitted)),admitted.stable_mode.to_numpy(float),
           admitted.global_jsd.to_numpy(float),
           np.log(admitted.n_classified_original.to_numpy(float)),
           np.log1p(admitted.pairs_50km.to_numpy(float))]
        if sensitive:
            M.extend([admitted.photographic_pc1_explained_variance.to_numpy(float),
                     admitted.photographic_coarse_label_count.to_numpy(float)])
        return np.column_stack(M)
    y=admitted.excess_jsd_50km_month.to_numpy(float)
    names=["main","additional_pc1_variance_and_category_count"]
    estimates={}
    genus_ids=admitted.genus.astype(str).to_numpy()
    units=[np.flatnonzero(genus_ids==g) for g in np.unique(genus_ids)]
    rng=np.random.default_rng(2026102322)
    boot_draws=rng.integers(0,len(units),size=(reps,len(units)))
    for name,sensitive in zip(names,[False,True]):
        X=vector(sensitive)
        assert np.linalg.matrix_rank(X)==X.shape[1]
        coefficient=float(np.linalg.lstsq(X,y,rcond=None)[0][1])
        vals=[]
        for draw in boot_draws:
            idx=np.concatenate([units[k] for k in draw])
            x=X[idx]
            if np.linalg.matrix_rank(x)==X.shape[1]:
                vals.append(float(np.linalg.lstsq(x,y[idx],rcond=None)[0][1]))
        arr=np.array(vals)
        estimates[name]={
            "coefficient_mode_minus_not_detected":coefficient,
            "bootstrap_ci95":[float(v) for v in np.quantile(arr,[.025,.975])] if len(arr)>=450 else None,
            "bootstrap_replicates_successful":len(arr),
            "bootstrap_replicates_requested":reps,
        }
    raw=admitted.groupby("stable_mode").excess_jsd_50km_month.mean()
    return {
        "schema":"photo-mode_vs_50km_continuous_geography_adjusted_v1",
        "n_stable_mode":int(admitted.stable_mode.sum()),
        "n_not_resolved_as_multimodal":int((admitted.stable_mode==0).sum()),
        "n_excluded_unstable_candidates":43,
        "n_distinct_genera":len(units),
        "raw_mean_mode":float(raw[1]),
        "raw_mean_not_detected":float(raw[0]),
        "raw_mode_minus_not_detected":float(raw[1]-raw[0]),
        "models":estimates,
        "adjustment_main":["species-wide JSD","log classified-photo count",
                           "log1p near-pair count"],
        "sensitivity_additions":["PC1 variance fraction","coarse category richness"],
        "not_genetic_morphs_or_fitness":True,
        "posthoc":True,
    }


def self_test():
    bad=pd.DataFrame({"inat_taxon_id":[1,1]})
    try:
        audit(bad)
    except AssertionError:
        print("FCP_MODE_AUDIT_SYNTHETIC_GATES_PASS")
        return
    raise RuntimeError("mode audit accepted invalid input")


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--artifact",type=Path)
    p.add_argument("--out",type=Path)
    p.add_argument("--continuous",type=Path)
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    if a.self_test:
        self_test();return
    assert a.artifact and a.out
    x=pd.read_csv(a.artifact,low_memory=False)
    report,candidates=audit(x)
    if a.continuous is not None:
        extra=pd.read_csv(a.continuous,low_memory=False)
        report["mode_vs_continuous_spatial_adjusted"]=adjust_mode_vs_spatial(x,extra)
    a.out.mkdir(parents=True,exist_ok=True)
    (a.out/"mode_geo_source_audit.json").write_text(json.dumps(
        report,sort_keys=True,indent=2,allow_nan=False)+"\n")
    candidates.to_csv(a.out/"top_mode_candidates.csv",index=False)
    print("FCP_MODE_GEO_INDEPENDENT_ARTIFACT_AUDIT_PASS",json.dumps(report,sort_keys=True))


if __name__=="__main__":
    main()
