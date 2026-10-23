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
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    if a.self_test:
        self_test();return
    assert a.artifact and a.out
    x=pd.read_csv(a.artifact,low_memory=False)
    report,candidates=audit(x)
    a.out.mkdir(parents=True,exist_ok=True)
    (a.out/"mode_geo_source_audit.json").write_text(json.dumps(
        report,sort_keys=True,indent=2,allow_nan=False)+"\n")
    candidates.to_csv(a.out/"top_mode_candidates.csv",index=False)
    print("FCP_MODE_GEO_INDEPENDENT_ARTIFACT_AUDIT_PASS",json.dumps(report,sort_keys=True))


if __name__=="__main__":
    main()
