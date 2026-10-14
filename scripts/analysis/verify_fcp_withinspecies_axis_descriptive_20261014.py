#!/usr/bin/env python3
"""FCP-only source rederivation of the 3-cohort photographed white-axis summary.

Source retrospective discovery/validation and prospective third cohort remain
separate tests. This script performs only arithmetic, no joint permutation/p.
"""
from __future__ import annotations
import json
import math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
A=ROOT/"results/polymorphism_white_axis_targeted_test_20260912/result.json"
B=ROOT/"results/polymorphism_h2_third_cohort_prospective_white_axis_20260917/result.json"
OUT=ROOT/"results/polymorphism_white_axis_three_cohorts_descriptive_20261014.json"

def verify():
    old=json.loads(A.read_text(encoding="utf8"))
    new=json.loads(B.read_text(encoding="utf8"))
    manifest=json.loads(OUT.read_text(encoding="utf8"))
    assert old["status"]=="targeted_post_audit_decomposition_not_untouched_confirmatory_test"
    assert new["status"]=="untouched_prospective_test_of_previously_frozen_axis"
    assert old["fixed_axis_palette_loadings"]==new["fixed_axis_palette_loadings"]
    for code,name,expected in [("primary_0_10","primary",439),("strict_0_20","strict",226)]:
        def old_row(key):
            x=old["thresholds"][code][key]
            return (x["species"],x["observed_mean_squared_white_axis_alignment"],x["structured_null_summary"]["mean"])
        def third():
            x=new["thresholds"][code]
            return (x["species"],x["observed_W"],x["structured_null_summary"]["mean"])
        rows=[old_row("discovery"),old_row("reserve"),third()]
        total=sum(row[0] for row in rows)
        W=sum(x*y for x,y,_ in rows)/total
        null=sum(x*n for x,_,n in rows)/total
        receipt=next(d for d in manifest["cohorts"] if d["tier"]==name)
        assert total==expected==receipt["n_eligible"]
        assert math.isclose(W,receipt["observed_species_equal_mean_W"],abs_tol=1e-13)
        assert math.isclose(null,receipt["weighted_original_separate_cohort_null_means"],abs_tol=1e-13)
        assert receipt["pooled_p_value"] is None
        assert len(receipt["cohorts"])==3
        assert receipt["cohorts"][-1]["source"]=="prospective"
    assert manifest["research_boundary"]["main_independent_prospective_H2_n"]==158
    assert manifest["research_boundary"]["not_global_species_photo_colour_moisture_analysis"] is True
    return manifest

if __name__=="__main__":
    r=verify()
    print("PASS fcp within-species H2: primary",r["cohorts"][0]["n_eligible"],
          "strict",r["cohorts"][1]["n_eligible"],"prospective source H2 only")
