"""Posthoc derived hue-pair tests; original FCP source inference remains unchanged."""
from __future__ import annotations
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/"scripts/analysis/audit_fcp_same_genus_east_elsewhere_rain_tradeoff_20261011.py"
sys.path.insert(0,str(SRC.parent))
spec=importlib.util.spec_from_file_location("fcp_white_blue_posthoc_source",SRC)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def synthetic(n=43,seed=17):
    rng=np.random.default_rng(seed)
    east={}; other={}
    for i in range(n):
        X=rng.normal(size=(16,8));X-=X.mean(axis=0)
        M=X.T@X
        beta=np.zeros((8,4))
        z=np.array([.06,-.035,-.01,-.015])+rng.normal(0,.018,4)
        z-=z.mean()
        beta[0,:]=z
        east[f"G{i}"]=(M,M@beta,16)
        other[f"G{i}"]=(M,np.zeros((8,4)),16)
    return east,other


def test_historical_default_is_exactly_unchanged():
    a,b=synthetic()
    parent=mod.paired_genus_bootstrap(a,b,repetitions=111,seed=99)
    extra=mod.paired_genus_bootstrap(a,b,repetitions=111,seed=99,
                                     include_white_blue_contrast=True)
    assert "posthoc_white_blue_region_rain_contrast" not in parent
    assert parent=={k:v for k,v in extra.items()
                   if k!="posthoc_white_blue_region_rain_contrast"}


def test_joint_contrast_is_identical_draw_coupled_and_exact_probability_units():
    a,b=synthetic()
    value=mod.paired_genus_bootstrap(a,b,repetitions=111,seed=99,
                                     include_white_blue_contrast=True)
    hue=value["response_four_original_photo_hues"]
    q=value["posthoc_white_blue_region_rain_contrast"]
    assert q["point"]==pytest.approx(
        hue["white"]["east_minus_elsewhere_rainfall_slope"]
        -hue["blue_purple"]["east_minus_elsewhere_rainfall_slope"],abs=1e-12)
    assert q["east_white_minus_blue"]>0
    assert q["unit"]=="photo_colour_probability_per_global_BIO12_SD"
    assert q["six_hue_pair_maxT_critical_95"]>=1
    assert q["six_pair_simultaneous_bootstrap_approx_CI95"][0] <= q["point"]
    assert q["six_pair_simultaneous_bootstrap_approx_CI95"][1] >= q["point"]
    assert q["bonferroni_conservative_one_sided_95pct_lower"] <= q["paired_genus_bootstrap_percentile_CI95"][0]+1e-12


def test_minimum_support_status_unaltered():
    a,b=synthetic(n=29)
    result=mod.paired_genus_bootstrap(a,b,repetitions=39,seed=99,
                                      include_white_blue_contrast=True)
    assert result["status"]=="HOLD_FEWER_THAN_30_SHARED_GENERA"
    assert "posthoc_white_blue_region_rain_contrast" not in result


def test_all_six_comparisons_are_not_the_same_fitted_variable():
    a,b=synthetic()
    out=mod.paired_genus_bootstrap(a,b,repetitions=111,seed=99,
                                    include_white_blue_contrast=True)
    assert len(mod.HUES)==4
    assert out["posthoc_white_blue_region_rain_contrast"]["retrospective_selected_east_sector_and_colour_pair"]
    assert out["posthoc_white_blue_region_rain_contrast"]["six_pair_family_only_region_not_multiplicity_adjusted"]
    assert sum(out["response_four_original_photo_hues"][h]["east_rainfall_slope"]
               for h in mod.HUES)==pytest.approx(0,abs=1e-10)
