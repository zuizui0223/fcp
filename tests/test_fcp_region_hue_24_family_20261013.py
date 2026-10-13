"""FCP whole-geography x whole-colour-pair family tests (no empirical source needed)."""
import importlib.util
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

HERE=Path(__file__).resolve().parents[1]
SOURCE=HERE/"scripts/analysis/audit_fcp_region_hue_24_family_20261013.py"
sys.path.insert(0,str(SOURCE.parent))
spec=importlib.util.spec_from_file_location("fcp_allregion_allhue_fwer",SOURCE)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def synthetic_regions(n=52,signal=.04,seed=23):
    rng=np.random.default_rng(seed)
    v={}
    for sector in range(4):
        a={}; b={}
        for g in range(n):
            mtx=np.eye(8)*(1 + (g%7)/3)
            effects=rng.normal(0,.025,4)
            effects-=effects.mean()
            if sector==3:
                effects+=signal*np.array([1.5,-.5,-.5,-.5])
            bet=np.zeros((8,4));bet[0,:]=effects
            a[str(g)]=(mtx,mtx@bet,10)
            b[str(g)]=(mtx,np.zeros((8,4)),10)
        v[sector]={"status":"OK","inside":a,"outside":b,
                   "meta":{"sector":sector,"n_genera":n,"n_photos":n*20,
                           "n_unique_named_species":n*10,
                           "n_region_named_species":n*10,
                           "n_other_named_species":n*10}}
    return v


def test_sector_boundaries_unchanged():
    x=pd.DataFrame({"longitude":[-180,-179.99,-90,0,89.999,90,179.999,180]})
    assert list(m.source_sectors(x))==[0,0,1,2,2,3,3,3]


def test_paired_region_and_hue_family_complete_and_deterministic():
    regions=synthetic_regions()
    v=m.simultaneous_source_bootstrap(regions,399,31)
    copy=m.simultaneous_source_bootstrap(regions,399,31)
    assert v==copy
    assert v["status"]=="OK"
    assert v["n_test_family"]==24
    assert v["n_union_genera"]==52
    assert len(v["all_directional_contrasts"])==4*6
    wb=v["historical_white_blue_east"]
    assert wb["point_pp_per_BIO12_SD"]>0
    assert wb["sector"]==3
    assert (wb["hue_A"],wb["hue_B"])==("white","blue_purple")
    assert v["maxT_critical_95"] > 1.96


def test_shared_genus_draws_are_not_independent_across_regions():
    # For identical genus-specific estimates across two regions, the same
    # multinomial weights must reconstruct identical bootstrap CI values.
    regions=synthetic_regions(signal=0.0)
    for sector in range(1,4):
        regions[sector]["inside"]=regions[0]["inside"]
        regions[sector]["outside"]=regions[0]["outside"]
    v=m.simultaneous_source_bootstrap(regions,399,31)
    assert v['status']=='OK'
    by_sector=[v['all_directional_contrasts'][s*6:(s+1)*6] for s in range(4)]
    for j in range(6):
        assert all(by_sector[s][j]['nominal_bootstrap_CI95_pp']==by_sector[0][j]['nominal_bootstrap_CI95_pp'] for s in range(4))
    assert v['maxT_critical_95']>=1.5


def test_exact_zero_variance_directional_pair_fails_closed():
    regions=synthetic_regions()
    for reg in regions.values():
        for side in ('inside','outside'):
            for genus,(mat,y,n) in reg[side].items():
                reg[side][genus]=(mat,np.zeros_like(y),n)
    v=m.simultaneous_source_bootstrap(regions,399,31)
    assert v['status']=='HOLD_DEGENERATE_DIRECTIONAL_PAIR'
    assert v['n_degenerate']==24


def test_insufficient_genera_are_not_reported_as_colour_absence():
    regions=synthetic_regions(n=29)
    out=m.simultaneous_source_bootstrap(regions,399,31)
    assert out['status']=='HOLD_FEWER_THAN_30_UNION_GENERA'


def test_historical_eastern_photo_point_contract_is_arithmetic():
    v=m.simultaneous_source_bootstrap(synthetic_regions(),399,31)
    sr=v['sector_metadata'][3]
    inside=sr['four_hue_inside_slopes_per_sd']; other=sr['four_hue_outside_slopes_per_sd']
    implied=100*((inside['white']-inside['blue_purple'])-(other['white']-other['blue_purple']))
    assert v['historical_white_blue_east']['point_pp_per_BIO12_SD']==pytest.approx(implied,abs=1e-12)
