"""Test photo-only locality signals do not imply genetic polymorphism or speciation."""
import numpy as np
import pandas as pd
import pytest
from scripts.analysis import fcp_local_coexistence_feasibility_20261018 as f


def sample():
    morphs=["white","blue_purple","white","blue_purple","white","blue_purple","blue_purple"]
    return pd.DataFrame({
        "inat_taxon_id":[4]*7,"species":["Campanula synthetic"]*7,
        "photo_id":np.arange(200,207),"observation_id":np.arange(300,307),
        "observer_id":np.arange(20,27),
        "latitude":[35,35.001,35.002,35.003,35.004,35.005,36],
        "longitude":[135]*7,"positional_accuracy_m":[12]*7,
        "observed_on":["2026-05-10"]*6+["2026-06-10"],
        "morph":morphs,"global_classifiable":[True]*7,
    })


def test_distinct_observation_and_observer_month_pairs():
    d=f.classify(sample())
    u,v,dist=f.great_circle_pairs(d)
    q=f.strict_pairs(d,u,v,dist,spatial_km=5,date_unit="year_month",max_accuracy_m=1000)
    assert q.sum()==15
    assert (d.morph.to_numpy()[u[q]]!=d.morph.to_numpy()[v[q]]).sum()==9
    g=d.copy()
    g.loc[1,"observer_id"]=g.loc[0,"observer_id"]
    u,v,dist=f.great_circle_pairs(g)
    assert f.strict_pairs(g,u,v,dist,spatial_km=5,date_unit="year_month",max_accuracy_m=1000).sum()==14
    g=d.copy()
    g.loc[1,"positional_accuracy_m"]=3000.
    u,v,dist=f.great_circle_pairs(g)
    assert f.strict_pairs(g,u,v,dist,spatial_km=5,date_unit="year_month",max_accuracy_m=1000).sum()==10


def test_within_month_permutations_keep_colour_composition():
    d=f.classify(sample())
    u,v,km=f.great_circle_pairs(d)
    good=f.strict_pairs(d,u,v,km,spatial_km=5,date_unit="year_month",max_accuracy_m=1000)
    observed,z=f.conditional_month_null(d,u[good],v[good])
    assert observed==pytest.approx(9/15)
    assert z.shape==(199,)
    assert np.allclose(z,observed)
    assert np.array_equal(z,f.conditional_month_null(d,u[good],v[good])[1])


def test_year_month_and_position_gates_reduce_support():
    d=sample()
    d.loc[1,"observed_on"]="2026-07-10"
    d.loc[2,"positional_accuracy_m"]=5000
    d.loc[3,"observer_id"]=d.loc[4,"observer_id"]
    d=f.classify(d)
    u,v,km=f.great_circle_pairs(d)
    pair=f.strict_pairs(d,u,v,km,spatial_km=5,date_unit="year_month",max_accuracy_m=1000)
    assert pair.sum()<15
    assert not np.any(pair[(u==1)|(v==1)])


def test_posthoc_null_and_source_contract():
    rows=[{"taxon_id":j,"strict_pairs_5km_same_month":5,
           "local_discordance":v} for j,v in ((1,.1),(2,.2),(3,.3))]
    null={j:np.zeros(199) for j in (1,2,3)}
    result=f.summarise_test(rows,null,5)
    assert result["n_evaluable_species"]==3
    assert result["is_confirmatory"] is False
    assert result["upper_p_199_month_conditioned_nulls"]==1.
    assert f.PHOTO_SHA256=="8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6"
    assert f.SOURCE_ROWS==447500 and f.SOURCE_COLOUR_COUNT==227107


def test_no_population_or_genetic_inference():
    assert "genetically distinct morphs" in f.__doc__
    assert "speciation" in f.__doc__
    assert f.MONTH_RADIUS_KM==5 and f.MONTH_MAX_POSITION_ACCURACY_M==1000
