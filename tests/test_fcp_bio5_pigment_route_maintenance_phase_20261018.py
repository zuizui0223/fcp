"""Rare-morph invasion and static FCP photo-frequency identifiability boundaries."""
import pytest
from scripts.analysis import fcp_bio5_pigment_route_maintenance_phase_20261018 as f


def test_invasion_requires_both_rare_genetically_inherited_morphs():
    d=f.phase(s=0.0,c_route=0.0,k=.2)
    assert d["local_outcome"]=="protected_intrapopulation_polymorphism"
    assert d["stable_p_true_genetic_morph_A"]==pytest.approx(.5)
    assert d["rare_A_invasion_growth"]>0 and d["rare_B_invasion_growth"]>0
    assert f.phase(s=.23,c_route=0.0,k=.2)["stable_p_true_genetic_morph_A"]==1
    assert f.phase(s=-.25,c_route=0.0,k=.2)["stable_p_true_genetic_morph_A"]==0
    assert f.phase(s=.2,c_route=0.0,k=.2)["local_outcome"]=="locally_A_fixed_or_dominant"
    with pytest.raises(ValueError,match="k>0"):
        f.phase(s=0,c_route=0,k=0)


def test_any_same_true_morph_frequency_can_have_many_mechanisms():
    p,bar_p=.35,.55
    for k,m in ((.1,0),(.2,.01),(.6,.1)):
        net=f.same_static_frequency_multiple_mechanisms(p,k,m,bar_p)
        derivative=p*(1-p)*(net+k*(1-2*p))+m*(bar_p-p)
        assert derivative==pytest.approx(0.,abs=1e-12)
    with pytest.raises(ValueError,match="invalid"):
        f.same_static_frequency_multiple_mechanisms(1.,.2,.1,.5)


def test_tissue_cost_shifts_climate_maintenance_window_but_not_a_universal_gradient():
    d=f.workflow()
    assert d["toy_parameter_values_not_fitted"] is True
    assert d["genetic_morph_A_not_necessarily_white"] is True
    assert d["not_a_reproductive_isolation_model"] is True
    petal=d["route_scenarios"][0]["results"]
    systemic=d["route_scenarios"][1]["results"]
    assert petal[1]["local_outcome"]=="protected_intrapopulation_polymorphism" # 20 C
    assert systemic[1]["local_outcome"]=="locally_B_fixed_or_dominant"
    assert petal[2]["stable_p_true_genetic_morph_A"]>systemic[2]["stable_p_true_genetic_morph_A"]
