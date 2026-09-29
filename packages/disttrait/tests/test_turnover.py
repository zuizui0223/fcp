import numpy as np
import pytest

from disttrait import interaction_beta_partition, null_centered_turnover_rho, turnover_rho


def test_turnover_rho_common_direction():
    separation=[1,2,3,4]
    assert turnover_rho(separation,[0.1,0.2,0.5,0.9]) == pytest.approx(1.0)
    assert turnover_rho(separation,[0.9,0.5,0.2,0.1]) == pytest.approx(-1.0)
    assert turnover_rho(separation,[0.4,0.4,0.4,0.4]) == 0.0


def test_interaction_beta_all_rewiring():
    a=np.array([[1,0],[1,1]])
    b=np.array([[1,1],[0,1]])
    r=interaction_beta_partition(a,b)
    assert r.beta_wn == pytest.approx(1/3)
    assert r.beta_os == pytest.approx(1/3)
    assert r.beta_st == pytest.approx(0.0)
    assert r.beta_wn == pytest.approx(r.beta_st+r.beta_os)


def test_interaction_beta_all_species_turnover():
    a=np.array([[1,0],[0,0]])
    b=np.array([[0,0],[0,1]])
    r=interaction_beta_partition(a,b)
    assert r.beta_wn == pytest.approx(1.0)
    assert r.beta_os == pytest.approx(0.0)
    assert r.beta_st == pytest.approx(1.0)


def test_interaction_beta_mixed_additive_components():
    a=np.array([[1,1],[0,0]])
    b=np.array([[1,0],[0,1]])
    r=interaction_beta_partition(a,b)
    assert r.beta_wn == pytest.approx(0.5)
    assert r.beta_os == pytest.approx(0.25)
    assert r.beta_st == pytest.approx(0.25)
    assert r.beta_wn == pytest.approx(r.beta_st+r.beta_os)


def test_interaction_beta_requires_aligned_species_universe():
    with pytest.raises(ValueError):
        interaction_beta_partition(np.ones((2,2)),np.ones((3,2)))


def test_null_centered_turnover_uses_effect_not_null_zscore():
    separation=[1,2,3,4]
    observed=[0.1,0.2,0.5,0.9]
    null=np.array([
        [0.1,0.5,0.2,0.9],
        [0.9,0.2,0.5,0.1],
        [0.2,0.9,0.1,0.5],
    ])
    r=null_centered_turnover_rho(separation,observed,null)
    assert r.observed_rho == pytest.approx(1.0)
    assert r.delta_rho == pytest.approx(r.observed_rho-r.null_mean_rho)
    assert r.null_rho.shape == (3,)
    assert 0 < r.p_upper <= 1
