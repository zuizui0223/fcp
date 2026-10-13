import sys
from pathlib import Path
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from field_partition import (load_sisyrinchium, load_boechera, within_site,
    permute_conditional, select_strata, summarize_species, expected_mixed,
    hold_reason)


def test_hypergeom_two_equal_quadrats_expected_one_third():
    n = np.array([2, 2]); v = np.array([0, 2])
    a = within_site(n, v)
    assert a['B'] == pytest.approx(1.0)
    assert a['expected_B'] == pytest.approx(1.0/3.0)
    assert a['delta'] == pytest.approx(2.0/3.0)
    assert a['observed_mixed'] == 0
    assert expected_mixed(n, 2).sum() == pytest.approx(4.0/3.0)


def test_hypergeom_fixed_inventory_and_seed():
    n=np.array([2, 1, 3]); v=np.array([0, 1, 2])
    a=permute_conditional(n,v,1000,99)
    b=permute_conditional(n,v,1000,99)
    assert np.array_equal(a['B'],b['B'])
    assert len(a['B'])==1000
    assert a['B'].min()>=0
    assert a['B'].mean() == pytest.approx(within_site(n,v)['expected_B'],abs=0.08)


def test_zeros_not_white_and_monochrome_hold():
    a=within_site(np.array([0, 2, 2]),np.array([0, 1, 1]))
    assert a['n_occupied']==2
    assert a['n_zero_target']==1
    assert hold_reason(20,3,0,20)=='HOLD_MONOMORPHIC'
    assert hold_reason(20,3,5,15) is None
    assert hold_reason(19,2,5,14)=='HOLD_LOW_N'


def test_sisyrinchium_true_m_species_only():
    x=load_sisyrinchium(ROOT/'source/Sisyrinchium_data.zip')
    assert len(x)==512
    assert x.duplicated(['year','population','line','quad']).sum()==0
    assert (x[['v','w','l','s']]>=0).all().all()
    assert (x['n']==x['v']+x['w']).all()
    assert set(x['year'])=={'Y11','Y13'}
    assert x.groupby(['year','population']).ngroups == 16


def test_boechera_individual_plot_exact_identity():
    x,y,qa=load_boechera(ROOT/'source/Boechera_natural_populations_individual.txt',ROOT/'source/Boechera_natural_populations_plot.txt')
    assert len(x)==547 and len(y)==26
    assert qa['joined_plantid_typos']==3
    assert qa['all_26_plots_match_original_individual_counts']
    assert x['plantID'].nunique()==547
    assert y['n'].sum()==547 and y['v'].sum()==122
    assert y['Site_name'].nunique()==9


def test_permutation_morph_inventory_and_site_weights():
    # A fixed field site gets exactly one vote despite larger population.
    a=[{'site':'small','delta':0.3},{'site':'big','delta':0.1}]
    assert summarize_species(a)['mean_delta']==pytest.approx(0.2)


def test_latest_year_sensitivity_does_not_double_count_site():
    x=load_sisyrinchium(ROOT/'source/Sisyrinchium_data.zip')
    m=select_strata(x,'sisyrinchium',latest_year_only=True)
    names=[q['population'] for q in m]
    assert len(names)==len(set(names))


def test_mixed_plot_probability_never_unrealistic():
    for n,V in [([1,1],1),([5,7,10],3),([2,3,4],6)]:
        p=expected_mixed(np.array(n),V)
        assert np.all((p>=-1e-12)&(p<=1+1e-12))
        assert np.allclose(p[np.array(n)==1],0)


def test_quality_receipt_is_json_serializable():
    import json
    _, _, qa = load_boechera(ROOT/'source/Boechera_natural_populations_individual.txt',ROOT/'source/Boechera_natural_populations_plot.txt')
    assert isinstance(qa['all_26_plots_match_original_individual_counts'],bool)
    json.dumps(qa)
