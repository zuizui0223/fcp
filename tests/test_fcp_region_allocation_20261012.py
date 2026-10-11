"""Guard tests of the 2026-10-12 fixed species-sector calculation."""
import numpy as np
import pandas as pd
import pytest
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import fcp_region_allocation_20261012 as m


def toy(n=2, mixed=True):
    # 2-3 sectors, true fractions and model forecasts all legal four-state shares.
    probs=np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0]],float)[:n]
    if not mixed:probs[:]=[1,0,0,0]
    out=[]
    for j in range(n):
        a={'inat_taxon_id':1,'genus':'Example','scheme':'buffer500','fold':0,
            'sector':j,'n_photos':3,'environment_support':True}
        for k,c in enumerate(m.CLASSES):
            a['truth__'+c]=probs[j,k]
            a['ridge_geo__p__'+c]=.25
            a['ridge_full__p__'+c]=probs[j,k]
            a['nonlinear_geo__p__'+c]=.25
            a['nonlinear_full__p__'+c]=probs[j,k]
        out.append(a)
    return pd.DataFrame(out)


def test_gate_no_outcome_dependency():
    d=toy(3)
    a=m.opportunity_gate(d,2)[['inat_taxon_id','sector']]
    for c in ('truth__white','ridge_full__p__red_pink'):
        d[c]=100
    b=m.opportunity_gate(d,2)[['inat_taxon_id','sector']]
    pd.testing.assert_frame_equal(a,b)
    d.loc[d.sector.eq(0),'n_photos']=1
    selected=m.opportunity_gate(d,2)
    assert len(selected)==2 and 0 not in selected.sector.tolist()


def test_centred_full_predictor_hits_exact_true_shifts():
    d=toy(3)
    reg=m.all_model_losses(d)
    s=m.species_ledger(reg)
    assert s['nonlinear_full__total'].iloc[0]==pytest.approx(0,abs=1e-14)
    assert s['zero__total'].iloc[0]==pytest.approx(2/3,abs=1e-14)
    assert s['environment_nonlinear__total'].iloc[0]==pytest.approx(2/3,abs=1e-14)
    assert s['environment_ridge__total'].iloc[0]==pytest.approx(2/3,abs=1e-14)
    assert s['full_nonlinear_vs_ridge__total'].iloc[0]==pytest.approx(0,abs=1e-14)


def test_constant_colour_true_regional_difference_is_zero_not_missing():
    d=toy(2,False)
    reg=m.all_model_losses(d)
    s=m.species_ledger(reg)
    assert s['zero__total'].iloc[0]==pytest.approx(0,abs=1e-14)
    assert s['nonlinear_full__total'].iloc[0]==pytest.approx(0,abs=1e-14)
    assert s['environment_nonlinear__total'].iloc[0]==pytest.approx(0,abs=1e-14)


def test_sector_reordering_is_invariant():
    original=toy(3)
    a=m.species_ledger(m.all_model_losses(original))
    b=m.species_ledger(m.all_model_losses(original.sample(frac=1,random_state=55)))
    for k in m.CONTRASTS:
        assert a[k+'__total'].iloc[0]==pytest.approx(b[k+'__total'].iloc[0])


def test_white_axis_orthogonal_partition():
    d=toy(3)
    a=m.species_ledger(m.all_model_losses(d))
    for k in ('zero',*m.MODELS,*m.CONTRASTS):
        assert a[k+'__total'].iloc[0]==pytest.approx(a[k+'__white_axis'].iloc[0]+a[k+'__remaining_hue'].iloc[0])


def test_summary_holds_inadequate_sample_and_refuses_fake_significance():
    d=toy(3)
    a=m.summarize_ledger(m.species_ledger(m.all_model_losses(d)),10,with_ci=False)
    assert a['status']=='HOLD_LOW_SUPPORT'
    assert a['bootstrap_reps']==0
    assert 'CI95' not in a['genus_equal']['environment_nonlinear__total']


def test_exact_primary_null_requires_true_multispecies_support():
    d=toy(3)
    r=m.all_model_losses(d)
    led=m.species_ledger(r)
    with pytest.raises(ValueError,match='insufficient primary support'):
        m.photo_sector_null(r,led)


def test_source_sha_length_and_model_keys_locked():
    assert len(m.SOURCE_SHA)==64
    assert m.CONTRASTS['environment_nonlinear']==('nonlinear_geo','nonlinear_full')
    assert len(m.CLASSES)==4 and len(m.MODELS)==4 and len(m.SCHEMES)==3


def test_observed_inventory_count_roundtrip():
    d=toy(3)
    p=d[['truth__'+c for c in m.CLASSES]].to_numpy()
    counts=np.rint(p*d.n_photos.to_numpy()[:,None]).astype(int)
    assert (counts.sum(1)==d.n_photos.to_numpy()).all()
    assert (counts.sum(0)==[3,3,3,0]).all()