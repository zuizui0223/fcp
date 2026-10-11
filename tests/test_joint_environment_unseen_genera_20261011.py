from pathlib import Path
import sys
import numpy as np
import pandas as pd
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/analysis'))
import fcp_joint_environment_unseen_genera_20261011 as m

def source_frame(n_gen=10,n_species=6,n_photo=3):
    rng=np.random.default_rng(75); rows=[]
    for g in range(n_gen):
        for s in range(n_species):
            for p in range(n_photo):
                x=rng.normal(size=4)
                y=np.array([.25+.05*x[0],.25-.05*x[0],.25,.25])
                rows.append(dict(genus=f'G{g}',inat_taxon_id=g*100+s,photo_id=len(rows),n_photos=1,
                    a=x[0],b=x[1],c=x[2],d=x[3],**{f'hue_{h}':float(y[i]) for i,h in enumerate(m.CLASSES)}))
    return pd.DataFrame(rows)

def test_genus_fold_has_no_taxon_or_photo_leak():
    d=source_frame(); folds=m.make_folds(d)
    seen=[]
    for tr,te in folds:
        assert not set(d.iloc[tr].genus)&set(d.iloc[te].genus)
        assert not set(d.iloc[tr].inat_taxon_id)&set(d.iloc[te].inat_taxon_id)
        seen.extend(te)
    assert sorted(seen)==list(range(len(d)))

def test_centering_exact_and_does_not_change_input():
    x=np.array([[1.,3.],[5.,7.],[2.,9.]])
    before=x.copy(); z=m.center(x,np.array(['A','A','B']))
    assert np.allclose(z,[[-2,-2],[2,2],[0,0]])
    assert np.array_equal(x,before)

def test_predictions_never_read_heldout_colour():
    d=source_frame(); tr=d[d.genus!='G0']; te=d[d.genus=='G0'].copy()
    first=m.predict(tr,te,('a','b','c','d'),'within_genus')
    for h in m.CLASSES:te['hue_'+h]=np.arange(len(te))*500
    second=m.predict(tr,te,('a','b','c','d'),'within_genus')
    assert np.allclose(first,second)

def test_training_unknown_genus_guard_rejects_leak():
    d=source_frame()
    with pytest.raises(ValueError,match='genus leakage'):
        m.predict(d,d.iloc[:6],('a','b'),'absolute')

def test_within_genus_fit_invariant_to_genus_intercept():
    d=source_frame(); tr=d[d.genus!='G0'].copy(); te=d[d.genus=='G0'].copy()
    one=m.predict(tr,te,('a','b','c','d'),'within_genus')
    intercept=tr.genus.str[1:].astype(int)*.01
    tr['hue_white']+=intercept;tr['hue_yellow_orange']-=intercept
    two=m.predict(tr,te,('a','b','c','d'),'within_genus')
    assert np.allclose(one,two)

def test_within_species_fit_invariant_to_species_mean():
    d=source_frame();tr=d[d.genus!='G0'].copy();te=d[d.genus=='G0']
    one=m.predict(tr,te,('a','b','c','d'),'within_species')
    shift=(tr.inat_taxon_id%7)*.01
    tr['hue_white']+=shift;tr['hue_red_pink']-=shift
    two=m.predict(tr,te,('a','b','c','d'),'within_species')
    assert np.allclose(one,two)

def test_absolute_predictions_are_probabilities():
    d=source_frame(); p=m.predict(d[d.genus!='G0'],d[d.genus=='G0'],('a','b'),'absolute')
    assert np.all(p>=0);assert np.allclose(p.sum(1),1)

def test_orthogonal_scores_sum_to_full_loss():
    rng=np.random.default_rng(1);r=rng.normal(size=(50,4));r-=r.mean(1)[:,None]
    loss=m.loss_components(r)
    assert np.allclose(loss[:,0],loss[:,1]+loss[:,2]); assert (loss>=-1e-12).all()

def test_equal_group_bootstrap_mathematical_constant():
    base=np.tile([.4,.2,.2],(60,1));full=base*.9
    out=m.summarize(base,full,np.repeat([f'G{i}' for i in range(10)],6),reps=99)
    assert np.isclose(out['genus_equal']['total']['gain'],.04)
    assert np.allclose(out['genus_equal']['total']['relative_gain_percent_CI95'],[10,10])

def test_full_environment_signal_recovers_across_unseen_genera():
    d=source_frame(n_gen=15,n_species=10)
    b=[];f=[]
    for tr,te in m.make_folds(d):
        target=d.iloc[te]; y=target[[f'hue_{h}' for h in m.CLASSES]].to_numpy()
        geo=m.predict(d.iloc[tr],target,('b','c','d'),'absolute')
        env=m.predict(d.iloc[tr],target,('a','b','c','d'),'absolute')
        b.append(((y-geo)**2).sum());f.append(((y-env)**2).sum())
    assert sum(f)<.05*sum(b)

def test_original_archive_hash_must_match(tmp_path):
    p=tmp_path/'bad';p.write_bytes(b'wrong')
    with pytest.raises(ValueError,match='SHA256'):m.verify_file(p,'0'*64)

def test_mean_deviation_decomposition_exact_with_unequal_genus_sizes():
    rng=np.random.default_rng(87)
    groups=np.repeat(['a','b','c'],[3,7,2])
    residual=rng.normal(size=(12,4));residual-=residual.mean(1)[:,None]
    loss=m.genus_error_partition(residual,groups)
    for g in set(groups):
        v=loss[groups==g].mean(0)
        assert np.isclose(v[0],v[1]+v[2])
    assert np.isclose(loss[:,0].mean(),loss[:,1:].sum(1).mean())

def test_general_summary_labels_for_genus_composition_decomposition():
    b=np.tile([.5,.3,.2],(12,1));f=b*.95
    result=m.summarize(b,f,np.repeat(['a','b','c'],4),reps=99,
                       components=('total','genus_mean','within_genus_deviation'))
    assert set(result['genus_equal'])=={'total','genus_mean','within_genus_deviation'}
