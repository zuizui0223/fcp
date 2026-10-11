"""Fixed nonlinear transport contract: synthetic data only, no source outcomes."""
from pathlib import Path
import importlib.util
import sys
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts/analysis'))
NAME = 'fcp_nonlinear_environment_transport_20261011'

@pytest.fixture(scope='module')
def m():
    if importlib.util.find_spec(NAME) is None:
        pytest.skip('Extension not yet implemented')
    return __import__(NAME)


def test_fixed_transport_extension_is_available():
    assert importlib.util.find_spec(NAME) is not None, 'Fixed nonlinear transport extension is missing'


def training(m):
    rng = np.random.default_rng(611)
    tr = pd.DataFrame(rng.normal(size=(240, 5)), columns=['a', 'b', *m.GEO])
    tr['inat_taxon_id'] = np.arange(len(tr))
    tr['genus'] = ['Train' + str(i // 10) for i in range(len(tr))]
    y = rng.dirichlet(np.ones(4), len(tr))
    for j, h in enumerate(m.CLASSES):
        tr['hue_' + h] = y[:, j]
    te = pd.DataFrame(rng.normal(size=(7, 5)), columns=['a', 'b', *m.GEO])
    te['inat_taxon_id'] = np.arange(len(te)) + 1000
    te['genus'] = ['Unseen' + str(i) for i in range(len(te))]
    return tr, te


def records(m):
    rng = np.random.default_rng(62)
    rows = []
    for s in range(12):
        for sector in range(1 + s % 3):
            r = {'inat_taxon_id':s, 'genus':'G' + str([0,0,0,1,2,2,3,3,4,5,5,5][s]), 'sector':sector,
                 'fold':[0,0,0,1,2,2,3,3,4,5,5,5][s] % 5, 'scheme':'buffer500', 'n_photos':1 + 10 * sector}
            truth = rng.dirichlet(np.ones(4))
            for j, h in enumerate(m.CLASSES):r['truth__' + h] = truth[j]
            for model in m.MODELS:
                pred = rng.dirichlet(np.ones(4))
                loss = m.loss_components((truth - pred)[None])[0]
                for j,h in enumerate(m.CLASSES):r[model + '__p__' + h] = pred[j]
                for j,c in enumerate(m.COMPONENTS):r[model + '__' + c] = loss[j]
            rows.append(r)
    return pd.DataFrame(rows)


def test_fixed_tree_fits_are_deterministic_simplex_and_ignore_heldout_outcomes(m):
    tr, te = training(m)
    first = m.predict_nonlinear(tr, te, ('a', 'b'))
    for h in m.CLASSES:te['hue_' + h] = 1e9
    second = m.predict_nonlinear(tr, te, ('a', 'b'))
    assert np.allclose(first, second, atol=1e-12, rtol=0)
    assert np.isfinite(first).all() and (first >= 0).all()
    assert np.allclose(first.sum(1), 1, atol=1e-12)
    # Each tree leaf averages training compositions; predictions stay in their range.
    assert np.all(first >= tr[['hue_' + h for h in m.CLASSES]].min().to_numpy() - 1e-12)
    assert m.TREE_PARAMS == dict(n_estimators=128, max_depth=16, min_samples_leaf=20,
        max_features=1.0, bootstrap=False, criterion='squared_error', random_state=2026101161, n_jobs=2)


def test_four_model_splits_reject_genus_species_leakage_and_duplicate_training(m):
    tr,te = training(m)
    pred = m.predict_models(tr,te,('a','b'))
    assert set(pred) == set(m.MODELS)
    for p in pred.values():assert p.shape == (7,4) and np.allclose(p.sum(1),1)
    bad = te.copy();bad.loc[0,'genus'] = tr.genus.iloc[0]
    with pytest.raises(ValueError, match='genus'):m.predict_models(tr,bad,('a','b'))
    bad = te.copy();bad.loc[0,'inat_taxon_id'] = tr.inat_taxon_id.iloc[0]
    with pytest.raises(ValueError, match='species'):m.predict_models(tr,bad,('a','b'))
    with pytest.raises(ValueError, match='unique'):m.predict_models(pd.concat([tr,tr.iloc[:1]]),te,('a','b'))


def test_summary_pairs_all_gains_after_equal_species_sector_averaging(m):
    d = records(m)
    result = m.summarize_records(d, reps=99, seed=18)
    cols = [model+'__'+c for model in m.MODELS for c in m.COMPONENTS]
    sp = d.groupby('inat_taxon_id').agg({**dict.fromkeys(cols,'mean'),'genus':'first'})
    gm = sp.groupby('genus')[cols].mean()
    for weight,values in [('species_equal',sp[cols].mean()),('genus_equal',gm.mean())]:
        s = result[weight]['total']
        losses = {model:values[model+'__total'] for model in m.MODELS}
        for model in m.MODELS:assert s['losses'][model]['estimate'] == pytest.approx(losses[model])
        ng = losses['nonlinear_geo'] - losses['nonlinear_full']
        rg = losses['ridge_geo'] - losses['ridge_full']
        assert s['contrasts']['nonlinear_environment_gain']['estimate'] == pytest.approx(ng)
        assert s['contrasts']['nonlinear_minus_ridge_environment_gain']['estimate'] == pytest.approx(ng-rg)
        assert s['contrasts']['full_model_improvement']['estimate'] == pytest.approx(losses['ridge_full']-losses['nonlinear_full'])
        assert s['contrasts']['geography_model_improvement']['estimate'] == pytest.approx(losses['ridge_geo']-losses['nonlinear_geo'])
    # Whole-genus bootstrap gives unequal-size genus clusters their own species weights.
    genera = sp.groupby('genus',sort=True)
    cnt = genera.size().to_numpy()
    means = genera[cols].mean().to_numpy()
    ix = np.random.default_rng(18).integers(len(cnt),size=(99,len(cnt)))
    draws = means[ix].mean(1)
    a=cols.index('nonlinear_geo__total');b=cols.index('nonlinear_full__total')
    wanted = np.quantile(draws[:,a]-draws[:,b],[.025,.975])
    assert np.allclose(result['genus_equal']['total']['contrasts']['nonlinear_environment_gain']['CI95'], wanted)
    assert m.summarize_records(d.assign(n_photos=999),reps=99,seed=18) == result


def test_reject_duplicate_missing_and_mismatched_comparison_rows(m):
    d=records(m);expected=d[['inat_taxon_id','sector','scheme','fold','genus']]
    m.validate_records(d,expected)
    with pytest.raises(ValueError,match='duplicate'):m.validate_records(pd.concat([d,d.iloc[:1]]),expected)
    with pytest.raises(ValueError,match='target'):m.validate_records(d.iloc[1:],expected)
    with pytest.raises(ValueError,match='target'):m.validate_records(d.assign(genus='changed'),expected)
    with pytest.raises(ValueError,match='finite'):m.validate_records(d.assign(nonlinear_full__total=np.nan),expected)
    with pytest.raises(ValueError,match='loss'):m.validate_records(d.assign(nonlinear_full__total=100),expected)
    with pytest.raises(ValueError,match='duplicate'):m.summarize_records(pd.concat([d,d.iloc[:1]]),99)


def test_archived_ridge_comparison_aligns_keys_and_rejects_missing_rows(m,tmp_path):
    d=records(m)
    parent=d[['inat_taxon_id','genus','sector','fold','scheme']].copy()
    for c in m.COMPONENTS:
        parent['geography__'+c]=d['ridge_geo__'+c]
        parent['environment36__'+c]=d['ridge_full__'+c]
    path=tmp_path/'soil_all_oof_losses.csv.gz';parent.sample(frac=1,random_state=4).to_csv(path,index=False)
    assert m.compare_archived_ridge(d,path)['status']=='MATCH'
    parent.iloc[1:].to_csv(path,index=False)
    with pytest.raises(ValueError,match='target'):m.compare_archived_ridge(d,path)


def test_pipeline_rejects_insufficient_training_instead_of_dropping_targets(m,tmp_path):
    # Use actual species aggregation/folds/exclusions/records; only lower source support
    # is synthetic. No altered model parameters, predictions or scoring helpers.
    rng=np.random.default_rng(71);rows=[]
    for i in range(40):
        for s in range(2):
            row={'inat_taxon_id':i,'genus':'G'+str(i),'photo_id':2*i+s,
                'latitude':float(rng.uniform(-30,30)),'longitude':[-135.,45.][s],
                'morph':m.CLASSES[(i+s)%4],'flower_effective_pixels':30000,
                'nuisance_pixel_fraction':.1}
            row.update({f:float(rng.normal()) for f in [*m.GEO,*m.FULL]})
            rows.append(row)
    raw=pd.DataFrame(rows)
    exclusions,_=m.spatial_exclusions(raw)
    # An empty buffer training set is an explicit failed design, never silently dropped.
    with pytest.raises(ValueError,match='training support'):
        m.evaluate_case(raw,False,exclusions,99,tmp_path)
    # Geometry where each target sector has disjoint species offers all three regimes.
    raw['longitude']=[-135.,-45.,45.,135.]*20
    exclusions,_=m.spatial_exclusions(raw)
    # Main source support thresholds deliberately fail rather than generate partial rows.
    with pytest.raises(ValueError,match='training support'):
        m.evaluate_case(raw,False,exclusions,99,tmp_path)


def test_fixed_nonlinear_learner_encodes_shared_curve_on_unseen_genus(m):
    # Frozen deterministic U shape, with no original-source outcomes or tuning.
    x=np.linspace(-1.,1.,1200)
    tr=pd.DataFrame({'a':x,'inat_taxon_id':np.arange(1200),
                     'genus':['Train'+str(i//40) for i in range(1200)]})
    white=.1+.8*x*x
    for j,h in enumerate(m.CLASSES):tr['hue_'+h]=white if j==0 else (1-white)/3
    xt=np.linspace(-.995,.995,240)
    te=pd.DataFrame({'a':xt,'inat_taxon_id':np.arange(240)+2000,'genus':'NeverSeen'})
    truth=np.column_stack([.1+.8*xt*xt,*[(.9-.8*xt*xt)/3]*3])
    nonlinear=m.predict_nonlinear(tr,te,('a',))
    ridge=m.predict(tr,te,('a',),'absolute')
    assert np.square(nonlinear-truth).sum(1).mean() < .1*np.square(ridge-truth).sum(1).mean()


def test_complete_synthetic_replay_saves_identical_targets_and_all_four_predictions(m,tmp_path):
    rng=np.random.default_rng(170)
    raw=pd.DataFrame({'inat_taxon_id':np.arange(600),
        'genus':['G'+str(i).zfill(4) for i in range(600)],
        'photo_id':np.arange(600),'latitude':rng.uniform(-20,20,600),
        'longitude':np.tile([-135.,-45.,45.,135.],150),
        'morph':[m.CLASSES[i%4] for i in range(600)],
        'flower_effective_pixels':30000,'nuisance_pixel_fraction':.1})
    for f in [*m.GEO,*m.FULL]:raw[f]=rng.normal(size=len(raw))
    # A subset occupies two sectors; complete prediction coverage still weights
    # each species once after averaging its two saved sector errors.
    extra=raw.iloc[:20].copy()
    extra['photo_id']+=1000
    extra['longitude']=-extra.longitude
    raw=pd.concat([raw,extra],ignore_index=True)
    exclusions,_=m.spatial_exclusions(raw)
    result=m.evaluate_case(raw,False,exclusions,99,tmp_path)
    d=pd.read_csv(tmp_path/'soil_all_oof_predictions.csv.gz')
    audit=pd.read_csv(tmp_path/'soil_all_fold_audit.csv')
    assert len(audit)==60 and audit.status.eq('OK').all()
    assert result['n_input_species']==600 and result['n_input_photos']==620
    keys=[]
    for scheme in m.SCHEMES:
        subset=d[d.scheme.eq(scheme)]
        m.validate_records(subset)
        assert subset.inat_taxon_id.nunique()==600
        assert subset.n_photos.sum()==620
        keys.append(subset[['inat_taxon_id','sector','genus','fold']].sort_values(['inat_taxon_id','sector']).reset_index(drop=True))
        assert result['schemes'][scheme]['overall']['n_species']==600
    assert all(k.equals(keys[0]) for k in keys)
