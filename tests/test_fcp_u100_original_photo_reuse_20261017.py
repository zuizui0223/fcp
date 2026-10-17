"""No-pixel synthetic protocol tests for preserving original FCP 100-photo/taxon rule."""
import pandas as pd
import pytest

from scripts.analysis import fcp_u100_original_photo_reuse_20261017 as m


def fake_group(tid, *, photo_start, obs_start):
    return pd.DataFrame({
        "photo_id":list(range(photo_start,photo_start+100)),
        "observation_id":list(range(obs_start,obs_start+100)),
        "inat_taxon_id":[tid]*100,
        "species":[f"Testus species{tid}"]*100,
        "global_classifiable":[True]*100,
        "morph":["white"]*100,
    })


def fake_source(monkeypatch):
    new=pd.concat(
        [fake_group(t,photo_start=100001+100*i,obs_start=700001+100*i)
         for i,t in enumerate([101,102,103,104])],ignore_index=True)
    old={
      "discovery":pd.DataFrame({
          "photo_id":[100005,100205,999981],
          "observation_id":[199999,199998,700102],
      }),
      "validation":pd.DataFrame({
          "photo_id":[88901],"observation_id":[98801],
      }),
      "prospective":pd.DataFrame({
          "photo_id":[88902],"observation_id":[98802],
      }),
      "wave4":pd.DataFrame({
          "photo_id":[88903],"observation_id":[98803],
      }),
    }
    monkeypatch.setattr(m,"EXPECTED_NEW_H9_SOURCE_PHOTOS",400)
    monkeypatch.setattr(m,"EXPECTED_NEW_H9_COMPLETE_TAXA",4)
    monkeypatch.setattr(m,"PINNED_OLD_DISCOVERY_PHOTO_COLLISIONS",frozenset({100005,100205}))
    monkeypatch.setattr(m,"PINNED_OLD_DISCOVERY_OBSERVATION_COLLISIONS",frozenset({700102}))
    return new,old


def test_entire_taxon_excluded_if_any_of_its_original_photos_or_observations_reused(monkeypatch):
    new,old=fake_source(monkeypatch)
    result,audit=m.inspect_reused_new_u100_source(new,old)
    assert audit["affected_new_H9_taxon_IDs"]==[101,102,103]
    assert audit["affected_new_H9_complete_species"]==3
    assert audit["original_reused_photo_ID_count"]==2
    assert audit["original_reused_observation_ID_count"]==1
    assert audit["original_source_photo_rows_removed_whole_taxon"]==300
    assert audit["effective_new_U100_complete100_taxa"]==1
    assert len(result)==100
    assert result.inat_taxon_id.unique().tolist()==[104]
    assert audit["original_50km_H1_D_W_estimators_unchanged"] is True


def test_no_99_photo_species_imputed_by_only_deleting_colliding_photo(monkeypatch):
    new,old=fake_source(monkeypatch)
    filtered,audit=m.inspect_reused_new_u100_source(new,old)
    assert filtered.groupby("inat_taxon_id").size().eq(100).all()
    assert len(filtered)==len(new)-100*audit["affected_new_H9_complete_species"]
    assert not set(filtered.photo_id).intersection({100005,100205})


def test_fail_closed_if_source_observation_overlap_not_discovery_only(monkeypatch):
    new,old=fake_source(monkeypatch)
    old["validation"].loc[0,"observation_id"]=700350
    with pytest.raises(RuntimeError,match="fingerprints drifted"):
        m.inspect_reused_new_u100_source(new,old)


def test_fail_closed_on_unexpected_extra_photo_collision(monkeypatch):
    new,old=fake_source(monkeypatch)
    old["wave4"].loc[0,"photo_id"]=100350
    with pytest.raises(RuntimeError,match="fingerprints drifted"):
        m.inspect_reused_new_u100_source(new,old)


def test_fail_closed_if_one_original_photo_missing_without_filling_from_api(monkeypatch):
    new,old=fake_source(monkeypatch)
    new=new.iloc[:-1].copy()
    with pytest.raises(RuntimeError,match="exactly 100 photos per taxon"):
        m.inspect_reused_new_u100_source(new,old)


def test_fail_closed_if_source_photo_id_duplicated_within_new_h9(monkeypatch):
    new,old=fake_source(monkeypatch)
    new.loc[1,"photo_id"]=new.loc[0,"photo_id"]
    with pytest.raises(RuntimeError,match="duplicate photographic"):
        m.inspect_reused_new_u100_source(new,old)


def test_source_with_no_collision_retains_all_only_without_pinned_real_anomaly(monkeypatch):
    new,old=fake_source(monkeypatch)
    old["discovery"]["photo_id"]=[1,2,3]
    old["discovery"]["observation_id"]=[4,5,6]
    retained,report=m.inspect_reused_new_u100_source(new,old,strict_original=False)
    assert len(retained)==400
    assert report["affected_new_H9_complete_species"]==0
    with pytest.raises(RuntimeError,match="fingerprints drifted"):
        m.inspect_reused_new_u100_source(new,old,strict_original=True)
