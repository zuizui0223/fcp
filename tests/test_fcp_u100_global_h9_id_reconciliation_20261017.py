"""Original FCP H9 four-shard reconciliation: synthetic, never biological data."""
import pandas as pd
import pytest
from scripts.acquisition import fcp_full_u100_4730_20261017 as m


def candidate(tid,ids,*,observer_ids=None):
    return pd.DataFrame({
        "inat_taxon_id":[tid]*len(ids),
        "species":[f"synthetic_{tid}"]*len(ids),
        "observation_id":observer_ids or [tid*100+k for k in range(len(ids))],
        "photo_id":ids,
        "h9_selection_order":list(range(len(ids))),
    })


def run(frames,taxa):
    return m.reconcile_historical_global_seen(frames,set(taxa),target_n=3)


def test_historical_taxon_order_reserves_photo_across_all_four_shards():
    frames=[
       candidate(200,[3,4,5]),
       candidate(100,[1,2,3]),
       candidate(400,[7,8,9]),
       candidate(300,[9,10,11]),
    ]
    photo,audit,ledger=run(frames,[100,200,300,400])
    # Original H9 processes 100,200,300,400 not shard presentation order.
    # Later taxa 200 and 400 share one earlier photo and cannot fabricate
    # a fourth replacement from iNaturalist's discarded raw API page.
    assert set(photo.inat_taxon_id)=={100,300}
    assert set(photo.photo_id)=={1,2,3,9,10,11}
    assert len(photo)==6
    assert ledger["globally_seen_cross_shard_rejected_source_candidate_rows"]==2
    assert ledger["post_source_global_100photo_species"]==2
    assert audit.set_index("inat_taxon_id").loc[200,"global_retained_original_H9_rows"]==2
    assert ledger["original_four_shard_inputs_untouched"] is True


def test_incomplete_earlier_taxon_also_reserves_original_selected_photo_IDs():
    frames=[
      candidate(200,[9,10,11]),
      candidate(100,[9]), # full source short, but original freeze_H9 DOES reserve ID 9
      candidate(300,[12,13,14]),
      candidate(400,[15,16,17]),
    ]
    keep,audit,receipt=run(frames,[100,200,300,400])
    assert set(keep.inat_taxon_id)=={300,400}
    assert receipt["globally_seen_cross_shard_rejected_source_candidate_rows"]==1
    assert audit.set_index("inat_taxon_id").loc[200,"global_retained_original_H9_rows"]==2


def test_shared_observation_id_even_with_distinct_image_ids_causes_conservative_exclusion():
    frames=[
      candidate(100,[1,2,3],observer_ids=[91,92,93]),
      candidate(200,[10,11,12],observer_ids=[93,94,95]),
      candidate(300,[20,21,22]),
      candidate(400,[30,31,32]),
    ]
    keep,audit,_=run(frames,[100,200,300,400])
    assert len(keep)==9
    assert 200 not in set(keep.inat_taxon_id)
    assert keep.observation_id.nunique()==len(keep)


def test_original_source_does_not_fill_missing_100_using_untested_supplemental_photos():
    frames=[
       candidate(200,[3,4,5]),
       candidate(100,[1,2,3]),
       candidate(300,[6,7,8]),
       candidate(400,[9,10,11]),
    ]
    before={int(tid):set(df.photo_id) for tid,df in zip([200,100,300,400],frames)}
    observed,_audit,_ledger=run(frames,[100,200,300,400])
    for tid,g in observed.groupby("inat_taxon_id"):
        assert set(g.photo_id).issubset(before[tid])
        assert len(g)==3


def test_reject_duplicate_source_within_species_not_silently_cleaned():
    frames=[
      candidate(100,[1,2,2]),
      candidate(200,[5,6,7]),
      candidate(300,[9,10,11]),
      candidate(400,[13,14,15]),
    ]
    with pytest.raises(RuntimeError,match="within-species"):
        run(frames,[100,200,300,400])


def test_incomplete_four_source_shards_never_gain_pixel_permission():
    with pytest.raises(RuntimeError,match="four original H9 candidate"):
        m.reconcile_historical_global_seen([candidate(100,[1,2,3])],{100},target_n=3)
