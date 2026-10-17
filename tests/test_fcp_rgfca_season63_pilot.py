"""Only synthetic source metadata. No original photo downloaded."""
import pytest
from scripts.acquisition.run_fcp_rgfca_season63_pilot import photo_ids,verify_row,cell,month_qc
def test_precolour_63_fixed_source_ids():
    assert len(photo_ids())==63
def test_source_date_and_accuracy_cannot_be_forged():
    row=dict(inat_taxon_id="740971",species="Ranunculus formosa-montanus",
        latitude="24.01",longitude="121.27",observed_on="2026-06-04",
        positional_accuracy_m="200",source_snapshot_date="2026-09-27",
        observation_id="o1",observer_id="u1")
    assert cell(row)==(456,1205)
    verify_row(row,"740971","2026-06")
    for k,v in [("longitude","-100"),("observed_on","2026-07-04"),
                ("positional_accuracy_m",""),("positional_accuracy_m","7000")]:
        with pytest.raises(ValueError):verify_row({**row,k:v},"740971","2026-06")
def test_same_roi_two_organs_never_establish_colour_polymorphism():
    b=dict(inat_taxon_id="740971",year_month="2026-06",observer_id="u1",
           observation_id="o1",morph="white",measurement_status="classified_four_state_morph",
           flower_fraction_white="1",flower_fraction_yellow="0")
    records=[b,{**b,"observation_id":"o2","observer_id":"u2","morph":"yellow_orange",
                 "flower_fraction_white":".3","flower_fraction_yellow":".3"},
             {**b,"year_month":"2026-07","observation_id":"o3","observer_id":"u3"},
             {**b,"year_month":"2026-07","observation_id":"o4","observer_id":"u4"}]
    months,contrasts=month_qc(records)
    assert months[0]["n_usable_independent_observations"]==1
    assert months[0]["n_same_roi_white_yellow_mixed_excluded"]==1
    assert contrasts[0]["both_months_two_obs_and_observers"] is False
