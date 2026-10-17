"""Synthetic tests for 4730 complete U100 source coverage and fail-closed gates."""
import pandas as pd
import pytest

from scripts.acquisition import fcp_full_u100_4730_20261017 as f


def setup():
    # 7 true original U100 parent taxa, 3 historical high-depth outside P100.
    p100=[{"inat_taxon_id":i,"species":f"Species_{i}","after_observer_cap":101} for i in range(7)]
    old={7:"Species_7",8:"Species_8",9:"Species_9"}
    wave={5:"Species_5",6:"Species_6"}
    return p100,old,wave


def test_all_original_source_species_included_and_only_unmeasured_selected():
    p100,old,wave=setup()
    allx,pending,ledger=f.construct_all_species_inventory(
        p100,old,wave,u100=10,old_n=3,wave_n=2)
    assert len(allx)==10 and allx.inat_taxon_id.nunique()==10
    assert pending.inat_taxon_id.tolist()==[0,1,2,3,4]
    assert ledger["new_unmeasured_species"]==5
    assert ledger["new_max_metadata_photos"]==500
    assert ledger["full_U100_if_every_species_authorized_photo_rows"]==1000
    assert ledger["species_selection_conditional_on_flower_colour"] is False
    assert set(allx.u100_status)=={
        "original_149900_photo_measured",
        "wave4_49800_photo_id_authorized",
        "unmeasured_original_u100_acquire",
    }


def test_drop_or_duplicate_original_identity_fails_closed():
    p100,old,wave=setup()
    with pytest.raises(RuntimeError,match="wave4 not inside original"):
        f.construct_all_species_inventory(
            p100[:-1],old,wave,u100=10,old_n=3,wave_n=2)


def test_wave4_original_overlap_rejected():
    p100,old,wave=setup()
    wave=dict(wave)
    wave[7]="Species_7"
    with pytest.raises(RuntimeError,match="wave4 not inside original"):
        f.construct_all_species_inventory(p100,old,wave,u100=10,old_n=3,wave_n=3)


def test_taxon_name_mismatch_halts_without_relabel():
    p100,old,wave=setup()
    wave=dict(wave)
    wave[5]="Different label"
    with pytest.raises(RuntimeError,match="wave4 taxon identity/name"):
        f.construct_all_species_inventory(p100,old,wave,u100=10,old_n=3,wave_n=2)


def test_pending_has_all_failed_wave4_and_third_as_target_not_excluded():
    p100,old,wave=setup()
    wave.pop(6)
    _,pending,ledger=f.construct_all_species_inventory(
        p100,old,wave,u100=10,old_n=3,wave_n=1)
    assert 6 in set(pending.inat_taxon_id)
    assert ledger["new_unmeasured_species"]==6
    assert ledger["all_species_included_even_if_prior_H9_attempt_failed"] is True


def test_original_4730_denominator_fixed():
    assert f.U100==4730
    assert f.REMAINING==2733
    assert f.PER_SPECIES==100
    assert f.SHARDS==4
    assert f.ORIGINAL_SPECIES+f.WAVE4_SPECIES+f.REMAINING==f.U100


def test_no_morph_or_environment_consulted_in_inventory():
    p100,old,wave=setup()
    assert set(f.construct_all_species_inventory(
        p100,old,wave,u100=10,old_n=3,wave_n=2)[1].columns)=={
            "inat_taxon_id","species"
        }
