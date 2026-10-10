"""All abiotic blocks preserve real Brownian and spatial kernel baseline."""
from __future__ import annotations
import sys
from pathlib import Path
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
import audit_fcp_42111_all_environment_beyond_true_phylo_spatial_20261010 as M


def test_all_seven_environment_groups_and_36_fields_retained():
    z=M.BLOCKS
    assert list(z)==["elevation","temperature","precipitation","solar_radiation","wind","vapor_pressure","soil"]
    assert sum(len(v) for v in z.values())==36
    assert z["soil"]==("soil_pH","soil_SOC","soil_N","soil_clay",
                         "soil_available_water_proxy","soil_cec_0_30cm_source_raw",
                         "soil_sand_0_30cm_source_raw","soil_silt_0_30cm_source_raw",
                         "soil_bdod_0_30cm_source_raw","soil_cfvo_0_30cm_source_raw")


def test_wrapper_supplies_exact_blocks_and_preserves_original_derived_result(monkeypatch):
    called=[]
    def fake(source,ledger,trees,*,strict,nboot,blocks):
        called.append({"strict":strict,"nboot":nboot,"blocks":blocks})
        return {"schema":"prior_source","cohorts":{"250km_climate":{"models":{}}},
                "original_42111_photo_species_denominator":42111}
    monkeypatch.setattr(M.P,"run",fake)
    out=M.explore(None,None,{},strict=False,nboot=9)
    assert out["schema"]==M.SCHEMA
    assert out["all_36_predictors_scanned_without_outcome_selected_drop"]
    assert out["full_42111_phylogenetically_controlled_result_NOT_identified"]
    assert called[0]["blocks"]==M.BLOCKS
    assert called[0]["nboot"]==9
