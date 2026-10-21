"""Synthetic and reproducibility guards: artificial records are not biological evidence."""
from __future__ import annotations

import importlib.util
import json
import sqlite3
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from fcp_pipeline.fcp_allphotos import (
    opportunity_frame, opportunity_schedule, init_index, ingest_photos,
    species_ledger, coverage_receipt, main
)


def synthetic_frame(tmp_path):
    p=tmp_path/'frame.csv'
    pd.DataFrame({"inat_taxon_id":[15,12,17],"species":["Lathyrus spicatus","Lathyrus rubra","Lathyrus alba"],"after_observer_cap":[200,30,0]}).to_csv(p,index=False)
    return p


def measured_table(tmp_path):
    p=tmp_path/'measured.csv'
    rows=[]
    # 70 valid photos from two observers, one observer dominating with 50 white
    # pictures. Repeated photos from same observation remain separate inputs.
    for i in range(70):
        rows.append({"inat_taxon_id":15,"species":"Lathyrus spicatus",
                     "photo_id":f"A{i}","observation_id":f"A{i//2}","observer_id":"big" if i<50 else "small",
                     "morph":"white" if i<50 else "red_pink","global_classifiable":True})
    for i in range(26):
        rows.append({"inat_taxon_id":12,"species":"Lathyrus rubra",
                     "photo_id":f"B{i}","observation_id":f"B{i}","observer_id":f"r{i%8}",
                     "morph":"blue_purple","global_classifiable":True})
    pd.DataFrame(rows).to_csv(p,index=False)
    return p


def test_no_photo_depth_or_species_selection_cap(tmp_path):
    source=synthetic_frame(tmp_path)
    frame=opportunity_frame(source,3)
    schedule=opportunity_schedule(frame,3)
    assert set(schedule.inat_taxon_id)=={15,12,17}
    assert schedule.acquisition_cap.eq("none").all()
    assert schedule.set_index('inat_taxon_id').loc[17,'priority_tier']=='currently_lt20'
    assert schedule.set_index('inat_taxon_id').loc[12,'priority_tier']=='available_20_39'
    assert len(set(schedule.partition)) <= 3
    with pytest.raises(ValueError,match='changed from 4'):
        opportunity_frame(source,4)


def test_all_photos_and_three_distinct_statuses(tmp_path):
    frame=opportunity_frame(synthetic_frame(tmp_path),3)
    path=measured_table(tmp_path)
    with sqlite3.connect(tmp_path/'staging.sqlite') as db:
        init_index(db)
        h,n=ingest_photos(db,{"fcp_new":path},frame,chunksize=7)
        assert n==96
        led=species_ledger(db,frame).set_index('inat_taxon_id')
    assert led.loc[15,'n_photos']==70  # not first 40 or first 100
    assert led.loc[15,'n_classifiable']==70
    assert led.loc[15,'n_observations']==35
    assert led.loc[15,'n_classifiable_observers']==2
    assert led.loc[15,'measured_status']=='multiple_photo_states_detected'
    assert led.loc[12,'measured_status']=='insufficient_classifiable'
    assert np.isnan(led.loc[12,'D_raw'])
    assert led.loc[17,'measured_status']=='not_measured'
    assert np.isnan(led.loc[17,'D_raw'])
    assert led.loc[15,'D_raw'] > 0
    assert led.loc[15,'D_observer_equal'] == pytest.approx(0.5)
    assert led.loc[15,'D_observer_equal'] > led.loc[15,'D_raw']
    assert led.loc[15,'H2_coarse_0_10_pass']
    rep=coverage_receipt(frame,led.reset_index(),h,n)
    assert rep['n_D_evaluable_species']==1
    assert rep['n_H2_coarse_gate_species']==1
    assert rep['frozen_original_H2_recomputed'] is False


def test_more_photos_same_observation_not_pseudoreplicates(tmp_path):
    frame=opportunity_frame(synthetic_frame(tmp_path),3)
    rows=[]
    # 100 photographs do not guarantee 100 distinct observations / plants.
    for j in range(100):
        rows.append({"inat_taxon_id":15,"species":"Lathyrus spicatus",
                     "photo_id":str(j+1),"observation_id":"one-observation",
                     "observer_id":"one-observer","global_classifiable":True,
                     "morph":"white" if j<50 else "red_pink"})
    path=tmp_path/'sameobs.csv';pd.DataFrame(rows).to_csv(path,index=False)
    with sqlite3.connect(':memory:') as db:
        init_index(db)
        ingest_photos(db,{"new":path},frame)
        d=species_ledger(db,frame).set_index('inat_taxon_id')
    assert d.loc[15,'n_photos']==100
    assert d.loc[15,'n_classifiable_observations']==1
    assert d.loc[15,'D_raw']==pytest.approx(0.5)
    # This is photo variation, not replicated, genetic local FCP.
    assert d.loc[15,'n_classifiable_observers']==1


def test_duplicate_photo_id_fails_across_inputs(tmp_path):
    frame=opportunity_frame(synthetic_frame(tmp_path),3)
    path=measured_table(tmp_path)
    with sqlite3.connect(':memory:') as db:
        init_index(db)
        with pytest.raises(ValueError, match='Duplicate photo_id'):
            ingest_photos(db,{"a":path,"b":path},frame,chunksize=30)


def test_readiness_nclass40_not_polymorphism_prevalence(tmp_path):
    frame=opportunity_frame(synthetic_frame(tmp_path),3)
    path=tmp_path/'photos.csv'
    pd.DataFrame([{'inat_taxon_id':15,'species':'Lathyrus spicatus','photo_id':f'P{i}',
                   'observation_id':f'Q{i}','observer_id':f'O{i}',
                   'morph':'white','global_classifiable':True} for i in range(100)]).to_csv(path,index=False)
    with sqlite3.connect(':memory:') as db:
        init_index(db);ingest_photos(db,{'fresh':path},frame)
        d=species_ledger(db,frame).set_index('inat_taxon_id')
    assert d.loc[15,'D_raw']==0.0
    assert not bool(d.loc[15,'H2_coarse_0_10_pass'])
    assert d.loc[15,'measured_status']=='single_photo_state_detected'


def test_main_outputs_manifest_not_artificial_photo_cap(tmp_path):
    frame=synthetic_frame(tmp_path);photos=measured_table(tmp_path);out=tmp_path/'output'
    assert main(['--opportunity',str(frame),'--expected-species','3',
                 '--measured',f'new={photos}','--outdir',str(out)])==0
    a=pd.read_csv(out/'all_species_acquisition_schedule.csv')
    b=json.loads((out/'all_photo_FCP_receipt.json').read_text())
    assert len(a)==3 and a.acquisition_cap.eq('none').all()
    assert b['n_unique_measured_photo_ids']==96
    assert b['n_D_evaluable_species']==1
    assert b['source_frame_not_global_flora_or_FCP_prevalence_denominator']
