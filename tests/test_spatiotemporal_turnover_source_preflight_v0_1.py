from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]
P=ROOT/'development'/'spatiotemporal_turnover_coupling_v0_1'/'source_preflight_result_v0_1.json'

def test_source_preflight_is_outcome_blind_and_frozen():
    x=json.loads(P.read_text())
    assert x['status']=='SPATIOTEMPORAL_SOURCE_ELIGIBILITY_PASS_BOTH_ARMS'
    assert x['biological_turnover_outcomes_opened'] is False
    assert x['trait_arm']['rows']==91970
    assert x['trait_arm']['traits_seen']==18
    assert x['trait_arm']['traits_passing_paired_source_gate']==14
    assert [t['spatial_eligible_species'] for t in x['trait_arm']['primary_traits']]==[88,64,79,34]
    fam=x['interaction_arm']['fixed_primary_families']
    assert [(z['interaction_type'],z['focal_role'],z['reference_taxa_in_ge3_networks']) for z in fam]==[
      ('mutualism','either',305),('predation','from',132),('parasitism','from',40)
    ]
    assert x['interaction_arm']['all_primary_families_pass'] is True
    assert x['active_fcp_manuscript_changed'] is False
    assert x['chun_el_v0_3_changed'] is False
