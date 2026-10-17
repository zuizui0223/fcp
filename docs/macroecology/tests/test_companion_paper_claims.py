"""Audit claims and documented source identities in the FCP macroecology manuscript.

These guards check this *draft's claim contract*, not substantive new statistical
results. Parent analyses have independent frozen source-replay workflows.
"""
from pathlib import Path
import json
import re

ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / 'MANUSCRIPT_V0_1.md'
LEDGER = ROOT / 'EVIDENCE_AND_CLAIM_LEDGER.md'
STRATEGY = ROOT / 'PAPER_STRATEGY_AND_GATES.md'
TABLES = ROOT / 'MAIN_TABLES.md'
SRC = ROOT / 'SOURCE_LINKS.json'


def test_all_documents_exist_and_scope_is_companion_not_nph_main():
    for p in [MANUSCRIPT, LEDGER, STRATEGY, TABLES, SRC]:
        assert p.exists() and p.stat().st_size > 1_000
    s = MANUSCRIPT.read_text()
    assert 'retrospective' in s.lower()
    assert 'independent main FCP manuscript' in s
    assert 'New Phytologist manuscript' in s
    assert 'not submitted' in s


def test_frozen_source_denominators_and_independent_measurement_limit():
    s = MANUSCRIPT.read_text()
    for q in ['100,543', '42,111', '45,953', '54,590', '45,836', '25,028',
              '31,925', '20,786', '2,476', '670', '58 genera', '9,370', '2,886']:
        assert q in s,q
    assert 'not independent species-occurrence' in s or 'independence of **phenotype measurement**' in s
    assert 'environmental sampling locations still came from' in s


def test_directional_primary_and_negative_external_are_jointly_reported():
    s = MANUSCRIPT.read_text()
    for q in ['+2.255', '+0.319', '+4.190', '+3.952', '+1.199', '+6.704',
              '+1.229', '−2.046', '+0.808', '−3.556', '+1.326', '+0.518',
              '−1.253', '+3.178', '−0.258', '+3.56', '−15.34', '+22.60']:
        assert q in s,q
    assert 'could not' not in s.split('## Conclusions')[-1].lower()  # avoid careless causal negatives
    assert 'a universal climatic colour response and its evolutionary mechanism are not' in s


def test_prediction_gain_is_not_pigment_or_fitness_effect():
    s = MANUSCRIPT.read_text()
    assert '0.002943' in s and '0.001432' in s
    assert 'prediction error' in s and 'relative survival or reproduction' in s
    assert 'one phenotype measurement' not in s.lower()
    assert 'genetic maintenance' in s.lower()
    assert '54,590' in s


def test_journal_submission_blockers_are_explicit_and_no_invented_metadata():
    s = MANUSCRIPT.read_text()
    assert '[to be confirmed]' in s
    for q in ['overlapping photographs', 'rights', 'DOI archive', 'not yet finalized']:
        assert q in s,q
    t = STRATEGY.read_text()
    assert '投稿' in t and 'M2' in t and 'M3' in t and 'M4' in t


def test_evidence_source_contract_points_to_exact_analyses():
    d = json.loads(SRC.read_text())
    for q in ['macro_direction','same_species_matched_spatial','joint_36_prediction',
              'genus_longitude_transport','independent_botanical_colour_measurement']:
        assert d[q].startswith('https://github.com/zuizui0223/fcp/pull/')
    assert d['source_photo_sha256'] == '2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0'
    assert d['source_abiotic_sha256'] == '6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca'


def test_no_unsupported_natural_selection_published_claim():
    s = MANUSCRIPT.read_text().lower()
    forbidden = ['demonstrates adaptive white evolution', 'universal selection for white',
                 'we identify the genetic maintenance mechanism',
                 'the first global study of flower colour and climate']
    assert not any(q in s for q in forbidden)


def test_main_figures_do_not_declare_a_confirmed_world_hotspot_map():
    s = (ROOT / 'FIGURE_TABLE_PLAN.md').read_text().lower()
    assert 'hotspot' in s
    assert 'do not' in s


def test_tables_use_disjoint_source_measures_and_explicit_ci_families():
    s = TABLES.read_text()
    assert 'different species-level binary colour endpoint' in s.lower()
    assert 'world' in s.lower() and 'paired' in s.lower()
    assert 'simultaneous maxT' in s
