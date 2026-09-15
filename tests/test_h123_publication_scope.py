from pathlib import Path
import json
import re
import csv
import hashlib

import yaml

ROOT = Path(__file__).resolve().parents[1]


def test_legacy_reproduction_is_manual_only():
    workflow = yaml.load((ROOT / '.github/workflows/34species-paper.yml').read_text(), Loader=yaml.BaseLoader)
    assert set(workflow['on']) == {'workflow_dispatch'}
    assert 'reproduce-paper' in workflow['jobs']


def test_current_publication_scope_preserves_closed_pixels():
    status = (ROOT / 'docs/FCP_H123_PUBLICATION_STATUS_20260915.md').read_text(encoding='utf-8')
    assert 'not submission-ready' in status
    assert 'P500 opening remains unauthorized' in status
    assert '34-species literature comparison is retired' in status


def test_manuscript_h1_values_match_direct_receipt():
    result = json.loads((ROOT / 'results/polymorphism_h1_observer_disjoint_d_20260913/result.json').read_text())
    reserve = result['analyses']['reserve']['D']
    manuscript = (ROOT / 'docs/FCP_H123_MANUSCRIPT.md').read_text(encoding='utf-8')
    for value in [reserve['spearman_rho'], reserve['lin_ccc'], *reserve['bootstrap_95_percentile_ci']]:
        assert f'{value:.10f}' in manuscript
    assert 'failed its 0.80 criterion' in manuscript
    assert 'not an untouched prospective confirmation' in manuscript
    assert 'not submission-ready' in manuscript


def test_manuscript_h2_values_match_targeted_receipt():
    result = json.loads((ROOT / 'results/polymorphism_white_axis_targeted_test_20260912/result.json').read_text())
    manuscript = (ROOT / 'docs/FCP_H123_MANUSCRIPT.md').read_text(encoding='utf-8')
    section = manuscript.split('### H2:')[1].split('### H3:')[0]
    for threshold in result['thresholds'].values():
        for cohort in ('discovery', 'reserve'):
            values = threshold[cohort]
            assert f"W={values['observed_mean_squared_white_axis_alignment']:.6f}" in section
            assert f"{values['species']} species" in section
            assert f"p={values['structured_null_upper_p']:.3f}" in section
    assert 'not prospective confirmation' in section


def test_active_manuscript_links_and_measurement_limits():
    for name in ('FCP_H123_MANUSCRIPT.md', 'FCP_H123_SUPPLEMENT.md', 'FCP_H123_PUBLICATION_STATUS_20260915.md'):
        path = ROOT / 'docs' / name
        content = path.read_text(encoding='utf-8')
        for target in re.findall(r'\]\(([^)]+)\)', content):
            if not target.startswith('https://'):
                assert (path.parent / target).is_file(), target
    manuscript = (ROOT / 'docs/FCP_H123_MANUSCRIPT.md').read_text(encoding='utf-8')
    for text in ('Monarda region-agreement gate failed', 'not true biological range size',
                 'No pigment pathway', 'P500 is currently\nmetadata-only'):
        assert text in manuscript


def test_h3_artifact_summaries_and_manuscript_values():
    base = ROOT / 'docs/supporting/h123'
    expected = {
        'h3a_signal_by_scenario.csv': '8b37b0edcdf0823351aad82d443d2ff6686d3684294855f20900e1770614f986',
        'h3b_span_summary.csv': 'a8f58d2cf7282ec8aa55fd5741eb34e9440112fd319710654c4f03c12551460b',
    }
    manuscript = (ROOT / 'docs/FCP_H123_MANUSCRIPT.md').read_text(encoding='utf-8')
    section = manuscript.split('### H3:')[1].split('## Discussion')[0]
    for name, digest in expected.items():
        raw = (base / name).read_bytes().replace(b'\r\n', b'\n')
        assert hashlib.sha256(raw).hexdigest() == digest
        rows = list(csv.DictReader(raw.decode().splitlines()))
        if name.startswith('h3a'):
            reserve = [r for r in rows if r['cohort'] == 'reserve']
            assert {r['scenario'] for r in reserve} == {'S1', 'S2', 'S3'}
            for row in reserve:
                assert f"{float(row['p_K_raw']):.3f}" in section
        else:
            assert {r['cohort'] for r in rows} == {'discovery', 'reserve'}
            for row in rows:
                assert f"rho={float(row['rho_D_span']):.6f}" in section
                assert f"p={float(row['p_D_span']):.6f}" in section
