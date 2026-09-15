from pathlib import Path
import json

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
