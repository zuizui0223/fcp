from pathlib import Path

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
