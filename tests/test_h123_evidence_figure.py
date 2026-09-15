import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('h123_figure', ROOT / 'scripts/analysis/make_h123_evidence_figure.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_figure_data_preserves_comparisons_and_interval_types():
    data = module.figure_data()
    assert len(data['h1']) == 3 and len(data['h2']) == 4
    assert len(data['h3a']) == 6 and len(data['h3b']) == 2
    assert data['h1'][0]['estimate'] > .8 > data['h1'][1]['estimate']
    assert data['h1'][2]['interval_type'] == 'partition 5th–95th percentiles'
    assert data['h1'][2]['n'] == 329
    assert {r['scenario'] for r in data['h3a']} == {'S1', 'S2', 'S3'}
    assert all(float(r['p_K_raw']) > .05 for r in data['h3a'] if r['cohort'] == 'reserve')
    assert float(data['h3b'][1]['rho_D_span']) < 0
    assert data['p500_opened'] is data['scientific_results_recomputed'] is False


def test_committed_plot_data_matches_current_frozen_sources():
    saved = json.loads((ROOT / 'docs/figures/h123/h123_evidence_data.json').read_text())
    assert saved == module.figure_data()


def test_render_is_repeatable_without_source_writes(tmp_path):
    before = module.figure_data()['source_sha256_lf']
    module.render(tmp_path / 'one')
    module.render(tmp_path / 'two')
    for name in ('h123_evidence.png', 'h123_evidence.pdf', 'h123_evidence_data.json'):
        assert (tmp_path / 'one' / name).read_bytes() == (tmp_path / 'two' / name).read_bytes()
    assert (tmp_path / 'one/h123_evidence.png').read_bytes().startswith(b'\x89PNG')
    assert (tmp_path / 'one/h123_evidence.pdf').read_bytes().startswith(b'%PDF')
    assert before == module.figure_data()['source_sha256_lf']
