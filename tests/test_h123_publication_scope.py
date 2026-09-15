from pathlib import Path
import json
import re
import csv
import hashlib
from statistics import median
from decimal import Decimal

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


def test_monarda_table_recomputed_from_unchanged_pixel_counts():
    base = ROOT / 'data/validation/monarda_region_agreement_v1'
    result_raw = (base / 'monarda_region_agreement_result_v1.json').read_bytes()
    rows_raw = (base / 'monarda_region_agreement_rows_v1.csv').read_bytes()
    assert hashlib.sha256(result_raw).hexdigest() == '49b2f017a6accb305867d62e7b4afe4e8ed96fcd03f8aafde86bda2da915fa4c'
    assert hashlib.sha256(rows_raw).hexdigest() == '67cbeaa2cd0ea1f1d3f98516b4fe9f7bb1abc8648ac80268e6dc4ec3e754dc22'
    result = json.loads(result_raw)
    rows = list(csv.DictReader(rows_raw.decode().splitlines()))
    assert len(rows) == 110
    assert len({(r['split'], r['coco_image_id']) for r in rows}) == 110
    assert all(r['model_status'] == 'model_success' for r in rows)
    positive = [r for r in rows if r['reference_role'] == 'positive_generic_flower_region']
    assert len(positive) == 109
    assert result['reference_unknown_images'] == 1
    assert result['unknown_zero_annotation_image']['counted_as_verified_negative'] is False

    def count(row, field):
        value = Decimal(row[field])
        assert value.is_finite() and value >= 0 and value == int(value)
        return int(value)

    counts = [(count(r, 'reference_pixels'), count(r, 'predicted_pixels'),
               count(r, 'intersection_pixels')) for r in positive]
    assert all(0 <= intersection <= min(reference, prediction)
               for reference, prediction, intersection in counts)
    assert sum(prediction == 0 for _, prediction, _ in counts) == 15
    reference, prediction, intersection = map(sum, zip(*counts))
    assert (reference, prediction, intersection) == (28428697, 21316820, 12113123)
    values = [intersection / prediction, intersection / reference,
              median(i / p if p else 0.0 for _, p, i in counts)]
    keys = ['pooled_prediction_precision', 'pooled_reference_recall',
            'median_positive_image_prediction_precision']
    assert values == [result['aggregate'][key] for key in keys]
    assert [v >= floor for v, floor in zip(values, [.70, .35, .70])] == [False, True, False]
    assert result['limited_gate_pass'] is False
    supplement = (ROOT / 'docs/FCP_H123_SUPPLEMENT.md').read_text(encoding='utf-8')
    labels = ['Pooled prediction precision', 'Pooled reference recall',
              'Median annotated-image prediction precision']
    for label, value, floor, decision in zip(labels, values, [.70, .35, .70], ['Failed', 'Passed', 'Failed']):
        assert f'| {label} | {value:.8f} | {floor:.2f} | {decision} |' in supplement
    manuscript = (ROOT / 'docs/FCP_H123_MANUSCRIPT.md').read_text(encoding='utf-8')
    for value in values:
        assert f'{value:.8f}' in manuscript
    assert 'not a biological absence' in manuscript


def test_h2_labeled_table_matches_frozen_targeted_result():
    result = json.loads((ROOT / 'results/polymorphism_white_axis_targeted_test_20260912/result.json').read_text())
    supplement = (ROOT / 'docs/FCP_H123_SUPPLEMENT.md').read_text(encoding='utf-8')
    for key, label in [('primary_0_10', 'Primary 0.10'), ('strict_0_20', 'Strict 0.20')]:
        for cohort in ('discovery', 'reserve'):
            r = result['thresholds'][key][cohort]
            observed = r['observed_mean_squared_white_axis_alignment']
            median_null = r['structured_null_summary']['q50']
            excess = observed - median_null
            assert abs(excess - r['observed_minus_null_median']) < 1e-14
            row = (f"| {label} | {cohort} | {r['species']} | {observed:.6f} | "
                   f"{median_null:.6f} | {excess:.6f} | {r['structured_null_upper_p']:.3f} |")
            assert row in supplement
    assert result['structured_null_replicates'] == 999
    assert '**not reapplied**' in supplement
    assert 'not uncertainty intervals' in supplement


def test_h3_sensitivity_tables_match_each_cohort_and_scenario():
    supplement = (ROOT / 'docs/FCP_H123_SUPPLEMENT.md').read_text(encoding='utf-8')
    base = ROOT / 'docs/supporting/h123'
    with (base / 'h3a_signal_by_scenario.csv').open(newline='') as f:
        rows = list(csv.DictReader(f))
    assert {(r['cohort'], r['scenario']) for r in rows} == {
        (c, s) for c in ('discovery', 'reserve') for s in ('S1', 'S2', 'S3')}
    for r in rows:
        values = ' | '.join(f"{float(r[k]):.{digits}f}" for k, digits in [
            ('K_raw', 6), ('p_K_raw', 4), ('lambda_raw', 6), ('p_lambda0', 8),
            ('K_opportunity_residual', 6), ('p_K_opportunity_residual', 4)])
        assert f"| {r['cohort']} | {r['scenario']} | {r['n_tips']} | {values} |" in supplement
    with (base / 'h3b_span_summary.csv').open(newline='') as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        for label, rho, p in [('Raw D', 'rho_D_span', 'p_D_span'),
                              ('Corrected D', 'rho_Dunbiased_span', 'p_Dunbiased_span'),
                              ('Partial ranks', 'partial_rho_D_span', 'p_partial_D_span')]:
            assert (f"| {r['cohort']} | {r['n_species']} | {label} | "
                    f"{float(r[rho]):.6f} | {float(r[p]):.6f} |") in supplement
    assert 'Non-support does not establish equivalence' in supplement


def test_h3_coverage_matches_frozen_tree_manifest():
    manifest = json.loads((ROOT / 'results/polymorphism_h3a_phylogeny_preflight_20260912/frozen_tree_manifest.json').read_text())
    manuscript = (ROOT / 'docs/FCP_H123_MANUSCRIPT.md').read_text(encoding='utf-8')
    for r in manifest['cohorts'].values():
        assert f"{r['retained_tips_each_scenario']}/{r['eligible_species']}" in manuscript
        assert r['retained_tips_each_scenario'] / r['eligible_species'] >= .9
    assert 'unmatched species are not assigned zero signal' in manuscript


def test_h1_agreement_table_matches_direct_analysis_only():
    result = json.loads((ROOT / 'results/polymorphism_h1_observer_disjoint_d_20260913/result.json').read_text())
    supplement = (ROOT / 'docs/FCP_H123_SUPPLEMENT.md').read_text(encoding='utf-8')
    for cohort in ('discovery', 'reserve'):
        for outcome in ('D', 'D_unbiased'):
            r = result['analyses'][cohort][outcome]
            values = ' | '.join(f'{r[k]:.6f}' for k in ('spearman_rho', 'lin_ccc', 'median_abs_difference', 'q90_abs_difference'))
            assert f"| {cohort} | {outcome} | {r['n']} | {values} |" in supplement
            assert r['bootstrap_replicates_finite'] == 5000
            assert r['permutation_replicates'] == 20000
    manuscript = (ROOT / 'docs/FCP_H123_MANUSCRIPT.md').read_text(encoding='utf-8')
    assert 'within species, not as a global partition' in manuscript
    assert 'not an\nobserver-cluster bootstrap' in manuscript


def test_initial_statistical_citations_have_bibliography_and_access_audit():
    manuscript = (ROOT / 'docs/FCP_H123_MANUSCRIPT.md').read_text(encoding='utf-8')
    audit = (ROOT / 'docs/FCP_H123_STATISTICAL_REFERENCES.md').read_text(encoding='utf-8')
    prose, references = manuscript.split('## References — verified initial set')
    for citation, doi in [('Lin, 1989', '10.2307/2532051'),
                          ('Phipson & Smyth, 2010', '10.2202/1544-6115.1585'),
                          ('Blomberg et al., 2003', '10.1111/j.0014-3820.2003.tb00285.x'),
                          ('Jin & Qian, 2022', '10.1016/j.pld.2022.05.005')]:
        assert citation in prose
        assert references.count(doi) == 1
        assert doi in audit
    assert 'were not fully text-audited' in manuscript
    assert 'does not establish the validity' in prose
