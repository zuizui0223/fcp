import json
import subprocess
import sys
import tempfile
from pathlib import Path

SCRIPT = Path('scripts/literature/evaluate_v22_calibration_gate.py')
FIELDS = [
    'record_relevance',
    'natural_intraspecific_variation',
    'floral_display_colour',
    'full_text_required',
]


def metric(n=384, raw=1.0, kappa=1.0, labels=None):
    return {
        'n_double_coded': n,
        'raw_agreement': raw,
        'cohen_kappa': kappa,
        'labels': ['a', 'b'] if labels is None else labels,
    }


def run_gate(overrides):
    agreement = {
        'agreement': {field: metric() for field in FIELDS},
        'focal_taxon_text': {'n_double_coded': 384, 'normalized_exact_agreement': 0.8},
    }
    for field, values in overrides.items():
        agreement['agreement'][field].update(values)
    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        src = root / 'agreement.json'
        out = root / 'gate.json'
        src.write_text(json.dumps(agreement), encoding='utf-8')
        subprocess.run(
            [sys.executable, str(SCRIPT), '--agreement', str(src), '--out', str(out)],
            check=True,
            capture_output=True,
            text=True,
        )
        return json.loads(out.read_text(encoding='utf-8'))


def test_not_ready_if_any_field_has_383_records():
    result = run_gate({'record_relevance': {'n_double_coded': 383}})
    assert result['status'] == 'not_ready'


def test_exact_prespecified_thresholds_pass():
    result = run_gate({
        'record_relevance': {'raw_agreement': 0.90, 'cohen_kappa': 0.60},
        'natural_intraspecific_variation': {'raw_agreement': 0.90, 'cohen_kappa': 0.60},
        'floral_display_colour': {'raw_agreement': 0.90, 'cohen_kappa': 0.60},
        'full_text_required': {'raw_agreement': 0.85, 'cohen_kappa': 0.60},
    })
    assert result['status'] == 'pass'


def test_raw_agreement_just_below_threshold_fails():
    result = run_gate({'record_relevance': {'raw_agreement': 0.899999}})
    assert result['status'] == 'fail'
    assert result['checks']['record_relevance']['raw_ok'] is False


def test_kappa_just_below_threshold_fails():
    result = run_gate({'natural_intraspecific_variation': {'cohen_kappa': 0.599999}})
    assert result['status'] == 'fail'
    assert result['checks']['natural_intraspecific_variation']['kappa_ok'] is False


def test_single_category_undefined_kappa_can_pass_with_raw_agreement():
    result = run_gate({
        'record_relevance': {
            'raw_agreement': 1.0,
            'cohen_kappa': None,
            'labels': ['include'],
        }
    })
    assert result['status'] == 'pass'
    assert result['checks']['record_relevance']['kappa_not_estimable_single_category'] is True


def test_multicategory_null_kappa_fails():
    result = run_gate({
        'record_relevance': {
            'raw_agreement': 1.0,
            'cohen_kappa': None,
            'labels': ['include', 'exclude'],
        }
    })
    assert result['status'] == 'fail'
    assert result['checks']['record_relevance']['kappa_not_estimable_single_category'] is False
