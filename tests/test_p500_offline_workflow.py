from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]


def assert_unique_keys(node):
    if isinstance(node, yaml.MappingNode):
        keys = [key.value for key, _ in node.value]
        assert len(keys) == len(set(keys)), keys
        for _, value in node.value:
            assert_unique_keys(value)
    elif isinstance(node, yaml.SequenceNode):
        for value in node.value:
            assert_unique_keys(value)


def test_recovery_workflow_duplicate_env_regression():
    path = ROOT / '.github/workflows/rgfca-independent-global-axis-metadata-recovery-v2.yml'
    assert_unique_keys(yaml.compose(path.read_text()))
    doc = yaml.safe_load(path.read_text())
    step = next(s for s in doc['jobs']['shard']['steps']
                if s.get('name') == 'Run one deterministic metadata recovery shard v2')
    assert step['env'] == {'PYTHONUNBUFFERED': '1', 'SHARD_INDEX': '${{ matrix.shard }}'}


def test_offline_workflow_has_no_write_permissions_or_acquisition_runner():
    path = ROOT / '.github/workflows/p500-white-control-offline.yml'
    assert_unique_keys(yaml.compose(path.read_text()))
    doc = yaml.safe_load(path.read_text())
    assert doc['permissions'] == {'contents': 'read'}
    assert set(doc['jobs']) == {'offline'}
    steps = doc['jobs']['offline']['steps']
    checkout = next(s for s in steps if s.get('uses', '').startswith('actions/checkout@'))
    assert checkout['with']['persist-credentials'] is False
    commands = '\n'.join(s.get('run', '') for s in steps)
    for forbidden in ('scripts/acquisition/', 'workflow run', 'git push', 'requests.', 'urlopen'):
        assert forbidden not in commands
    assert "report['opening_authorized'] is False" in commands
