"""Exercise the exact inherited reassembler with artificial IDs only.

Requires the recorded Git object locally; never fetches or opens real results.
This is a component audit, not execution of the P500 wrapper or H2 test.
"""
from pathlib import Path
import subprocess
import sys
import types

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SOURCE = '9fae6ccdf684a46026f72ba12e98de2c5c54bf2a'
PATH = 'fcp_pipeline/photo_first_measurement_execution.py'
BLOB = 'b6e76efb4db8e7529b5ce5a91c60b6ae68637fb5'


def main():
    actual = subprocess.check_output(
        ['git', 'rev-parse', f'{SOURCE}:{PATH}'], cwd=ROOT, text=True).strip()
    if actual != BLOB:
        raise RuntimeError('Inherited source identity mismatch')
    source = subprocess.check_output(['git', 'show', f'{SOURCE}:{PATH}'], cwd=ROOT)
    module = types.ModuleType('_fcp_exact_reassembly_audit')
    sys.modules[module.__name__] = module
    exec(compile(source, f'{SOURCE}:{PATH}', 'exec'), module.__dict__)

    ids = [f'ARTIFICIAL-{i:05d}' for i in range(49900)]
    worker = pd.DataFrame({'measurement_id': ids})
    key = pd.DataFrame({'measurement_id': ids,
                        'species': [f'ARTIFICIAL_SPECIES_{i // 100}' for i in range(49900)]})
    rows = pd.DataFrame({'measurement_id': ids,
                         'measurement_status': ['image_acquisition_failed'] * 49900,
                         'morph': ['mixed_uncertain'] * 49900})
    partitions = [rows.iloc[p::256].copy() for p in range(256)]

    def assemble(parts=partitions, workers=worker, metadata=key):
        return module.reassemble_complete_measurement(
            parts, workers, metadata, expected_partition_receipts=256)

    result = assemble()
    assert len(result.joined_photos) == 49900
    assert result.result_manifest['mixed_uncertain_rows'] == 49900
    assert result.result_manifest['classified_rows'] == 0
    print('PASS: all 49,900 artificial failure rows retained, zero biological classifications')

    def rejects(label, fn, expected):
        try:
            fn()
        except ValueError as error:
            assert expected in str(error), (label, str(error))
            print(f'PASS: {label} rejected')
        else:
            raise AssertionError(f'{label} was accepted')

    rejects('missing partition', lambda: assemble(partitions[:-1]), 'incomplete_measurement_partitions')
    missing = [p.copy() for p in partitions]
    missing[0] = missing[0].iloc[1:]
    rejects('missing terminal ID', lambda: assemble(missing), 'incomplete_terminal_measurement_coverage')
    duplicate = [p.copy() for p in partitions]
    duplicate[0] = pd.concat([duplicate[0], duplicate[1].iloc[:1]], ignore_index=True)
    rejects('duplicate terminal ID', lambda: assemble(duplicate), 'duplicate IDs')
    foreign = [p.copy() for p in partitions]
    foreign[0].iloc[0, 0] = 'ARTIFICIAL-FOREIGN'
    rejects('foreign terminal ID', lambda: assemble(foreign), 'incomplete_terminal_measurement_coverage')
    rejects('duplicate worker ID', lambda: assemble(workers=pd.concat([worker, worker.iloc[:1]])),
            'duplicate measurement IDs')
    rejects('duplicate metadata ID', lambda: assemble(metadata=pd.concat([key, key.iloc[:1]])),
            'duplicate IDs')
    rejects('missing metadata ID', lambda: assemble(metadata=key.iloc[1:]), 'does not match worker denominator')

    bad = rows.iloc[:1].copy()
    bad['measurement_status'] = 'ARTIFICIAL_UNKNOWN_STATUS'
    rejects('unknown status at sealing',
            lambda: module.validate_terminal_partition_results(bad, bad.measurement_id.tolist()),
            'unknown morph/status')
    # Characterize the dependency boundary, not a safe behavior to rely on.
    unvalidated = [p.copy() for p in partitions]
    unvalidated[0].iloc[0, 1] = 'ARTIFICIAL_UNKNOWN_STATUS'
    accepted = assemble(unvalidated)
    assert accepted.result_manifest['terminal_status_counts']['ARTIFICIAL_UNKNOWN_STATUS'] == 1
    print('LIMITATION: aggregate helper trusts prior sealing status validation; does not repeat it')
    print(f'SOURCE: {SOURCE}:{PATH} blob={BLOB}; synthetic component audit only')


if __name__ == '__main__':
    main()
