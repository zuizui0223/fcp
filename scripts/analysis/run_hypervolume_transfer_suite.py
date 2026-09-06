#!/usr/bin/env python3
"""Run or resume the 57 frozen hypervolume benchmark arms, then validate summaries.

Existing arm files are retained, not overwritten. Use a new output directory for
an independent reproduction. This driver does not change scientific settings.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import os
from pathlib import Path
import subprocess
import sys


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--geometry', type=Path, help='Optional geometry-only CSV with adjacent audit.json')
    args = parser.parse_args()
    if not 1 <= args.workers <= 4:
        parser.error('--workers must be in 1..4')
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    script = Path(__file__).with_name('run_hypervolume_transfer_benchmark.py')
    nuisance = ('geographic_specific', 'environmental_specific', 'mixed_specific', 'no_structure')
    jobs: list[list[str]] = []
    for mapping in ('harmonic', 'confounded'):
        evaluation = (('geographic_shared', 'environmental_shared') if mapping == 'harmonic' else ('confounded_shared',)) + nuisance
        for amplitude in (0.5, 1.0, 2.0):
            for stage, worlds in (('calibration', nuisance), ('evaluation', evaluation)):
                for world in worlds:
                    destination = out / f'{mapping}__{stage}__{world}__a{amplitude:g}.csv'
                    if not destination.exists():
                        jobs.append([sys.executable, str(script), '--output-dir', str(out), '--mapping', mapping,
                                     '--stage', stage, '--world', world, '--amplitude', str(amplitude)])
    env = os.environ.copy()
    for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
        env[key] = '1'

    def execute(command: list[str]) -> None:
        subprocess.run(command, env=env, check=True)

    print(f'Missing frozen arms: {len(jobs)}; existing outputs retained.', flush=True)
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        list(pool.map(execute, jobs))
    execute([sys.executable, str(script), '--output-dir', str(out), '--summarize'])
    if args.geometry:
        execute([sys.executable, str(script), '--output-dir', str(out), '--geometry', str(args.geometry.resolve())])
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
