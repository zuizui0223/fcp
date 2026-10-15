"""Byte-verified, memory-bounded positional index for frozen FCP photo workers.

A full source partition can contain enough original photo IDs that repeatedly
reading the whole worker and acquisition CSV for each 64-photo ROI work unit
would require quadratic I/O. One metadata-only linear pass builds 64-bit row
offsets. Each unit then reads only its own immutable 64 rows.

No original pixels, species names, coordinates, old colours, or source outcomes
are opened by the index. The final all-128 global planner remains independent.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from array import array
from pathlib import Path

import pandas as pd

from fcp_pipeline.fcp_allphotos import sha256_file
from fcp_pipeline.fcp_fresh_study import WORKER

ACQUISITION = ('measurement_id', 'image_filename', 'photo_url_large', 'photo_license')
JOIN = ('measurement_id', 'image_filename', 'photo_license')


def _row(raw: bytes, fields: tuple[str, ...]) -> dict[str, str]:
    """Source generated CSV must have exactly one physical line per record."""
    try:
        cells = next(csv.reader([raw.decode('utf-8')], strict=True))
    except (csv.Error, UnicodeDecodeError) as exc:
        raise ValueError('Malformed single-line frozen photo source packet') from exc
    if len(cells) != len(fields) or any('\n' in cell or '\r' in cell for cell in cells):
        raise ValueError('Frozen photo source must have exactly one complete CSV record per line')
    return dict(zip(fields, cells))


class FrozenPartitionSourceIndex:
    """Once-per-partition metadata scan; O(batch-size) photo preflight afterward."""

    def __init__(self, firewall: Path, partition: int):
        freeze = json.loads((firewall / 'fresh_study_freeze.json').read_text(encoding='utf-8'))
        if (freeze.get('status') != 'FRESH_BLIND_METADATA_FROZEN'
                or not isinstance(partition, int)
                or not 0 <= partition < freeze['worker_packets']
                or freeze.get('historical_measured_CSV_read') is not False):
            raise ValueError('Original-photo index requires valid frozen source partition')
        if freeze.get('blind_firewall_sha256') != sha256_file(firewall / 'firewall_receipt.json'):
            raise ValueError('Frozen source firewall changed before indexed photo access')
        self.partition = partition
        self.worker = firewall / f'worker_partition_{partition:03}.csv'
        self.acquisition = firewall / f'acquisition_partition_{partition:03}.csv'
        self.worker_offsets = array('Q')
        self.acquisition_offsets = array('Q')
        seen = set()
        hashes = (hashlib.sha256(), hashlib.sha256())
        with self.worker.open('rb') as w, self.acquisition.open('rb') as a:
            header_w, header_a = w.readline(), a.readline()
            hashes[0].update(header_w)
            hashes[1].update(header_a)
            if (tuple(next(csv.reader([header_w.decode('utf-8')], strict=True))) != WORKER
                    or tuple(next(csv.reader([header_a.decode('utf-8')], strict=True))) != ACQUISITION):
                raise ValueError('Original photo source packet headers changed')
            while True:
                pos_w, pos_a = w.tell(), a.tell()
                line_w, line_a = w.readline(), a.readline()
                if not line_w and not line_a:
                    break
                if not line_w or not line_a:
                    raise ValueError('Source-blind worker and acquisition rows have different counts')
                hashes[0].update(line_w)
                hashes[1].update(line_a)
                worker = _row(line_w, WORKER)
                acq = _row(line_a, ACQUISITION)
                ident = worker['measurement_id']
                if not ident or ident in seen:
                    raise ValueError('Duplicate or missing original-photo source ID in indexed partition')
                seen.add(ident)
                if (any(worker[key] != acq[key] for key in JOIN)
                        or not acq['photo_url_large']):
                    raise ValueError('Original photo URL/worker identity not source-bound')
                self.worker_offsets.append(pos_w)
                self.acquisition_offsets.append(pos_a)
        for path, digest in ((self.worker, hashes[0].hexdigest()),
                             (self.acquisition, hashes[1].hexdigest())):
            if digest != freeze['frozen_file_sha256'][path.name]:
                raise ValueError(f'Frozen original photo source SHA256 changed: {path.name}')
        self.total_photo_ids = len(self.worker_offsets)

    def get_batch(self, partition: int, batch: int, batch_size: int):
        if (partition != self.partition or batch < 0 or batch_size < 1 or
                batch_size > 1024):
            raise ValueError('Indexed worker query outside frozen photo source/compute bounds')
        begin = batch * batch_size
        if begin >= self.total_photo_ids:
            raise ValueError('Requested batch past last source photograph')
        end = min(begin + batch_size, self.total_photo_ids)
        worker_rows, acquisition_rows = [], []
        with self.worker.open('rb') as w, self.acquisition.open('rb') as a:
            for i in range(begin, end):
                w.seek(self.worker_offsets[i])
                a.seek(self.acquisition_offsets[i])
                worker = _row(w.readline(), WORKER)
                acq = _row(a.readline(), ACQUISITION)
                if any(worker[key] != acq[key] for key in JOIN):
                    raise ValueError('Frozen indexed photo ID/URL changed before pixel download')
                worker_rows.append(worker)
                acquisition_rows.append(acq)
        blind = pd.DataFrame(worker_rows, columns=list(WORKER))
        keyed = pd.DataFrame(acquisition_rows, columns=list(ACQUISITION))
        # Keep the worker's 3-column interface completely independent from
        # source URLs; acquisition lookup is never presented to the evaluator.
        selected = blind.merge(keyed, on=list(JOIN), how='left',
                               validate='one_to_one')
        if selected['photo_url_large'].isna().any():
            raise ValueError('Frozen original photo URL disappeared from indexed packet')
        return blind, selected, math.ceil(self.total_photo_ids / batch_size)
