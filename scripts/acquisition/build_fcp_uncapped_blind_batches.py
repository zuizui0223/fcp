#!/usr/bin/env python3
"""Freeze ALL completed metadata photos into the original 3-field blind ROI worker interface.

Never drops photos based on taxon photo count, apparent morph or site climate.
Image URL and species/location/observer identities are hidden from colour worker.
This does not open image pixels and does not authorize biological execution.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sqlite3
from contextlib import ExitStack
from pathlib import Path

LEGACY_PHOTO_SHA256 = {
    "discovery": "ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4",
    "validation": "0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6",
    "prospective": "57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186",
}

FIELDS = ("inat_taxon_id", "species", "observation_id", "photo_id", "photo_url_large",
          "photo_license", "attribution", "observer_id", "observed_on", "latitude", "longitude", "positional_accuracy_m")
WORKER = ("measurement_id", "image_filename", "photo_license")
ACQUISITION = ("measurement_id", "image_filename", "photo_url_large", "photo_license")
SALT = "FCP_ALL_OBSERVABLE_PHOTO_BLIND_20261014_V1"


def new_id(photo_id: str) -> str:
    b=f"{SALT}\x1fphoto\x1f{photo_id}".encode()
    return "FCPR-"+hashlib.sha256(b).hexdigest().upper()[:24]


def freeze_inputs(directories: list[Path], output: Path, *, partitions: int = 128,
                  already_measured: dict[str, Path] | None = None) -> dict:
    if partitions < 1:
        raise ValueError("partitions must be positive")
    output.mkdir(parents=True, exist_ok=True)
    counter = 0
    sources = []
    original_checksums = {}
    previously_measured = 0
    repeated_checkpoint_rows = 0
    already_measured_photo_hits = 0
    # Disk-backed distinct-photo audit keeps the acquisition memory bounded even
    # when one iNat species has millions of photos. No photo count limit.
    db_path = output / "_photo_identity_staging.sqlite"
    if db_path.exists():
        db_path.unlink()
    try:
        with sqlite3.connect(db_path) as db, ExitStack() as stack:
            db.execute("CREATE TABLE ids (photo_id TEXT PRIMARY KEY, row_fingerprint TEXT)")
            # Exclusions preserve ORIGINAL measured data as authority. Photo IDs
            # already measured must not be downloaded / inferred for a second time.
            for label, old_path in (already_measured or {}).items():
                if label not in LEGACY_PHOTO_SHA256:
                    raise ValueError(f"Unknown frozen measured source label: {label}")
                digest = hashlib.sha256(old_path.read_bytes()).hexdigest()
                if digest != LEGACY_PHOTO_SHA256[label]:
                    raise ValueError(f"Frozen {label} measured source SHA256 mismatch")
                original_checksums[label] = digest
                with old_path.open("r", encoding="utf-8", newline="") as old_handle:
                    reader = csv.DictReader(old_handle)
                    if not {"photo_id", "inat_taxon_id"}.issubset(reader.fieldnames or []):
                        raise ValueError("Frozen measured source has no photo identity")
                    for row in reader:
                        photo_id = str(row["photo_id"]).strip()
                        if not photo_id:
                            raise ValueError("Frozen measured source missing photo ID")
                        try:
                            db.execute("INSERT INTO ids(photo_id,row_fingerprint) VALUES (?,NULL)", (photo_id,))
                        except sqlite3.IntegrityError as exc:
                            raise ValueError("Original measured sources overlap, cannot double count") from exc
                        previously_measured += 1
            workers = []
            acq_partitions=[]
            for part in range(partitions):
                fh=stack.enter_context((output / f"worker_partition_{part:03}.csv").open("w",encoding="utf-8",newline=""))
                writer=csv.DictWriter(fh,fieldnames=WORKER)
                writer.writeheader()
                workers.append(writer)
                af=stack.enter_context((output / f"acquisition_partition_{part:03}.csv").open("w",encoding="utf-8",newline=""))
                aw=csv.DictWriter(af,fieldnames=ACQUISITION)
                aw.writeheader()
                acq_partitions.append(aw)
            acq_file=stack.enter_context((output/'acquisition_key.csv').open('w',encoding='utf-8',newline=''))
            meta_file=stack.enter_context((output/'metadata_join_key.csv').open('w',encoding='utf-8',newline=''))
            acq_writer=csv.DictWriter(acq_file,fieldnames=ACQUISITION)
            acq_writer.writeheader()
            join_fields=('measurement_id', *[c for c in FIELDS if c!='photo_url_large'])
            join_writer=csv.DictWriter(meta_file,fieldnames=join_fields)
            join_writer.writeheader()
            for directory in sorted(directories):
                checkpoint=json.loads((directory/'checkpoint.json').read_text(encoding='utf-8'))
                if checkpoint.get('status')!='complete' or checkpoint.get('photo_cap') is not None or checkpoint.get('colour_pixels_opened') is not False:
                    raise ValueError(f"Partial/nonblind photo metadata collection: {directory}")
                source=directory/'photo_manifest.ndjson'
                if not source.exists():
                    raise ValueError(f"Missing completed photo metadata source: {directory}")
                digest=hashlib.sha256()
                source_rows=0
                with source.open('rb') as fh:
                    for raw in fh:
                        digest.update(raw)
                        if not raw.strip():
                            continue
                        v=json.loads(raw)
                        if not set(FIELDS).issubset(v):
                            raise ValueError("Missing source identity/coordinate fields")
                        if int(v['inat_taxon_id']) != int(checkpoint['taxon_id']):
                            raise ValueError('Source taxon changed within completed photo ledger')
                        photo_id=str(v['photo_id'])
                        fingerprint=hashlib.sha256(json.dumps(v,sort_keys=True,ensure_ascii=False).encode()).hexdigest()
                        old=db.execute("SELECT row_fingerprint FROM ids WHERE photo_id=?", (photo_id,)).fetchone()
                        if old is not None:
                            if old[0] is None:  # frozen original, already measured
                                already_measured_photo_hits += 1
                                continue
                            if old[0] == fingerprint:  # retry after checkpoint write interrupted
                                repeated_checkpoint_rows += 1
                                continue
                            raise ValueError("Repeated photo ID with conflicting source metadata")
                        db.execute("INSERT INTO ids(photo_id,row_fingerprint) VALUES (?,?)",(photo_id,fingerprint))
                        measurement=new_id(photo_id)
                        part=int(hashlib.sha256(measurement.encode()).hexdigest()[:16],16)%partitions
                        filename=measurement+'.jpg'
                        workers[part].writerow(dict(measurement_id=measurement,
                                                    image_filename=filename,photo_license=v['photo_license']))
                        acquisition_row=dict(measurement_id=measurement,
                                             image_filename=filename,
                                             photo_url_large=v['photo_url_large'],
                                             photo_license=v['photo_license'])
                        acq_writer.writerow(acquisition_row)
                        acq_partitions[part].writerow(acquisition_row)
                        join_writer.writerow({'measurement_id':measurement,
                                              **{k:v.get(k) for k in FIELDS if k!='photo_url_large'}})
                        counter += 1
                        source_rows += 1
                sources.append({'taxon_id':checkpoint['taxon_id'],
                                'metadata_bytes_sha256':digest.hexdigest(),
                                'metadata_raw_rows':source_rows})
            db.commit()
    finally:
        # Temporary identity DB only; no extra duplicate of photo URL is saved.
        for suffix in ('','-shm','-wal'):
            path=Path(str(db_path)+suffix)
            if path.exists():
                path.unlink()
    summary={"schema":"fcp_uncapped_blind_photo_firewall_v1","photo_cap":None,
             "eligible_total_photo_ids":counter,"partitions":partitions,
             "source_species_manifests":len(directories),"metadata_sources":sources,
             "n_original_photo_ids_preloaded":previously_measured,
             "n_new_metadata_photo_ids_already_in_original":already_measured_photo_hits,
             "duplicate_retry_rows_reconciled":repeated_checkpoint_rows,
             "original_measured_sha256":original_checksums,
             "statuses":"METADATA_FROZEN_PIXEL_OPENING_NOT_AUTHORIZED",
             "worker_visible_fields":list(WORKER),"no_colour_or_environment_selection":True,
             "all_completed_photo_ids_included":True}
    (output/'firewall_receipt.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    return summary


def main(argv:list[str]|None=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--metadata-dir",action="append",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--partitions",type=int,default=128)
    p.add_argument("--already-measured",action="append",default=[],
                   help="Frozen discovery/validation/prospective=path; skip those photo IDs")
    args=p.parse_args(argv)
    old={}
    for item in args.already_measured:
        if "=" not in item:raise ValueError("--already-measured expects label=path")
        label,path=item.split("=",1)
        if label in old:raise ValueError("Duplicated original measurement source")
        old[label]=Path(path)
    print(json.dumps(freeze_inputs(args.metadata_dir,args.output,partitions=args.partitions,
                                   already_measured=old)))
    return 0


if __name__=='__main__':
    raise SystemExit(main())