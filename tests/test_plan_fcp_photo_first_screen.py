"""No network or photo pixels: audit 42111-species pre-colour screening design."""
from __future__ import annotations

import csv
import json
import sqlite3
from pathlib import Path

import pytest

from fcp_pipeline.photo_polymorphism_screen import observer_split, sha256
from scripts.acquisition.plan_fcp_photo_first_screen import plan, SELECTION, _choice_hash


def setup_source(tmp_path, *, tamper_stage=False):
    frame = tmp_path / "original_frame.csv"
    names = {11:"Spec alba", 12:"Spec rubra", 13:"Spec zero", 14:"Spec species"}
    with frame.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["inat_taxon_id","species"])
        w.writerows([[i, name] for i, name in names.items()])
    rows = []
    # A: multiple source observations and observers, both folds represented.
    for fold in (0, 1):
        owners = []
        k = 0
        while len(owners) < 2:
            obs = f"observerA_{k}"
            if observer_split(11, obs) == fold:
                owners.append(obs)
            k += 1
        for j in range(4):
            pid = str(1100 + fold*10 + j)
            rows.append((pid, f"photo_{pid}",owners[j%2], 35.0+j, 139.0+j,11))
    # B: one observer: ALL selected to one fold, not falsely negative.
    for j in range(5):
        pid = str(2200+j)
        rows.append((pid,f"photo_{pid}","one_observer",35.1,136.1,12))
    # D: geographically dispersed, separate observations.
    for j in range(8):
        pid = str(4400+j)
        rows.append((pid,f"photo_{pid}",f"owner{j}",-20.0+j,100.0+j,14))
    photo_counts = {ident:sum(x[-1] == ident for x in rows) for ident in names}
    counts = tmp_path / "original_source_photo_counts.csv"
    with counts.open("w",newline="") as f:
        w=csv.writer(f)
        w.writerow(["inat_taxon_id","species","n_unique_source_photo_ids","snapshot_date"])
        w.writerows([[i,names[i],photo_counts[i],"2026-09-27"] for i in names])
    stage = tmp_path / "staged_original_source"
    stage.mkdir()
    db_path = stage / "_working_fcp_snapshot.sqlite"
    with sqlite3.connect(db_path) as db:
        db.execute("CREATE TABLE observations(uuid TEXT PRIMARY KEY, observer TEXT, lat REAL, lon REAL)")
        db.execute("CREATE TABLE photos(photo_id TEXT PRIMARY KEY, uuid TEXT, license TEXT, extension TEXT, photographer TEXT, taxon INTEGER, ambiguous INTEGER)")
        for pid,uuid,observer,lat,lon,taxon in rows:
            db.execute("INSERT INTO observations VALUES (?,?,?,?)", (uuid,observer,lat,lon))
            db.execute("INSERT INTO photos VALUES (?,?,?,?,?,?,?)",
                       (pid,uuid,"cc-by","JPG",observer,taxon,0))
        db.execute("""CREATE INDEX photos_taxon_numeric_order
                      ON photos(taxon, CAST(photo_id AS INTEGER)) WHERE ambiguous=0""")
        db.commit()
    receipt={
        "schema": "fcp_staged_sqlite_photo_universe_v1",
        "status": "SOURCE_SQLITE_STAGED_SHARDS_NOT_YET_EMITTED",
        "source_snapshot_date": "2026-09-27",
        "source_species_frame_sha256": sha256(frame),
        "n_target_taxon_licensed_photos": len(rows),
        "full_source_census_complete": False,
        "full_study_original_images_measured": False,
        "historical_flower_colour_csv_used": False,
        "original_photo_ids_deduplicated_globally": True,
        "verified_original_dated_archive":{"source_archive_sha256":"test-fixture-only-not-production"},
        "sqlite_basename":db_path.name,
        "sqlite_bytes":db_path.stat().st_size,
        "sqlite_sha256":sha256(db_path),
    }
    if tamper_stage:
        receipt["sqlite_sha256"] = "0"*64
    (stage / "fcp_staged_sqlite_source_receipt.json").write_text(json.dumps(receipt))
    return stage, frame, counts, rows


def read_csv(path):
    with path.open("r", newline="") as f:
        return list(csv.DictReader(f))


def test_first_image_screen_uses_actual_source_ids_not_literature_or_color(tmp_path):
    stage, frame, counts, photo_rows = setup_source(tmp_path)
    out = tmp_path / "photo_screen"
    report = plan(stage,frame,counts,out,photos_per_fold=2,
                  cell_candidate_limit=4,max_per_observer=2,_synthetic=True)
    assert report["source_species"] == 4
    assert report["original_source_photo_ids_globally_considered"] == len(photo_rows)
    assert report["outcome_blind_image_selection_no_floral_colour_opened"] is True
    assert report["literature_polymorphism_labels_used"] is False
    assert report["selected_sample_photo_pixels_measured"] is False
    assert report["selected_sample_NOT_full_93million_photo_census"] is True
    chosen = read_csv(out/"photo_screen_original_metadata_selected.csv")
    source_ids = {r[0] for r in photo_rows}
    assert set(x["photo_id"] for x in chosen).issubset(source_ids)
    assert len({x["photo_id"] for x in chosen}) == len(chosen)
    assert all(x["photo_url_large"].endswith(".JPG") for x in chosen)
    assert all(x["source_snapshot_date"]=="2026-09-27" for x in chosen)
    all_taxa = read_csv(out/"photo_screen_all_species_opportunities.csv")
    assert len(all_taxa) == 4
    assert {int(x["inat_taxon_id"]) for x in all_taxa} == {11,12,13,14}
    assert next(x for x in all_taxa if x["inat_taxon_id"]=="13")["total_selected_for_image_screen"] == "0"
    assert all(x["screen_unmeasured_is_not_monomorphism"] == "True" for x in all_taxa)
    for taxon in (11,12,14):
        observers_by_fold={f:set(x["observer_id"] for x in chosen
                                 if x["inat_taxon_id"]==str(taxon) and x["screen_fold"]==str(f))
                           for f in (0,1)}
        assert not (observers_by_fold[0] & observers_by_fold[1])
        assert all(int(x["screen_fold"]) == observer_split(taxon,x["observer_id"])
                   for x in chosen if x["inat_taxon_id"]==str(taxon))


def test_selection_is_deterministic_across_independent_repeats(tmp_path):
    stage,frame,counts,_ = setup_source(tmp_path)
    a=tmp_path/"one"
    b=tmp_path/"two"
    pa=plan(stage,frame,counts,a,photos_per_fold=2,cell_candidate_limit=4,
            max_per_observer=2,_synthetic=True)
    pb=plan(stage,frame,counts,b,photos_per_fold=2,cell_candidate_limit=4,
            max_per_observer=2,_synthetic=True)
    assert pa["photo_selection_csv_sha256"] == pb["photo_selection_csv_sha256"]
    assert (a/"photo_screen_original_metadata_selected.csv").read_bytes() == (
        b/"photo_screen_original_metadata_selected.csv").read_bytes()
    assert _choice_hash(11,"1100") == _choice_hash(11,"1100")
    assert _choice_hash(11,"1100") != _choice_hash(11,"1101")


def test_mismatched_original_source_sqlite_rejected_before_any_colour(tmp_path):
    stage,frame,counts,_=setup_source(tmp_path,tamper_stage=True)
    with pytest.raises(ValueError,match="not fully frozen"):
        plan(stage,frame,counts,tmp_path/"fail",photos_per_fold=2,
             cell_candidate_limit=4,max_per_observer=2,_synthetic=True)


def test_unchanged_original_species_photo_count_required(tmp_path):
    stage,frame,counts,_=setup_source(tmp_path)
    lines=counts.read_text().replace("Spec alba,8,","Spec alba,7,")
    counts.write_text(lines)
    with pytest.raises(ValueError,match="photo counts drifted"):
        plan(stage,frame,counts,tmp_path/"fail",photos_per_fold=2,
             cell_candidate_limit=4,max_per_observer=2,_synthetic=True)


def test_zero_photo_species_not_falsely_declared_monomorphic(tmp_path):
    stage,frame,counts,_=setup_source(tmp_path)
    out=tmp_path/"photo_screen"
    plan(stage,frame,counts,out,photos_per_fold=2,cell_candidate_limit=4,
         max_per_observer=2,_synthetic=True)
    status=next(row for row in read_csv(out/"photo_screen_all_species_opportunities.csv")
                if row["inat_taxon_id"]=="13")
    assert int(status["source_photo_opportunities"])==0
    assert int(status["total_selected_for_image_screen"])==0
    assert status["screen_unmeasured_is_not_monomorphism"]=="True"


def test_no_overwrite_or_untrusted_unlimited_pilot(tmp_path):
    stage,frame,counts,_=setup_source(tmp_path)
    out=tmp_path/"photo_screen"
    plan(stage,frame,counts,out,photos_per_fold=2,cell_candidate_limit=4,
         max_per_observer=2,_synthetic=True)
    with pytest.raises(ValueError,match="Never overwrite"):
        plan(stage,frame,counts,out,photos_per_fold=2,cell_candidate_limit=4,
             max_per_observer=2,_synthetic=True)
    with pytest.raises(ValueError,match="Invalid first-stage"):
        plan(stage,frame,counts,tmp_path/"bad",photos_per_fold=1000,
             max_per_observer=2,_synthetic=True)
