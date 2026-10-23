"""No-network original Sept-2026 source provenance gate for one REAL-ROI worker."""
from __future__ import annotations
import json
import shutil
import pytest

from scripts.acquisition.freeze_fcp_photo_first_blind_workers import freeze
from scripts.acquisition.plan_rgfca_screen_from_targeted_original_monthly_photo_ids import create
from scripts.acquisition.prepare_first_rgfca_original_monthly_roi_unit import prepare
from test_plan_rgfca_screen_from_targeted_original_monthly_photo_ids import fixture


def source_tree(tmp_path):
    manifest,receipt,frame,counts,audit,roster,cohort_receipt=fixture(tmp_path)
    root=tmp_path/"original_monthly_source_artifact"
    proof=root/"rgfca_1000_exact_source_photo_IDs"
    proof.mkdir(parents=True)
    shutil.copy2(manifest,proof/manifest.name)
    shutil.copy2(receipt,proof/receipt.name)
    source=root/"rgfca_1000_photo_colour_blind_selection"
    create(proof/manifest.name,proof/receipt.name,frame,counts,audit,
           roster,cohort_receipt,source,photos_per_fold=2,
           max_per_observer=2,cell_candidate_limit=3,synthetic=True)
    workers=root/"rgfca_1000_original_photo_blind_worker_freeze"
    frozen=freeze(source/"photo_screen_original_metadata_selected.csv",
                  source/"photo_screen_source_selection_receipt.json",
                  frame,counts,workers,_synthetic=True)
    assert frozen["screen_selected_original_photo_IDs"]>0
    return root


def test_first_original_monthly_ROI64_unit_is_from_admitted_source_only(tmp_path):
    root=source_tree(tmp_path)
    out=tmp_path/"first_true_roi"
    result=prepare(root,out,synthetic=True)
    assert result["status"]=="ORIGINAL_SOURCE_64_PHOTO_ROI_UNIT_READY_NOT_MEASURED"
    assert result["n_rgfca_enrolled_source_taxa"]==2
    assert result["n_rgfca_all_original_source_photo_IDs"]==5
    assert result["n_photo_IDs_in_this_unit"]>=1
    assert result["batch_size"]==64 and result["batch"]==0
    assert 0<=result["partition"]<128
    assert result["original_blind_worker_3_fields_only"] is True
    assert result["original_image_pixels_opened"] is False
    assert result["photo_polymorphism_Ppoly_W_H2_climate_or_genetics_computed"] is False
    assert json.loads((out/"first_original_ROI64_ready_receipt.json").read_text())==result


def test_original_photo_worker_sha_mutation_is_caught_before_first_image(tmp_path):
    root=source_tree(tmp_path)
    workers=root/"rgfca_1000_original_photo_blind_worker_freeze"
    meta=workers/"metadata_join_key.csv"
    meta.write_text(meta.read_text()+"TAMPERED_PHOTOGRAPHER_BIOLOGY\n")
    with pytest.raises(ValueError,match="Immutable original photo worker"):
        prepare(root,tmp_path/"bad",synthetic=True)


def test_fake_exact_dated_source_or_missing_global_collision_proof_rejected(tmp_path):
    root=source_tree(tmp_path)
    receipt=root/"rgfca_1000_exact_source_photo_IDs/rgfca_1000_targeted_monthly_photo_ID_receipt.json"
    claim=json.loads(receipt.read_text())
    claim["all_enrolled_species_match_independent_full_source_photo_ID_counts"]=False
    receipt.write_text(json.dumps(claim))
    with pytest.raises(ValueError,match="No SHA-sealed original"):
        prepare(root,tmp_path/"impostor",synthetic=True)


def test_substitution_with_live_API_photo_ID_manifest_rejected_and_output_immutable(tmp_path):
    root=source_tree(tmp_path)
    out=tmp_path/"original"
    prepare(root,out,synthetic=True)
    with pytest.raises(ValueError,match="Never overwrite"):
        prepare(root,out,synthetic=True)
    source=root/"rgfca_1000_exact_source_photo_IDs/rgfca_1000_complete_original_source_photo_IDs.csv"
    source.write_text(source.read_text().replace("Species alba","Fake red-flowered species"))
    with pytest.raises(ValueError,match="No SHA-sealed original"):
        prepare(root,tmp_path/"wrong",synthetic=True)
