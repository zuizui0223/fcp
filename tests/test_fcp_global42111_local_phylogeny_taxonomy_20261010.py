"""Synthetic public family genus taxonomy guards for fixed original FCP 1761 taxa."""
from __future__ import annotations

import sys
from pathlib import Path
import pandas as pd
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"acquisition"))
import freeze_fcp_global42111_local_phylogeny_taxonomy_20261010 as M


def row(taxon_id=3):
    return {"inat_taxon_id":taxon_id,
            "original_scientific_name":f"Genus species{taxon_id}",
            "original_genus":"Genus",
            "genus_cell_id":"Genus|9",
            "photo_cell_162":9,
            "in_250km_cohort":taxon_id<=872}


def tax(taxon_id=3, name=None, family="Familyaceae", genus="Genus",rank="species"):
    return {"id":taxon_id,"rank":rank,"name":name or f"Genus species{taxon_id}",
            "ancestors":[
                {"rank":"kingdom","name":"Plantae"},
                {"rank":"family","name":family},
                {"rank":"genus","name":genus}]}


def test_verified_exact_source_family_and_genus():
    z=M.interpret(row(3),tax(3))
    assert z["status"]=="EXACT_SPECIES_FAMILY_GENUS_READY"
    assert z["current_family"]=="Familyaceae"
    assert z["taxon_genus_matches_original"] is True


def test_synonym_replacement_and_unknown_ancestor_do_not_pass():
    a=M.interpret(row(3),tax(3,name="Other species"))
    assert a["status"]=="CURRENT_TAXONOMY_DRIFT_OR_SUBSPECIES_REVIEW"
    b=M.interpret(row(3),{"id":3,"rank":"species","name":"Genus species3"})
    assert b["status"]=="TAXON_RETURNED_WITHOUT_RANKED_ANCESTOR_TAXONOMY"
    c=M.interpret(row(3),tax(3,family="Familyaceae",genus="Other"))
    assert c["status"]=="CURRENT_TAXONOMY_DRIFT_OR_SUBSPECIES_REVIEW"


def test_unresolved_no_family_does_not_invent_placement():
    missing=M.interpret(row(4),None)
    assert missing["status"]=="NOT_RETURNED"
    assert missing["current_family"] is None
    odd=M.interpret(row(5),{"id":5,"rank":"species","name":"Genus species5",
        "ancestors":[{"rank":"genus","name":"Genus"}]})
    assert odd["status"]=="AMBIGUOUS_OR_UNRESOLVED_FAMILY_GENUS_ANCESTORS"


def test_public_API_bound_before_network():
    with pytest.raises(ValueError,match="30"):
        M.get_taxa(list(range(31)))
    with pytest.raises(ValueError,match="30"):
        M.get_taxa([3,3])


def test_mocked_original_1761_source_taxonomy_is_query_bounded_and_cohort_frozen():
    manifest=pd.DataFrame([row(i) for i in range(1,1762)])
    calls=[]
    def client(ids):
        calls.append(ids)
        return {"results":[tax(i) for i in ids]}
    r,d=M.audit(manifest,client=client,pause=lambda _:None)
    assert r["original_local_500km_species_denominator"]==1761
    assert r["original_local_250km_species_denominator"]==872
    assert r["n_500km_ready_exact_original_species_family_genus"]==1761
    assert r["n_250km_ready_exact_original_species_family_genus"]==872
    assert r["n_current_taxonomy_api_batches"]==59
    assert max(len(q) for q in calls)<=30
    assert len(d)==1761


def test_api_source_errors_do_not_claim_family_zero_or_usable_tip():
    manifest=pd.DataFrame([row(i) for i in range(1,1762)])
    def fail(ids):raise TimeoutError("simulated")
    r,d=M.audit(manifest,client=fail,pause=lambda _:None)
    assert r["n_500km_ready_exact_original_species_family_genus"]==0
    assert r["n_API_error_batches"]==59
    assert set(d.status)=={"API_ERROR_UNRESOLVED"}
