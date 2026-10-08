"""Source-backed FCP mechanism/photo crosswalk failure modes."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SOURCE=Path(__file__).resolve().parents[1]/"scripts/analysis/run_polymorphism_mechanism_photo_overlap_20261008.py"
sp=importlib.util.spec_from_file_location("fcp_mechanism_bridge",SOURCE)
m=importlib.util.module_from_spec(sp)
sp.loader.exec_module(m)


def registry(n=8):
    return {
        "schema":"fcp_morph_maintenance_primary_study_registry_v1",
        "systems":[{
            "species":f"Genus{i} species",
            "family":"Syntheticaceae",
            "mechanism":"unit_test_not_biology",
            "evidence_role":"synthetic_phenotypic_observation",
            "source_doi":[f"10.1234/synthetic{i}"],
            "demonstrated":["synthetic outcome"],
            "unresolved":["genetics"],
        } for i in range(n)]
    }


def synthetic_photo(species,n=60,white=20,cohort="discovery",distance=False):
    morph=np.array(["white"]*white+["red_pink"]*(n-white),dtype=object)
    lat=np.linspace(10,10.001,n) if not distance else np.r_[np.full(n//2,10.),np.full(n-n//2,12.)]
    return pd.DataFrame({
        "species":[species]*n,
        "species_binomial":[m.canonical_binomial(species)]*n,
        "inat_taxon_id":np.repeat(100,n),
        "photo_id":np.arange(n),
        "morph":morph,
        "global_classifiable":[True]*n,
        "coordinate_valid":[True]*n,
        "latitude":lat,"longitude":np.zeros(n),
        "observer_id":["u"]*n,
        "observed_on":["2020-04-01"]*n,
        "classifiable":[True]*n,
    })


def test_canonical_scientific_binomial_not_substring():
    assert m.canonical_binomial("Silene littorea Brot.")=="Silene littorea"
    assert m.canonical_binomial("Silene littorea")=="Silene littorea"
    assert m.canonical_binomial("Silene littorea var. X")=="Silene littorea"
    assert m.canonical_binomial("Silene")==""
    assert m.canonical_binomial("Silene gallica")!="Silene littorea"


def test_registry_source_identity_quality(tmp_path):
    r=registry()
    f=tmp_path/"registry.json"
    f.write_text(json.dumps(r))
    result=m.registry_check(f)
    assert len(result["systems"])==8
    r["systems"][1]["species"]=r["systems"][0]["species"]
    f.write_text(json.dumps(r))
    with pytest.raises(ValueError,match="Duplicate biological system"):
        m.registry_check(f)


def test_photo_morph_gates_do_not_infer_pigment_genotype():
    sample=synthetic_photo("Genus0 species",n=60,white=12)
    s=m.one_species_photo_summary(sample)
    assert s["meets_40_photo_eligibility"]
    assert s["meets_5_white_and_5_nonwhite"]
    assert s["meets_5_white_and_one_specific_nonwhite_hue"]
    assert not s["meets_two_nonwhite_hues_each_5"]
    assert s["photo_neighborhood_50km_pairs"]==1770
    assert s["photo_50km_white_nonwhite_pairs"]==12*48
    assert not s["photo_white_state_chemically_confirmed"]
    assert not s["in_pop_fitness_genotype_mapped_to_source_images"]


def test_geographic_separation_can_create_two_morphs_without_local_mix():
    d=synthetic_photo("Genus0 species",n=60,white=30,distance=True)
    s=m.one_species_photo_summary(d)
    assert s["meets_5_white_and_5_nonwhite"]
    assert s["photo_50km_white_nonwhite_pairs"]==0
    assert s["photo_neighborhood_50km_pairs"]==2*(30*29//2)


def test_eligibility_and_unknown_photo_source_are_separated():
    old=synthetic_photo("Genus0 species",n=39,white=12)
    s=m.one_species_photo_summary(old)
    assert s["observed_in_photo_candidate_cohort"]
    assert not s["meets_40_photo_eligibility"]
    assert not s["meets_5_white_and_5_nonwhite"]
    assert s["photo_neighborhood_50km_pairs"]==0


def test_exact_registry_join_never_fills_missing_fitness_or_synonyms():
    r=registry()
    d1=pd.concat([
        synthetic_photo("Genus0 species",white=20),
        synthetic_photo("Genus1 species",white=1),
        synthetic_photo("Genus2 other",white=30),
    ],ignore_index=True)
    d2=synthetic_photo("Genus3 species",white=30,cohort="validation")
    out,table=m.join_panel(r,{"discovery":d1,"validation":d2,
                              "third":pd.DataFrame(columns=d1.columns)})
    assert out["n_independent_literature_species"]==8
    assert out["n_literature_species_any_measured_photo_cohort"]==3
    assert out["n_literature_species_ge40_classifiable"]==3
    assert out["n_literature_species_ge5_white_and_ge5_nonwhite"]==2
    assert out["n_true_morph_genotype_chemistry_photo_links"]==0
    assert not out["direct_biochemical_genotype_fitness_bridge_estimable"]
    assert "Genus2 species" not in out["source_species_any_sample"]
    assert len(table)==3


def test_species_disjoint_source_strict_fail_closed():
    reg=registry()
    d=synthetic_photo("Genus0 species")
    with pytest.raises(RuntimeError,match="species-disjoint sampling broken"):
        m.join_panel(reg,{"discovery":d,"validation":d.copy(),
                          "third":pd.DataFrame(columns=d.columns)})
