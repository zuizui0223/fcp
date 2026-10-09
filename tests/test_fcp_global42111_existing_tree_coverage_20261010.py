"""Exact tree species-tip matching is not genus matching and never sees photo colours."""
from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd
import pytest
from Bio import Phylo

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
import audit_fcp_global42111_existing_tree_coverage_20261010 as M


def test_source_canonical_identity_without_synonym_imputation():
    assert M.key("Genus_alpha")=="genus alpha"
    assert M.key(" Genus  alpha ")=="genus alpha"
    assert M.key("Genus_alpha_var_minor")!="genus alpha"
    assert M.key("Genus_beta")!="genus alpha"


def test_old_tree_exact_species_tip_counts_not_mere_genus_overlap(tmp_path,monkeypatch):
    fake="((Genus_alpha:1,Genus_beta:1):1,Another_giga:2,Else_else:3,Other_other:3,More_more:3,Plant_plant:3,Echo_echo:3,Eleven_one:3,Eleven_two:3,Eleven_three:3,Eleven_four:3);"
    paths={}
    for label in ("h3a_s1","h3a_s2","h3a_s3","jbi34"):
        p=tmp_path/f"{label}.tre";p.write_text(fake)
        paths[label]=p
    monkeypatch.setattr(M,"TREE_SHA256",{k:M.sha(v) for k,v in paths.items() if k!="jbi34"})
    leaves,receipt=M.trees_from_sources(paths)
    assert len(leaves["h3a_s1"])==11
    assert receipt["h3a_s1"]["n_unique_exact_source_tree_tip_names"]==11
    source=pd.DataFrame({
        "inat_taxon_id":range(1,101),
        "species":["Genus alpha","Genus beta","Genus other"]+[f"Genus sp{i}" for i in range(97)],
        "genus":["Genus"]*100,
        "genus_cell_id":["Genus|100"]*100,
    })
    coverage=M.evaluate_coverage(source,leaves)["h3a_s1"]
    assert coverage["n_exact_original_source_species_tree_tips"]==2
    assert coverage["n_original_source_species_with_genus_name_somewhere_in_tree"]==100
    assert coverage["n_genera_with_at_least_two_exact_tree_tips"]==1
    assert not coverage["exact_tip_coverage_passed_conservative_gate"]


def test_no_morphological_outcome_used_in_tree_coverage():
    d=pd.DataFrame({
        "inat_taxon_id":range(1,103),
        "species":["G alpha","G beta"]+[f"H sp{i}" for i in range(100)],
        "genus":["G","G"]+["H"]*100,
        "genus_cell_id":["G|8","G|8"]+["H|9"]*100,
        "morph":["white"]*102,
    })
    a=M.evaluate_coverage(d,{"test":{"g alpha","g beta","h sp5"}})
    d["morph"]="red_pink"
    b=M.evaluate_coverage(d,{"test":{"g alpha","g beta","h sp5"}})
    assert a==b
    assert a["test"]["n_exact_original_source_species_tree_tips"]==3
    assert a["test"]["n_original_local_genus_cell_groups_with_at_least_two_exact_tips"]==1


def test_ambiguous_duplicate_newick_tips_do_not_pass(tmp_path,monkeypatch):
    p=tmp_path/"h3a_s1.tre"
    p.write_text("(G_alpha:1,'G alpha':1,X_x:2,Y_y:2,Z_z:2,A_a:2,B_b:2,C_c:2,D_d:2,E_e:2,F_f:2);")
    ps={}
    for k in ("h3a_s1","h3a_s2","h3a_s3","jbi34"):
        path=tmp_path/f"{k}.tre"
        path.write_text(p.read_text())
        ps[k]=path
    monkeypatch.setattr(M,"TREE_SHA256",{k:M.sha(v) for k,v in ps.items() if k!="jbi34"})
    with pytest.raises(ValueError,match="Ambiguous"):
        M.trees_from_sources(ps)
