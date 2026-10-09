"""Synthetic original direct LCVP tip phylogenetic spatial coverage, no flower labels."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import Phylo
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
import audit_fcp_global42111_direct_phylo_nearby_photo_pairs_20261010 as M


@pytest.fixture
def source_tree():
    from io import StringIO
    tree=Phylo.read(StringIO("(((G_a:1,G_b:1):1,G_c:2):1,((H_a:1,H_b:1):1,H_c:2):1);"),"newick")
    source=pd.DataFrame({
       "inat_taxon_id":[1,2,3,4,5,6],
       "species":["G a","G b","G c","H a","H b","H c"],
       "original_direct_LCVP_tip":["G_a","G_b","G_c","H_a","H_b","H_c"],
       "genus":["G","G","G","H","H","H"],
       "genus_cell_id":["G|50"]*3+["H|50"]*3,
       "latitude":[0.,0.,0.1,4.,4.,4.1],
       "longitude":[10.,10.2,10.1,20.,20.1,20.2],
    })
    return source,tree


def test_true_backbone_local_pairs_not_artificial_genus_distance(source_tree,monkeypatch):
    source,tree=source_tree
    monkeypatch.setattr(M,"MIN_DIRECT_PHYLO_CONNECTED_TAXA",6)
    monkeypatch.setattr(M,"MIN_CONGENERIC_GENERA",2)
    monkeypatch.setattr(M,"MIN_LOCAL_MICROGROUPS",2)
    d=M.local_coverage(source,tree,50)
    assert d["n_original_direct_backbone_photo_species"]==6
    assert d["n_local_multispecies_direct_tip_groups"]==2
    assert d["n_original_congeneric_phylogenetic_tip_pairs_in_local_neighborhoods"]==6
    assert d["n_original_species_in_at_least_one_local_direct_phylo_pair"]==6
    assert d["n_source_original_congeneric_genera_with_local_direct_pairs"]==2
    assert d["all_phylogenetic_pairs_both_direct_real_backbone_tips"]
    assert d["readiness"]=="EXPLORATORY_DIRECT_TIP_SPATIAL_PHYLO_COVERAGE_PASS"


def test_source_colours_cannot_change_actual_phylogenetic_comparability(source_tree):
    d,tree=source_tree
    before=M.local_coverage(d,tree,100)
    after=d.assign(morph=["white","red_pink"]*3)
    assert M.local_coverage(after,tree,100)==before


def test_exact_name_mapping_and_never_match_genus_only(tmp_path,source_tree,monkeypatch):
    source,tree=source_tree
    original=pd.DataFrame([{"inat_taxon_id":i+1,
        "original_scientific_name":s,
        "in_fixed_250km_original_source":True,
        "direct_lcvp_backbone_tip":True} for i,s in enumerate(source.species)])
    file=tmp_path/"six_direct.tre"
    Phylo.write(tree,str(file),"newick")
    monkeypatch.setattr(M,"KNOWN_EXACT_BACKBONE",{"250":6,"500":6})
    x,phy=M.direct_taxa(source,original,file,250,strict=False)
    assert len(x)==6 and len(phy.get_terminals())==6
    bad=source.copy()
    bad.loc[0,"species"]="G x"
    with pytest.raises(ValueError,match="exact binomial"):
        M.direct_taxa(bad,original,file,250,strict=False)


def test_incomplete_phylogeny_pairs_hold_not_negative(source_tree):
    d,tree=source_tree
    local=M.local_coverage(d,tree,50)
    assert local["readiness"]=="HOLD_INSUFFICIENT_ORIGINAL_DIRECT_PHYLOGENETIC_LOCAL_COMPARISONS"


def test_tree_with_unmatched_photograph_tip_rejected(tmp_path,source_tree,monkeypatch):
    d,tree=source_tree
    ledger=pd.DataFrame([{"inat_taxon_id":i+1,
        "original_scientific_name":s,
        "in_fixed_250km_original_source":True,
        "direct_lcvp_backbone_tip":True} for i,s in enumerate(d.species)])
    tree.get_terminals()[0].name="Other_other"
    f=tmp_path/"bad.tre";Phylo.write(tree,str(f),"newick")
    monkeypatch.setattr(M,"KNOWN_EXACT_BACKBONE",{"250":6,"500":6})
    with pytest.raises(ValueError,match="exact binomial"):
        M.direct_taxa(d,ledger,f,250,strict=False)
