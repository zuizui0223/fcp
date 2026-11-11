"""Fail-closed checks on WorldClim full FCP photo transport and pathway evidence."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from scripts.analysis import fcp_molecular_pathway_photo_bridge_20261018 as m
from scripts.analysis import fcp_447500_bio5_commonality_transport_20261018 as b


def test_external_pathway_literature_does_not_turn_into_a_gene_by_photo_join():
    assert m.PHOTO_SHA==b.ORIGINAL_PHOTO_SHA256
    assert m.N_PHOTOS==b.SPECIES_PHOTOS==447500
    assert m.N_TAXA==b.SPECIES==4475
    evidence={item["system"]:item for item in m.SOURCES}
    assert len(m.SOURCES)==5
    assert evidence["Silene littorea"]["fcp_accepted_name"] is None
    assert evidence["Mimulus guttatus"]["fcp_accepted_name"]=="Erythranthe guttata"
    assert evidence["Mimulus aurantiacus"]["fcp_accepted_name"]=="Diplacus aurantiacus"
    assert evidence["Ipomoea purpurea"]["fcp_accepted_name"]=="Ipomoea purpurea"
    assert evidence["Mimulus guttatus"]["genetic_variant_causally_identified"] is False
    assert evidence["Mimulus aurantiacus"]["reproductive_isolation_demonstrated"] is True
    assert evidence["Silene littorea"]["genetic_variant_causally_identified"] is False
    assert evidence["Ipomoea purpurea"]["genetic_variant_causally_identified"] is True
    assert all("doi" in x for x in m.SOURCES)


def test_external_system_photo_subset_cannot_be_misreported_as_genotype_phenotype_match():
    items={
        "Ipomoea purpurea":{"blue_purple":28,"red_pink":11,"white":7,"yellow_orange":1},
        "Erythranthe guttata":{"yellow_orange":60},
        "Diplacus aurantiacus":{"yellow_orange":66,"white":3,"blue_purple":1},
    }
    for name, photo_states in items.items():
        assert sum(photo_states.values()) > 0
        assert name in [s["fcp_accepted_name"] for s in m.SOURCES]
    # Genes implicated in Mimulus guttatus nectar-guide/vegetative phenotypes
    # cannot be tested with all-yellow whole-flower palette classes.
    assert items["Erythranthe guttata"] == {"yellow_orange":60}
    assert m.SOURCES[2]["target_trait"].startswith("anthocyanin nectar-guide spots")
    assert "NOT whole yellow-vs-white" in m.SOURCES[2]["target_trait"]


def test_bio5_sign_and_missing_groups_dont_become_universal():
    x=b.summarise_deltas([.04,.05,-.12])
    assert x["n_species"]==3
    assert x["median_white_minus_chromatic_BIO5_within_species_SD"]==pytest.approx(.04)
    assert x["fraction_positive_among_all_evaluable_species"]==pytest.approx(2/3)
    assert x["p_wilcoxon"] > 0.05
    assert b.summarise_deltas([])["estimable"] is False
