"""Test direct reproductive traits without substituting self-compatibility proxies."""
import importlib.util
from pathlib import Path

import pandas as pd

SRC=Path(__file__).resolve().parents[1]/"scripts/analysis/run_polymorphism_direct_reproductive_assurance_overlap_20261008.py"
sp=importlib.util.spec_from_file_location("fcp_reproductive_direct_overlap",SRC)
m=importlib.util.module_from_spec(sp)
sp.loader.exec_module(m)


def photos():
    records=[]
    for cohort,name,nwhite in [
        ("discovery","Silene littorea",12),
        ("validation","Moricandia arvensis",30),
        ("third","Dactylorhiza sambucina",0)
    ]:
        records.append({"cohort":cohort,"species":name,
            "canonical_binomial":name,"genus":name.split(" ")[0],
            "inat_taxon_id":len(records)+1,
            "n_classifiable":60,"photo_white":nwhite,
            "photo_chromatic":60-nwhite,
            "white_chromatic_photo_eligible":5<=nwhite<=55,
            "two_chromatic_hues_at_least5":True})
    return pd.DataFrame(records)


def test_binomial_names_do_not_fuzzy_match():
    assert m.binomial("Silene littorea Brot.")=="Silene littorea"
    assert m.binomial("Silene_littorea")=="Silene littorea"
    assert m.binomial("Silene")==""
    assert m.binomial("Silene littoera")!="Silene littorea"
    assert m.binomial("")==""


def test_goodwillie_only_natural_population_direct_outcrossing():
    g=pd.DataFrame({
        "Genus species":["Silene littorea","Moricandia arvensis","Dactylorhiza sambucina"],
        "Mean-tm":["0.35","0.95","0.4"],
        "PopType (0=natural, 1= experimental, 2=seed orchard, 3=agricultural)":[0,1,0]
    })
    names,col=m.get_source_names(g,"goodwillie")
    info=m.direct_trait_info(g,"goodwillie")
    assert col=="Genus species"
    assert names.tolist()[0]=="Silene littorea"
    assert info["measured_mask"].tolist()==[True,False,True]
    assert info["n_natural_measured_rows"]==2


def test_razanajatovo_autofertility_does_not_impute_self_compatibility():
    d=pd.DataFrame({
        "Species":["Silene_littorea","Moricandia_arvensis","Dactylorhiza_sambucina"],
        "Autofertility_index_FS":["0.25","NA","0"],
        "Autofertility_index_SFL":["NA","NA","NA"],
        "Self-compatibility_index_FS":["0.5","0.92","1"],
        "Self-compatibility_index_SFL":["NA","NA","NA"]
    })
    taxa,_=m.get_source_names(d,"razanajatovo")
    info=m.direct_trait_info(d,"razanajatovo")
    assert taxa.iloc[0]=="Silene littorea"
    assert info["measured_mask"].tolist()==[True,False,True]
    assert info["n_self_compatibility_only_rows"]==1


def test_rodger_requires_explicit_exclusion_outcome_and_unambiguous_taxon():
    d=pd.DataFrame({"genus.species":["Silene_littorea","Moricandia_arvensis"],
                    "taxon":["authored name","another synonym"],
                    "auto.fs.x":["0","0.55"]})
    names,col=m.get_source_names(d,"rodger")
    a=m.direct_trait_info(d,"rodger")
    assert col=="genus.species"
    assert set(names)=={"Silene littorea","Moricandia arvensis"}
    assert a["measured_mask"].tolist()==[True,True]
    d["taxon_name"]=d["taxon"]
    x,col=m.get_source_names(d,"rodger")
    assert col=="genus.species"
    assert x.tolist()==names.tolist()
    d=d.drop(columns="genus.species")
    x,col=m.get_source_names(d,"rodger")
    assert col=="HOLD_MISSING_SOURCE_DEFINED_GENUS_SPECIES"
    assert x.eq("").all()
    bad=pd.DataFrame({"auto.fs.x":["0","0.55"],"auto.spfr.x":["-0.02","invalid"]})
    assert not m.direct_trait_info(bad,"rodger")["measured_mask"].any()


def test_direct_trait_overlap_does_not_promote_selection_and_requires_coverage():
    p=photos()
    g=pd.DataFrame({
        "Genus species":["Silene littorea","Moricandia arvensis"],
        "Mean-tm":["0.3","0.7"],
        "PopType (0=natural, 1= experimental, 2=seed orchard, 3=agricultural)":[0,0]
    })
    r=pd.DataFrame({
        "Species":["Silene_littorea","Moricandia_arvensis"],
        "Autofertility_index_FS":["0.2","NA"],
        "Autofertility_index_SFL":["NA","NA"],
        "Self-compatibility_index_FS":["1","1"],
        "Self-compatibility_index_SFL":["NA","NA"]
    })
    x=pd.DataFrame({"genus.species":["Moricandia_arvensis"],"auto.fs.x":[0.4]})
    out,ledger=m.summarize_overlap(p,{"rodger":x,"goodwillie":g,"razanajatovo":r})
    assert out["source_photo_species_total"]==3
    assert out["source_trait_panels"]["goodwillie"]["n_fcp_species_measured_trait_overlap"]==2
    assert out["source_trait_panels"]["rodger"]["n_fcp_white_chromatic_photo_species_measured_trait_overlap"]==1
    assert out["source_trait_panels"]["razanajatovo"]["n_fcp_species_measured_trait_overlap"]==1
    assert out["all_sources_direct_trait_overlap_status"]=="HOLD_NO_SUFFICIENT_INDEPENDENT_SPECIES_COHORT_COVERAGE"
    assert not any(z["photo_colour_trait_association_eligible"] for z in out["source_trait_panels"].values())
    assert (~ledger["source_plant_genotype_linked_to_photo"]).all()


def test_empty_direct_measurement_does_not_count_species():
    x=pd.DataFrame({"genus.species":["Silene_littorea"],"auto.fs.x":["NA"]})
    p=photos()
    r,tab=m.summarize_overlap(p,{"rodger":x})
    assert r["source_trait_panels"]["rodger"]["n_fcp_species_any_named_source_overlap"]==1
    assert r["source_trait_panels"]["rodger"]["n_fcp_species_measured_trait_overlap"]==0
    assert tab.empty
