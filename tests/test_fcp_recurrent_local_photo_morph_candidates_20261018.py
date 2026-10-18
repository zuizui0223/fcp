"""Photo recurrence is an audit lead, never proof of genetically maintained FCP."""
import numpy as np
import pandas as pd
import pytest

from scripts.analysis import fcp_recurrent_local_photo_morph_candidates_20261018 as m


def photos(*, same_month=True, four_observers=True, repeat_hues=True, purity=0.85):
    # Four actual distinct original iNat observations from two spring seasons.
    states=["white","red_pink","white","red_pink" if repeat_hues else "blue_purple"]
    dates=["2021-05-12","2021-05-16","2023-05-08",
           "2023-05-22" if same_month else "2023-06-22"]
    obs=[101,102,103,104] if four_observers else [101,102,101,102]
    records=[]
    for i,morph in enumerate(states):
        row={"inat_taxon_id":7,"species":"Testflower example","photo_id":1000+i,
             "observation_id":2000+i,"observer_id":obs[i],
             "latitude":35.0+i*.0004,"longitude":135.0,
             "positional_accuracy_m":9,"observed_on":dates[i],
             "morph":morph,"global_classifiable":True}
        for components in m.PALETTE.values():
            for col in components:
                row[col]=0.0
        component=m.PALETTE[morph][0]
        row[component]=purity
        row["flower_fraction_purple"]=row.get("flower_fraction_purple",0.)+(1-purity)
        records.append(row)
    return pd.DataFrame(records)


def test_repeated_cross_year_hue_pairs_with_same_calendar_month():
    d=m.read_and_validate_from_dataframe(photos()) if hasattr(m,"read_and_validate_from_dataframe") else photos()
    # Prepare the synthetic exactly as real read_and_validate does.
    d["year"]=pd.to_datetime(d.observed_on).dt.year.astype("int16")
    d["month"]=pd.to_datetime(d.observed_on).dt.month.astype("int8")
    d["year_month"]=d.year.astype("int32")*12+d.month
    for morph,cols in m.PALETTE.items():
        d["fraction_"+morph]=d[list(cols)].sum(axis=1)
    f=d[["fraction_"+x for x in m.MORPHS]].to_numpy()
    mask=np.column_stack([d.morph.to_numpy()!=x for x in m.MORPHS])
    main=~mask
    d["primary_fraction"]=np.sum(np.where(main,f,0.),axis=1)
    d["max_other_fraction"]=np.max(np.where(mask,f,-1.),axis=1)
    result,table=m.analyse(d)
    assert result["precise_coordinates_and_source_observation_identifiers_publicly_withheld"]
    assert result["radius_sensitivity"]["1.0"]["interannual_recurrent_four_observer_species"]==1
    assert table.loc[table.radius_km==1,"some_four_photo_event_purity_ge_65"].all()
    assert table.loc[table.radius_km==1,"calendar_month"].iloc[0]==5
    assert table.loc[table.radius_km==1,"coarse_colour_pair"].iloc[0]=="red_pink / white"
    assert not {"latitude","longitude","photo_id","observation_id","observer_id"}.intersection(table.columns)


def test_two_years_not_an_independent_population_census():
    d=photos(four_observers=False)
    d["year"]=pd.to_datetime(d.observed_on).dt.year
    d["month"]=pd.to_datetime(d.observed_on).dt.month
    d["year_month"]=d.year*12+d.month
    d["primary_fraction"]=.9
    d["max_other_fraction"]=.1
    xyz=m.xyz_from_lat_lon(d)
    _,_,u,v=m.strict_colour_pairs(d,xyz,1.)
    assert len(u)==2
    assert m.best_recurrence(d,xyz,u,v,1.)==[]


def test_recurrent_colour_pair_requires_same_month_and_same_two_colour_names():
    for option in ("different_month","different_morph"):
        d=photos(same_month=option!="different_month",
                 repeat_hues=option!="different_morph")
        d["year"]=pd.to_datetime(d.observed_on).dt.year
        d["month"]=pd.to_datetime(d.observed_on).dt.month
        d["year_month"]=d.year*12+d.month
        d["primary_fraction"]=.9
        d["max_other_fraction"]=.1
        xyz=m.xyz_from_lat_lon(d)
        _,_,u,v=m.strict_colour_pairs(d,xyz,1.)
        assert m.best_recurrence(d,xyz,u,v,1.)==[]


def test_four_photo_compositional_purity_does_not_misreport_bicolour():
    d=photos(purity=.60)
    d["year"]=pd.to_datetime(d.observed_on).dt.year
    d["month"]=pd.to_datetime(d.observed_on).dt.month
    d["year_month"]=d.year*12+d.month
    d["primary_fraction"]=.60
    d["max_other_fraction"]=.40
    xyz=m.xyz_from_lat_lon(d)
    _,_,u,v=m.strict_colour_pairs(d,xyz,1.)
    r=m.best_recurrence(d,xyz,u,v,1.)
    assert len(r)==1
    assert r[0]["some_four_photo_event_purity_ge_60"]
    assert not r[0]["some_four_photo_event_purity_ge_65"]
    assert not r[0]["some_four_photo_event_alt_le20"]
    assert r[0]["possible_photo_hue_mixing"]


def test_explicit_limits_and_source_requirements():
    assert m.SOURCE_PHOTOS==447500 and m.SOURCE_CLASSIFIED_PHOTOS==227107
    assert m.SOURCE_CSV_SHA256=="8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6"
    assert m.RADIUS_KM==(0.25,0.5,1.0,2.0,5.0)
    assert m.PURITY_THRESHOLDS==(0.50,0.55,0.60,0.65,0.70,0.75)
    assert "NOT evidence" in m.__doc__
    assert "genotypes" in m.__doc__
    with pytest.raises(ValueError,match="source-photo SHA256"):
        m.read_and_validate(__file__)


def test_real_recurrence_receipt_keeps_sample_gates_and_population_nonclaims():
    import json
    from pathlib import Path
    root=Path(__file__).resolve().parents[1]
    d=json.loads((root/"results"/"fcp_fivecohort_recurrent_local_morph_candidates_20261018.json").read_text())
    assert d["schema"]=="fcp_fivecohort_recurrent_observational_colour_contrast_audit_receipt_v1"
    assert d["original_measured_CSV_sha256"]==m.SOURCE_CSV_SHA256
    assert d["total_genuine_original_photo_terminals"]==447500
    assert d["four_state_classifiable_photo_terminals"]==227107
    assert d["exact_reproduction_gitHub_run_id"]==39097882762
    assert d["true_source_workflow_success"] is True
    assert d["source_run_output_JSON_sha256"]=="f5ac0c5e1b02bdec88275fc9eaecaef6f703610855470b46ab26f4610a273bce"
    assert d["source_run_anonymized_species_CSV_sha256"]=="2899742b2bb737d9f36dd1f7907bb0257de1ba3b722530a949ebb598d1e47fae"
    s={x["radius_km"]:x for x in d["opportunity"]}
    assert [s[r]["n_species_with_ge1_contrasting_photo_pair"] for r in (.25,.5,1,2,5)]==[28,49,91,145,300]
    assert [s[r]["n_species_with_four_observer_recurrence"] for r in (.25,.5,1,2,5)]==[2,3,8,8,17]
    assert s[1]["n_species_with_all_four_photos_primary_fraction_ge65pct"]==1
    assert s[5]["n_species_with_all_four_photos_primary_fraction_ge65pct"]==3
    assert d["quality_purity_threshold_species_sensitivity_at_1km"] == {
       "50":7,"55":6,"60":2,"65":1,"70":1,"75":0
    }
    assert len(d["one_km_interannual_species"])==8
    assert len(d["five_km_additional_species"])==9
    assert all(x["year_examples"][0] < x["year_examples"][1] for x in d["one_km_interannual_species"])
    assert d["no_exact_coordinates_photo_ids_observation_ids_or_observer_ids_in_source_summary"] is True
    for t in d["warnings"]:
        assert isinstance(t,str) and len(t)>12
    text=(root/"docs"/"FCP_RECURRENT_LOCAL_PHOTO_MORPH_QA_AND_MAINTENANCE_IDENTIFIABILITY_20261018.md").read_text()
    for phrase in ("8 unverified taxa","Collinsia multicolor",
                   "Castilleja rubicundula","0.300","0.275",
                   "parameter nonidentifiability","not a demarcated reproductive population"):
        assert phrase.lower() in text.lower(),phrase
