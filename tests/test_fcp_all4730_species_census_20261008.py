"""All U100 4730 species must be allocated exactly once, without colour inputs."""
import importlib.util
from pathlib import Path
import pandas as pd
import pytest

SRC=Path(__file__).resolve().parents[1]/"scripts/analysis/rebuild_fcp_all4730_species_census_20261008.py"
sp=importlib.util.spec_from_file_location("fcp_all_u100",SRC)
m=importlib.util.module_from_spec(sp)
sp.loader.exec_module(m)


def synthetic_objects():
    sp=pd.DataFrame({"inat_taxon_id":range(1,42112),
                     "species":[f"Genus{i} epitheton" for i in range(1,42112)]})
    cap=pd.DataFrame({"inat_taxon_id":range(1,42112),
                      "after_observer_cap":[100]*4730+[0]*(42111-4730),
                      "request_error":[""]*42111})
    ids=list(range(1,4731))
    def sample(ids):
        return pd.DataFrame({"inat_taxon_id":ids,
            "species":[f"Genus{i} epitheton" for i in ids]})
    # First 1000 old; remaining 3730 partition into 500 + 500 + 2000 + 730.
    return sp,cap,sample(ids[1000:]),sample(ids[1000:1500]),sample(ids[1500:2000]),sample(ids[2000:4000]),sample(ids[4000:])


def test_all_4730_exact_partition():
    r=m.group_allocation(*synthetic_objects())
    assert len(r)==4730
    assert r.inat_taxon_id.is_unique
    assert r.allocation_group.value_counts().to_dict()==m.GROUPS
    assert set(r.inat_taxon_id)==set(range(1,4731))
    assert not set(r.columns) & {"morph","flower_colour","H2_W"}


def test_missing_previous_taxon_is_fail_closed():
    objs=list(synthetic_objects())
    p500=objs[3].copy()
    p500.loc[p500.index[0],"inat_taxon_id"]=10000
    objs[3]=p500
    with pytest.raises(RuntimeError,match="Out-of-U100|stale"):
        m.group_allocation(*objs)


def test_repeat_selection_identity_is_blocked():
    objs=list(synthetic_objects())
    p500=objs[3].copy()
    p500.loc[p500.index[0],"inat_taxon_id"]=objs[4].iloc[0].inat_taxon_id
    p500.loc[p500.index[0],"species"]=objs[4].iloc[0].species
    objs[3]=p500
    with pytest.raises(RuntimeError,match="overlap"):
        m.group_allocation(*objs)


def test_no_colour_or_sex_column_used():
    objs=list(synthetic_objects())
    objs[0]["morph"]=["white"]*42111
    objs[1]["any_color_photo"]=["purple"]*42111
    # Raw source may contain extra columns, but the hard-coded output
    # reads only ID/name and metadata cap. Distinguishing source from
    # analysis outcome is a strict contract.
    ledger=m.group_allocation(*objs)
    assert not {"morph","any_color_photo"}.intersection(ledger.columns)


def test_annual_capacity_failure_does_not_substitute_species():
    objs=list(synthetic_objects())
    objs[1].loc[1000,"after_observer_cap"]=99
    with pytest.raises(RuntimeError,match="Out-of-U100|mismatch"):
        m.group_allocation(*objs)
