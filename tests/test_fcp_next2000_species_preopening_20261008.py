"""No-image prospective selection firewall and source-frozen U100 opportunity tests."""
import csv
import importlib.util
from pathlib import Path

import pytest

SRC=Path(__file__).resolve().parents[1]/"scripts/analysis/freeze_fcp_next2000_species_preopening_20261008.py"
spec=importlib.util.spec_from_file_location("preopen",SRC)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def pool(n=12,start=1):
    return [{"inat_taxon_id":i,"species":f"Genus{i} species",
             "after_observer_cap":100+(i%7)} for i in range(start,start+n)]


def write_p(path,rows,extra=False):
    names=m.POOL_COLUMNS+(["morph"] if extra else [])
    with path.open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=names,lineterminator="\n")
        writer.writeheader()
        for row in rows:
            r=row.copy()
            if extra:r["morph"]="white"
            writer.writerow(r)


def manifest(path,rows):
    with path.open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=m.MANIFEST_COLUMNS,delimiter="\t",lineterminator="\n")
        writer.writeheader()
        for idx,r in enumerate(rows,start=1):
            writer.writerow({"prospective_rank":idx,
                             "inat_taxon_id":r["inat_taxon_id"],"species":r["species"]})


def test_fixed_source_identity_and_capacity_only(tmp_path):
    parent=tmp_path/"p100.csv"
    write_p(parent,pool())
    x=m.read_parent(parent)
    assert len(x)==12
    assert set(x[0])==set(m.POOL_COLUMNS)
    assert min(z["after_observer_cap"] for z in x)>=100
    assert m.SALT=="FCP_NEXT2000_SPATIAL_20261008_FIXED_V1"


def test_plant_morph_or_geography_leak_is_rejected_before_sampling(tmp_path):
    bad=tmp_path/"bad.csv"
    write_p(bad,pool(),extra=True)
    with pytest.raises(ValueError,match="schema mismatch"):
        m.read_parent(bad)
    rows=pool()
    rows[0]["after_observer_cap"]=99
    write_p(bad,rows)
    with pytest.raises(ValueError,match="below-U100 capacity"):
        m.read_parent(bad)


def test_independent_selection_and_reserve_are_disjoint_and_reproducible():
    parent=pool(21)
    used_p500=[parent[0]]
    past_third=[{"prospective_rank":1,**parent[1]}]
    a=m.allocate(parent,used_p500,past_third)
    b=m.allocate(parent,used_p500,past_third)
    assert a==b
    rebuilt,new,reserve=a
    assert len(rebuilt)==20
    assert len(new)==19
    assert not reserve
    assert {r["inat_taxon_id"] for r in new}.isdisjoint({1,2})
    assert sorted(z["prospective_rank"] for z in new)==list(range(1,20))
    assert all(len(z["selection_hash"])==64 for z in new)
    assert all(z["after_observer_cap"]>=100 for z in new)


def test_wrong_previous_allocation_is_not_silently_replaced():
    p=pool()
    with pytest.raises(ValueError,match="unknown P100"):
        m.allocate(p,pool(1,start=1000),[{"prospective_rank":1,**p[0]}])
    with pytest.raises(ValueError,match="outside 3230 parent"):
        m.allocate(p,[p[0]],[{"prospective_rank":1,**pool(1,start=1000)[0]}])
    wrong=p[1].copy()
    wrong["species"]="Genus1 incorrect"
    with pytest.raises(ValueError,match="identity name has drifted"):
        m.allocate(p,[p[0]],[{"prospective_rank":1,**wrong}])


def test_duplicate_species_or_taxon_identity_blocked():
    p=pool()
    with pytest.raises(ValueError,match="Duplicate inat_taxon_id"):
        m.check_unique(p+[p[0]],"unit")
    x=[dict(y) for y in p]
    x[1]["species"]=x[0]["species"]
    with pytest.raises(ValueError,match="Duplicate species"):
        m.check_unique(x,"unit")


def test_frozen_third_rank_and_taxon_are_verified(tmp_path):
    f=tmp_path/"third.tsv"
    rows=pool(3)
    manifest(f,rows)
    got=m.read_third_manifest(f)
    assert len(got)==3
    assert got[0]["prospective_rank"]==1
    with f.open("w",newline="",encoding="utf-8") as o:
        o.write("prospective_rank\tinat_taxon_id\tspecies\n")
        o.write("1\t1\tGenus1 species\n")
        o.write("1\t2\tGenus2 species\n")
    with pytest.raises(ValueError,match="Third-cohort ranks"):
        m.read_third_manifest(f)


def test_canonical_candidate_csv_deterministic(tmp_path):
    p=tmp_path/"pool.csv"
    rows=pool(6)
    m.canonical_pool(p,rows)
    assert p.read_text().startswith("inat_taxon_id,species,after_observer_cap\n")
    assert len(p.read_text().strip().splitlines())==7
    s=tmp_path/"selected.tsv"
    _,selection,_=m.allocate(pool(6),[rows[0]],[{"prospective_rank":1,**rows[1]}])
    m.write_manifest(s,selection)
    assert s.read_text().splitlines()[0]=="\t".join(m.SAMPLED_COLUMNS)
    assert all(z["species"] not in (rows[0]["species"],rows[1]["species"])
               for z in selection)
