"""Review-queue contract: truly separate reviewer-facing IDs from source colours."""
import importlib.util
from pathlib import Path
import pandas as pd
import numpy as np

path=Path(__file__).resolve().parents[1]/"scripts/analysis/build_polymorphism_blinded_photo_audit_queue_20261010.py"
spec=importlib.util.spec_from_file_location("fcp_blindqueue",path)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_haversine_nearest_distance_and_hash_determinism():
    points=pd.DataFrame({"latitude":[0.,0.],"longitude":[0.,2.]})
    z=m.get_haversine(0.,0.,points)
    assert z[0]==0 and 200<z[1]<250
    assert m.token("x",123)==m.token("x",123)
    assert m.token("x",123)!=m.token("x",124)


def test_review_queue_contains_all_target_photos_and_blinds_labels(tmp_path,monkeypatch):
    counts={"discovery":27,"validation":23,"third":20}
    source={}
    swaps=[]
    for cohort,count in counts.items():
        photo=[str(i+100000) for i in range(300)]
        labels=["white"]*150+["red_pink"]*150
        source[cohort]=pd.DataFrame({
          "photo_id":photo, "inat_taxon_id":[7]*300,
          "species":["Plant exemplar"]*300, "morph":labels,
          "latitude":[0.]*300,"longitude":[i/1000 for i in range(300)]})
        for i in range(count):
            swaps.append({"cohort":cohort,"photo_id_i":photo[i],
                          "photo_id_j":photo[150+i]})
    table=tmp_path/"swaps.csv"
    pd.DataFrame(swaps).to_csv(table,index=False)
    monkeypatch.setattr(m.attack,"load",lambda path,cohort:source[cohort])
    monkeypatch.setattr(m.attack,"build_species",lambda frame:[type("S",(),{"taxon":7})()])
    out=tmp_path/"queue"
    result=m.run({c:tmp_path/f"{c}.csv" for c in counts},table,out)
    assert result["n_high_leverage_photo_cases"]==140
    assert result["n_species_colour_matched_controls"]==140
    assert result["n_random_source_controls"]==140
    assert result["n_total_unique_review_photo_ids"]==420
    blind=pd.read_csv(out/"BLINDED_reannotation_photo_queue.csv",dtype=str)
    key=pd.read_csv(out/"UNBLINDING_KEY_do_not_show_reviewers.csv",dtype=str)
    assert len(blind)==420 and blind.audit_case_id.is_unique
    assert len(key)==420 and key[["cohort","photo_id"]].drop_duplicates().shape[0]==420
    assert not set(["kind","morph","original_algorithm_colour_label","cohort","species"]).intersection(blind.columns)
    assert {"original_algorithm_colour_label","kind","reference_photo_id"}.issubset(set(key.columns))
