import pandas as pd

from fcp_pipeline.random_photo_h9_pool import (
    freeze_h9_metadata,
    geographic_maximin,
    h9_query_for_species,
    observer_cap,
    parse_h9_observation,
)


def obs(obs_id, photo_id, taxon_id, lat, lon, observer=1):
    return {
        "id": obs_id,
        "quality_grade": "research",
        "positional_accuracy": 100,
        "obscured": False,
        "geoprivacy": None,
        "geojson": {"coordinates": [lon, lat]},
        "taxon": {"id": taxon_id, "rank": "species", "name": f"Species {taxon_id}"},
        "user": {"id": observer, "login": f"u{observer}"},
        "photos": [{"id": photo_id, "license_code": "cc-by", "url": "https://x/square.jpg"}],
    }


def test_query_is_species_specific_and_random_page():
    q = h9_query_for_species(123, per_page=200)
    assert q["taxon_id"] == 123
    assert q["rank"] == "species"
    assert q["order_by"] == "random"
    assert q["per_page"] == 200 and q["page"] == 1


def test_parser_requires_exact_taxon():
    allowed = frozenset({"cc-by"})
    assert parse_h9_observation(obs(1, 11, 123, 10, 20), expected_taxon_id=123, maximum_positional_accuracy_m=5000, allowed_photo_licenses=allowed)
    assert parse_h9_observation(obs(1, 11, 124, 10, 20), expected_taxon_id=123, maximum_positional_accuracy_m=5000, allowed_photo_licenses=allowed) is None


def test_observer_cap_and_maximin_are_deterministic():
    rows = pd.DataFrame([
        {"observation_id": i, "photo_id": 100+i, "observer_id": "a" if i < 4 else str(i), "latitude": float(i), "longitude": float(i*10)}
        for i in range(8)
    ])
    capped = observer_cap(rows, 2)
    assert (capped["observer_id"] == "a").sum() == 2
    a = geographic_maximin(capped, 4)
    b = geographic_maximin(capped.sample(frac=1, random_state=1), 4)
    assert a["observation_id"].tolist() == b["observation_id"].tolist()


class FakeClient:
    def __init__(self):
        self.calls = []
    def observations(self, params):
        self.calls.append(dict(params))
        tid = int(params["taxon_id"])
        if tid == 2:
            raise RuntimeError("network")
        results = [obs(i, 1000+i, tid, -20+i, -100+4*i, observer=i) for i in range(1, 8)]
        results.append(obs(90, 1090, tid+999, 0, 0, observer=90))
        return {"results": results}


def test_freeze_uses_one_call_per_species_excludes_prior_and_records_error():
    frame = pd.DataFrame({"species": ["A a", "B b"], "inat_taxon_id": [1, 2]})
    client = FakeClient()
    frozen = freeze_h9_metadata(
        client=client,
        species_frame=frame,
        exclusion_observation_ids={1},
        exclusion_photo_ids=set(),
        per_page=200,
        observer_cap_n=2,
        fixed_raw_photos=3,
    )
    assert len(client.calls) == 2
    assert frozen.manifest["request_errors"] == 1
    assert frozen.species_audit.loc[frozen.species_audit["inat_taxon_id"] == 1, "retained"].iloc[0] == 3
    assert frozen.species_audit.loc[frozen.species_audit["inat_taxon_id"] == 2, "retained"].iloc[0] == 0
    assert 1 not in set(frozen.observations["observation_id"].astype(int))
    assert set(frozen.observations["inat_taxon_id"].astype(int)) == {1}
