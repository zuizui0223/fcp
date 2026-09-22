from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd

from fcp_pipeline.random_photo_h9_pool import (
    freeze_h9_metadata,
    geographic_maximin,
    observer_cap,
)

ROOT = Path(__file__).resolve().parents[1]


def _load_script():
    path = ROOT / "scripts" / "acquisition" / "run_fcp_v2_fresh_metadata_20260923.py"
    spec = importlib.util.spec_from_file_location("fcp_v2_fresh_metadata", path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


def _obs(obs_id: int, photo_id: int, taxon_id: int, lat: float, lon: float, observer: int = 1):
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


class FakeClient:
    def __init__(self):
        self.calls = []

    def observations(self, params):
        self.calls.append(dict(params))
        tid = int(params["taxon_id"])
        if tid == 2:
            raise RuntimeError("network")
        return {
            "results": [
                _obs(i, 1000 + i, tid, -20 + i, -100 + 4 * i, observer=i)
                for i in range(1, 8)
            ]
        }


def test_metadata_freeze_is_one_call_per_species_and_excludes_prior_ids():
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
    assert 1 not in set(frozen.observations["observation_id"].astype(int))
    assert set(frozen.observations["inat_taxon_id"].astype(int)) == {1}


def test_observer_cap_and_maximin_are_deterministic():
    rows = pd.DataFrame(
        [
            {
                "observation_id": i,
                "photo_id": 100 + i,
                "observer_id": "a" if i < 4 else str(i),
                "latitude": float(i),
                "longitude": float(i * 10),
            }
            for i in range(8)
        ]
    )
    capped = observer_cap(rows, 2)
    assert (capped["observer_id"] == "a").sum() == 2
    a = geographic_maximin(capped, 4)
    b = geographic_maximin(capped.sample(frac=1, random_state=1), 4)
    assert a["observation_id"].tolist() == b["observation_id"].tolist()


def test_queue_reader_rejects_hash_or_rank_drift(tmp_path):
    module = _load_script()
    q = pd.DataFrame(
        {
            "queue_rank": range(1, 301),
            "panel": ["P"] * 300,
            "inat_taxon_id": range(1, 301),
            "species": [f"S{i} sp" for i in range(1, 301)],
            "after_observer_cap": [100] * 300,
            "selection_hash": [f"{i:064x}" for i in range(300)],
        }
    )
    path = tmp_path / "queue.csv"
    q.to_csv(path, index=False)
    # The real frozen queue hash is intentionally required, so a synthetic queue
    # must fail rather than being silently accepted.
    try:
        module.read_queue(path, "P")
    except RuntimeError as exc:
        assert "queue SHA drift" in str(exc)
    else:
        raise AssertionError("synthetic queue unexpectedly passed frozen SHA guard")
