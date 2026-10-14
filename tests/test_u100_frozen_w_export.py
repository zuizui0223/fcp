"""Synthetic-only tests for frozen W_i export; never count as biological evidence."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from fcp_pipeline.u100_frozen_w_export import TARGET, export_rows, verified_export


def make_vectors():
    rows = []
    for name, component in (("Species one", "white"), ("Species two", "pink")):
        row = {"species": name, "n_classifiable": 50, "H1_coarse_second_fraction": 0.25,
               "continuous_minor_fraction": 0.2}
        for colour in ("white", "yellow", "orange", "red", "pink", "magenta", "purple", "blue", "bronze"):
            row[f"delta_{colour}"] = 1.0 if colour == component else 0.0
        row["delta_norm"] = 1.0
        rows.append(row)
    return pd.DataFrame(rows)


def make_frame():
    return pd.DataFrame([{"inat_taxon_id": 111, "species": "Species one"},
                         {"inat_taxon_id": 112, "species": "Species two"}])


def test_fixed_unit_white_axis_without_refitting_modes():
    vectors = make_vectors()
    rows, mean_w = export_rows(vectors, make_frame())
    assert rows["w_alignment_sq"].tolist() == pytest.approx([TARGET[0] ** 2, TARGET[4] ** 2])
    assert mean_w == pytest.approx((TARGET[0] ** 2 + TARGET[4] ** 2) / 2)
    assert set(rows["method"]) == {"frozen_study_specific_H2"}
    assert rows["inat_taxon_id"].tolist() == [111, 112]


def test_missing_species_and_invalid_vectors_fail():
    v = make_vectors()
    with pytest.raises(ValueError, match="missing in U100"):
        export_rows(v, make_frame().iloc[:1])
    v.loc[0, "delta_norm"] = 2.0
    with pytest.raises(ValueError, match="vector/norm"):
        export_rows(v, make_frame())
    v = make_vectors()
    v.loc[1, "continuous_minor_fraction"] = .09
    with pytest.raises(ValueError, match="qualification"):
        export_rows(v, make_frame())


def test_only_original_pinned_vector_file_is_accepted(tmp_path: Path):
    csv = tmp_path / "v.csv"
    u100 = tmp_path / "roster.csv"
    make_vectors().to_csv(csv, index=False)
    make_frame().to_csv(u100, index=False)
    with pytest.raises(ValueError, match="Git blob identity mismatch"):
        verified_export(u100, csv, tmp_path / "w.csv")
