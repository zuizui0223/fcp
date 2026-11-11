"""Synthetic 50-km spatial estimator regression checks."""
import numpy as np
import pandas as pd
import pytest
from fcp_pipeline.u100_frozen_spatial_export import spatial_row, frozen_check, observed_spatial_table


def group(labels, lat, lon):
    return pd.DataFrame({"inat_taxon_id": 17, "cohort": "prospective", "morph": labels,
                         "usable_colour": True, "latitude": lat, "longitude": lon})


def test_spatial_partitioning_uses_unique_unordered_photo_pairs():
    # At 0 and 20 degrees, within-site pairs have identical colours.
    g = group(["white"] * 20 + ["red_pink"] * 20,
              [0.] * 20 + [20.] * 20, [130.] * 40)
    r = spatial_row(g)
    assert r["D_pair"] == pytest.approx(400 / 780)
    assert r["local_pairs_50km"] == 380
    assert r["depletion_50km"] == pytest.approx(400 / 780)


def test_fully_white_photos_not_confused_with_unmeasured():
    g = group(["white"] * 40, [0.] * 40, [130.] * 40)
    result = observed_spatial_table(g)
    assert result.iloc[0]["D_pair"] == 0
    assert result.iloc[0]["depletion_50km"] == 0


def test_unsupported_scale_is_nan_and_untrusted_reference_fails():
    # Forty isolated latitudes: small number of local pairs.
    g = group(["white"] * 20 + ["blue_purple"] * 20,
              list(np.linspace(-80, 80, 40)), [140.] * 40)
    x = spatial_row(g)
    assert x["local_pairs_50km"] < 30
    assert np.isnan(x["depletion_50km"])
    with pytest.raises(ValueError, match="wrong frozen primary radius"):
        frozen_check(pd.DataFrame([x]), {"primary_radius_km": 100})
