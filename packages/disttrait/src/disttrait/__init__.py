"""Validated species-level distributional trait inference."""

from .diversity import categorical_diversity, gini_simpson
from .geometry import (
    TwoModeResult,
    alignment_statistic,
    hellinger_rows,
    one_vs_rest_contrast,
    permute_rows_within_strata,
    two_mode_axis,
)

__all__ = [
    "TwoModeResult",
    "alignment_statistic",
    "categorical_diversity",
    "gini_simpson",
    "hellinger_rows",
    "one_vs_rest_contrast",
    "permute_rows_within_strata",
    "two_mode_axis",
]
