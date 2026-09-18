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
from .reliability import ReliabilityResult, observer_disjoint_reliability
from .spatial import (
    AssociationResult,
    distribution_spatial_association,
    great_circle_pairwise_km,
    jensen_shannon_pairwise,
    matched_difference_spatial_permutation_null,
    matched_difference_spatial_rho,
    partial_distribution_spatial_association,
    partial_rank_correlation,
    spatial_permutation_null,
    spatial_rho,
)

__all__ = [
    "AssociationResult",
    "ReliabilityResult",
    "TwoModeResult",
    "alignment_statistic",
    "categorical_diversity",
    "distribution_spatial_association",
    "gini_simpson",
    "great_circle_pairwise_km",
    "hellinger_rows",
    "jensen_shannon_pairwise",
    "matched_difference_spatial_permutation_null",
    "matched_difference_spatial_rho",
    "observer_disjoint_reliability",
    "one_vs_rest_contrast",
    "partial_distribution_spatial_association",
    "partial_rank_correlation",
    "permute_rows_within_strata",
    "spatial_permutation_null",
    "spatial_rho",
    "two_mode_axis",
]
