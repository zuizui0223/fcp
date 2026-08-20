"""Frozen analysis constants for the Journal of Biogeography paper."""

METRICS = [
    "temperature_breadth",
    "moisture_breadth",
    "climatic_heterogeneity",
    "pca_dispersion",
    "pca_hull_area",
]

EXPECTED_COUNTS = {
    "species": 34,
    "families": 25,
    "within_population": 20,
    "among_population": 14,
}

MODEL_FORMULA = "among ~ metric_z + effort_z"
DEFAULT_PERMUTATIONS = 9999
DEFAULT_SEED = 20260719
