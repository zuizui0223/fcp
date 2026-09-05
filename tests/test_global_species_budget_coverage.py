from __future__ import annotations

import numpy as np
import pandas as pd

from fcp_pipeline.global_species_budget_coverage import (
    cell_count_correlation,
    gini_nonnegative,
    occupied_cell_retention,
    summarize_cell_coverage,
    unique_taxon_cell_links,
)


def test_unique_links_deduplicate_across_discovery_generations():
    v1 = pd.DataFrame({"inat_taxon_id": [1, 1, 2], "cell_id": [0, 0, 1]})
    v2 = pd.DataFrame({"inat_taxon_id": [1, 3], "cell_id": [0, 2]})
    links = unique_taxon_cell_links(v1, v2)
    assert len(links) == 3


def test_coverage_retention_and_counts_are_exact():
    links = pd.DataFrame({
        "inat_taxon_id": [1, 1, 2, 3, 4],
        "cell_id": [0, 1, 1, 2, 3],
    })
    full = summarize_cell_coverage(links, [1, 2, 3, 4], n_cells=5)
    subset = summarize_cell_coverage(links, [1, 2], n_cells=5)
    assert full.occupied_cells == 4
    assert subset.occupied_cells == 2
    assert occupied_cell_retention(full, subset) == 0.5
    assert full.total_taxon_cell_links == 5
    assert subset.total_taxon_cell_links == 3


def test_equal_cell_counts_have_zero_gini_and_perfect_correlation():
    assert np.isclose(gini_nonnegative([2, 2, 2, 2]), 0.0)
    links = pd.DataFrame({
        "inat_taxon_id": [1, 2, 3, 4],
        "cell_id": [0, 1, 2, 3],
    })
    a = summarize_cell_coverage(links, [1, 2, 3, 4], n_cells=4)
    b = summarize_cell_coverage(links, [1, 2, 3, 4], n_cells=4)
    # Constant counts make Pearson undefined by design.
    assert np.isnan(cell_count_correlation(a, b))


def test_nonconstant_cell_pattern_correlation_is_one_for_scaled_pattern():
    links = pd.DataFrame({
        "inat_taxon_id": [1, 2, 3, 4, 5, 6],
        "cell_id": [0, 0, 0, 1, 1, 2],
    })
    a = summarize_cell_coverage(links, [1, 2, 3, 4, 5, 6], n_cells=4)
    # Create a second summary with the same relative pattern directly.
    b = type(a)(
        taxon_count=12,
        occupied_cells=3,
        total_taxon_cell_links=12,
        cell_species_gini=gini_nonnegative([6, 4, 2, 0]),
        cell_species_counts=np.array([6, 4, 2, 0]),
    )
    assert np.isclose(cell_count_correlation(a, b), 1.0)
