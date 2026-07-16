#!/usr/bin/env python3
"""Run the validated discovery pipeline with a broader search vocabulary.

The underlying mapping and context filters remain in global_literature_discovery.py.
This wrapper expands recall without duplicating or weakening those safeguards.
"""
from __future__ import annotations

import global_literature_discovery as discovery

EXPANDED_QUERIES = (
    # Core polymorphism terminology.
    '"flower color polymorphism"',
    '"flower colour polymorphism"',
    '"floral color polymorphism"',
    '"floral colour polymorphism"',
    '"petal color polymorphism"',
    '"petal colour polymorphism"',
    '"corolla color polymorphism"',
    '"corolla colour polymorphism"',
    # Morph and dimorphism terminology.
    '"flower color morph"',
    '"flower colour morph"',
    '"floral color morph"',
    '"floral colour morph"',
    '"flower color dimorphism"',
    '"flower colour dimorphism"',
    '"floral color dimorphism"',
    '"floral colour dimorphism"',
    '"petal color morph"',
    '"petal colour morph"',
    # Broader intraspecific variation, constrained by ecology/evolution terms.
    '"flower color variation" population',
    '"flower colour variation" population',
    '"floral color variation" population',
    '"floral colour variation" population',
    '"corolla color variation" population',
    '"corolla colour variation" population',
    '"flower color variation" selection',
    '"flower colour variation" selection',
    '"flower color variation" pollinator',
    '"flower colour variation" pollinator',
    # Mechanistic and geographic wording often used without “polymorphism”.
    '"flower morph frequency" color',
    '"flower morph frequency" colour',
    '"geographic variation" "flower color"',
    '"geographic variation" "flower colour"',
    '"spatial variation" "flower color"',
    '"spatial variation" "flower colour"',
    '"anthocyanin polymorphism" flower',
    '"pigment polymorphism" flower',
    '"white-flowered morph" population',
    '"white flower morph" population',
    '"colour morphs" pollination plant',
    '"color morphs" pollination plant',
)


def main() -> None:
    discovery.QUERIES = EXPANDED_QUERIES
    discovery.main()


if __name__ == "__main__":
    main()
