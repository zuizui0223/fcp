"""Validate the empirical seed outcome table using only the standard library."""
from __future__ import annotations

import csv
from pathlib import Path

ALLOWED_OUTCOMES = {
    "monomorphic",
    "geographic_mosaic",
    "mixed",
    "local_coexistence",
}
ALLOWED_TRISTATE = {"yes", "no", "not_coded"}
ALLOWED_CONFIDENCE = {"low", "medium", "high"}
REQUIRED_COLUMNS = {
    "species",
    "outcome_class",
    "within_population_coexistence",
    "among_population_differentiation",
    "evidence_scope",
    "source_year",
    "source_doi",
    "confidence",
    "notes",
}


def validate(path: Path = Path("data/seed_outcomes.csv")) -> list[str]:
    errors: list[str] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            return ["missing header"]
        missing = REQUIRED_COLUMNS - set(reader.fieldnames)
        if missing:
            errors.append(f"missing columns: {sorted(missing)}")
        seen_species: set[str] = set()
        for line_no, row in enumerate(reader, start=2):
            species = row.get("species", "").strip()
            if not species:
                errors.append(f"line {line_no}: empty species")
            elif species in seen_species:
                errors.append(f"line {line_no}: duplicate species {species}")
            seen_species.add(species)

            outcome = row.get("outcome_class", "")
            if outcome not in ALLOWED_OUTCOMES:
                errors.append(f"line {line_no}: invalid outcome_class {outcome!r}")
            for field in ("within_population_coexistence", "among_population_differentiation"):
                value = row.get(field, "")
                if value not in ALLOWED_TRISTATE:
                    errors.append(f"line {line_no}: invalid {field} {value!r}")
            confidence = row.get("confidence", "")
            if confidence not in ALLOWED_CONFIDENCE:
                errors.append(f"line {line_no}: invalid confidence {confidence!r}")
            doi = row.get("source_doi", "")
            if not doi.startswith("10."):
                errors.append(f"line {line_no}: malformed DOI {doi!r}")

            # Internal logical consistency for the outcome hierarchy.
            within = row.get("within_population_coexistence", "")
            among = row.get("among_population_differentiation", "")
            if outcome == "local_coexistence" and within != "yes":
                errors.append(f"line {line_no}: local_coexistence requires within_population_coexistence=yes")
            if outcome == "geographic_mosaic" and among != "yes":
                errors.append(f"line {line_no}: geographic_mosaic requires among_population_differentiation=yes")
            if outcome == "mixed" and not (within == "yes" and among == "yes"):
                errors.append(f"line {line_no}: mixed requires both spatial evidence flags=yes")
    return errors


if __name__ == "__main__":
    problems = validate()
    if problems:
        raise SystemExit("\n".join(problems))
    print("seed outcome table: OK")
