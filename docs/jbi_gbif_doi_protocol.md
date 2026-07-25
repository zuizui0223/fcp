# GBIF DOI protocol for the Journal of Biogeography submission

## Why two GBIF citation records are needed

The paginated sensitivity analysis did not use an unmodified GBIF.org download. It retrieved records synchronously, capped each species at 3,000 records, retained specified basis-of-record classes, excluded geospatial issues and coordinate uncertainty above 20 km, retained missing uncertainty values and removed duplicate coordinate pairs rounded to 0.001° within each species.

A regular GBIF occurrence-download DOI can document the broad predicate-based retrieval, but the standard download service cannot reproduce the local per-species cap or coordinate deduplication. The exact analytical subset therefore also requires a GBIF **Derived Dataset** record after the frozen file is deposited at a permanent public URL.

## Prepared and verified bundle

Directory: `docs/supporting/jbi_gbif_doi_bundle/`

| File | Purpose |
|---|---|
| `jbi_gbif_exact_occurrence_subset.csv.gz` | Exact 58,455-row occurrence subset used before climate extraction. |
| `jbi_gbif_parent_dataset_counts.csv` | Contribution counts for all 389 parent GBIF datasets. |
| `jbi_gbif_exact_species_counts.csv` | Retained record counts for the 34 audited species. |
| `jbi_gbif_broad_download_request.json` | Standard GBIF occurrence-download request approximating the broad pre-filtered retrieval. |
| `jbi_gbif_derived_dataset_metadata_template.json` | Title, description, archive hash and placeholders for Derived Dataset registration. |
| `jbi_gbif_doi_bundle_manifest.json` | Source artifact, row counts, file hashes and citation boundary. |

Verified integrity:

- source workflow run: `30065308023`;
- source artifact: `8586269620`;
- source artifact digest: `sha256:4247fe3dae52f7723a25df6cd9956fa8b2f220a3fc101e088a3f0cccd768447a`;
- source CSV: 58,455 rows, 34 species, 389 parent datasets and 58,455 unique GBIF occurrence keys;
- source CSV SHA-256: `b0614a729acde5a1daab599d52c39ac4018583e1f73d83e5304ac0afa6f6e7ad`;
- deterministic gzip SHA-256: `f25ae0cf2c84c45ae461a932d6c6063edda64591913a2495e4a3da82d573f094`;
- DOI-bundle workflow run: `30081674018`;
- DOI-bundle artifact: `8592043987`;
- DOI-bundle artifact digest: `sha256:61f6d8549ce59a2d52e438af7430af7556ffdf0cd13210579b0b53e3854925a1`.

## Step 1: request the broad GBIF download DOI

A registered GBIF account is required. Credentials must never be committed.

```bash
export GBIF_PASSWORD='your-password'
python submit_jbi_gbif_download.py \
  --username YOUR_GBIF_USERNAME \
  --email YOUR_GBIF_ACCOUNT_EMAIL \
  --out jbi_gbif_download_key.txt
```

The submission script reads the prepared JSON request, replaces only the notification email and sends the request to the GBIF occurrence-download API using HTTP Basic authentication. The returned download key is stored locally. Wait until the download status is `SUCCEEDED`, then record the DOI shown by GBIF.

The prepared predicate includes:

- 34 numeric `TAXON_KEY` values from the audited GBIF resolution table;
- `HAS_COORDINATE = true`;
- `HAS_GEOSPATIAL_ISSUE = false`;
- `OCCURRENCE_STATUS = PRESENT`;
- the five retained basis-of-record classes; and
- coordinate uncertainty ≤20,000 m or a null uncertainty value.

No `checklistKey` is supplied because the original numeric keys were resolved through the legacy v1 GBIF Backbone species-match service. Changing taxonomy during DOI preparation would no longer reproduce the original query basis.

### Citation boundary for the broad DOI

The broad download may contain more records than the analysed subset because it does not apply the local 3,000-record-per-species cap or 0.001° coordinate deduplication. It may be cited as the source retrieval but must not be described as the exact 58,455-row analysis file.

## Step 2: create a permanent archive for the exact subset

Create a frozen release or archive containing at least:

- `jbi_gbif_exact_occurrence_subset.csv.gz`;
- `jbi_gbif_parent_dataset_counts.csv`;
- `jbi_gbif_exact_species_counts.csv`;
- `jbi_gbif_doi_bundle_manifest.json`;
- the analysis scripts and the manuscript-specific README; and
- the completed broad-download DOI, when available.

Deposit the release in a permanent public repository such as Zenodo. Verify the uploaded exact archive against SHA-256 `f25ae0cf2c84c45ae461a932d6c6063edda64591913a2495e4a3da82d573f094`. Copy the permanent landing-page URL into `jbi_gbif_derived_dataset_metadata_template.json`.

## Step 3: register the exact GBIF Derived Dataset

Sign in to the GBIF Derived Dataset tool and provide:

1. the title from `jbi_gbif_derived_dataset_metadata_template.json`;
2. the permanent public landing-page URL of the frozen exact subset;
3. `jbi_gbif_parent_dataset_counts.csv`, containing one `datasetKey,record_count` row for each of 389 contributing datasets;
4. the prepared description of the retrieval, cap, filters and deduplication; and
5. optionally, the broad GBIF download DOI from Step 1 as the original download DOI.

After submission, copy the Derived Dataset DOI and suggested citation into the release records and manuscript. This Derived Dataset DOI is the citation that most closely represents the exact analytical subset.

## Step 4: update the manuscript and submission files

Replace the remaining GBIF DOI placeholder in:

- `docs/jbi_manuscript_editorial_revision.md`;
- `docs/jbi_title_page_template.md`;
- `docs/jbi_cover_letter_template.md`;
- `docs/jbi_submission_completion_checklist.md`; and
- the final archived release README.

Recommended Data Accessibility wording after both records exist:

> The exact code and analysis files are archived at [PERMANENT RELEASE DOI]. GBIF-mediated occurrence data are represented by the broad occurrence download [GBIF DOWNLOAD DOI] and by the exact locally processed subset registered as a GBIF Derived Dataset [GBIF DERIVED DATASET DOI].

Use the exact citation strings supplied by GBIF rather than inventing author names, years or dataset titles.

## Items that still require account-holder action

- a GBIF username, password and verified account email;
- authenticated submission of the broad request;
- a permanent public archive URL for the exact subset;
- authenticated registration of the Derived Dataset; and
- insertion of the issued DOI strings into the frozen submission files.
