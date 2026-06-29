# Cancer Genomics Study Explorer API

## Overview

The Cancer Genomics Study Explorer API is a FastAPI application that retrieves public cancer genomics data from cBioPortal and converts it into structured, analysis-ready outputs.

The application is designed to simplify common research workflows that would otherwise require multiple API requests, manual data reshaping, and repeated joins between study, sample, patient, clinical, and mutation datasets.

## Project Objectives

This project was developed to:

* simplify access to cBioPortal study data,
* standardize clinical and mutation records,
* integrate sample-level and patient-level information,
* generate sample-by-gene mutation matrices,
* construct reproducible clinical-genomic cohorts,
* assess data completeness,
* and export results in CSV format.

The application does not generate new biological data. Its value lies in organizing and combining existing cBioPortal data into consistent, reusable tables.

## Core Functionality

### Study discovery

Users can search for studies by keyword and retrieve concise study summaries, including study identifiers, cancer types, sample counts, and available data categories.

### Clinical data processing

Clinical data can be retrieved in:

* long format, with one clinical observation per row,
* wide format, with one sample or patient per row.

Clinical data can also be exported directly as CSV.

### Mutation analysis

The API supports:

* mutation retrieval for individual genes,
* single-gene mutation summaries,
* multi-gene mutation panel summaries,
* sample-by-gene mutation matrices,
* mutation counts,
* mutation frequencies,
* mutation types,
* and protein-change annotations.

### Clinical-genomic integration

The analysis-table endpoint combines:

```text
sample metadata
+ sample-level clinical data
+ patient-level clinical data
+ mutation status
+ mutation counts
+ protein changes
```

Each output row represents one profiled sample.

### Cohort construction

Users can define filtered cohorts using:

* mutated or non-mutated status for a selected gene,
* exact clinical values,
* minimum numeric thresholds,
* maximum numeric thresholds,
* or combined mutation and clinical criteria.

Filtered cohorts can be exported as CSV files.

### Data-quality assessment

The quality-summary endpoint reports:

* row and column counts,
* missing-value counts,
* missing percentages,
* unique-value counts,
* sparse columns,
* numeric ranges,
* means,
* medians,
* and mutation frequencies.

### Reliability

The application includes:

* automatic pagination for large cBioPortal responses,
* configurable connection and read timeouts,
* retry handling for temporary HTTP failures,
* and automated tests for key filtering and quality-summary logic.

## Workflow

```text
Search for a study
→ inspect available metadata
→ retrieve clinical and mutation data
→ reshape and integrate records
→ define a cohort
→ assess data quality
→ export an analysis-ready table
```

## Technology Stack

* Python
* FastAPI
* Uvicorn
* Requests
* Pytest
* cBioPortal public API

## Project Structure

```text
cancer-genomics-api/
│
├── api.py
├── test_api.py
├── README.md
└── requirements.txt
```

## Installation

Create and activate a virtual environment:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
```

Install dependencies:

```powershell
python -m pip install fastapi uvicorn requests pytest
```

Alternatively:

```powershell
python -m pip install -r requirements.txt
```

## Running the Application

Start the development server:

```powershell
python -m uvicorn api:app --reload
```

Open the interactive API documentation:

```text
http://localhost:8000/docs
```

Health check:

```text
http://localhost:8000/alive
```

Expected response:

```json
{
  "status": "running"
}
```

## Example Usage

The following examples use the TCGA PanCancer Atlas melanoma study:

```text
skcm_tcga_pan_can_atlas_2018
```

Search for melanoma studies:

```text
http://localhost:8000/studies/search?keyword=melanoma
```

View a study overview:

```text
http://localhost:8000/studies/skcm_tcga_pan_can_atlas_2018/overview
```

Build a multi-gene mutation summary:

```text
http://localhost:8000/studies/skcm_tcga_pan_can_atlas_2018/mutation-panel/summary?genes=BRAF,NRAS,NF1,KIT
```

Build a combined clinical-genomic analysis table:

```text
http://localhost:8000/studies/skcm_tcga_pan_can_atlas_2018/analysis-table?genes=BRAF,NRAS,NF1,KIT
```

Inspect data completeness:

```text
http://localhost:8000/studies/skcm_tcga_pan_can_atlas_2018/analysis-table/quality?genes=BRAF,NRAS
```

Filter for BRAF-mutated samples:

```text
http://localhost:8000/studies/skcm_tcga_pan_can_atlas_2018/cohort?genes=BRAF&mutation_gene=BRAF&mutation_status=mutated
```

Export the filtered cohort:

```text
http://localhost:8000/studies/skcm_tcga_pan_can_atlas_2018/cohort/export.csv?genes=BRAF&mutation_gene=BRAF&mutation_status=mutated
```

## Selected Endpoints

| Endpoint                                     | Description                                 |
| -------------------------------------------- | ------------------------------------------- |
| `/studies/search`                            | Search studies by keyword                   |
| `/studies/{study_id}/overview`               | Summarize study metadata and available data |
| `/studies/{study_id}/clinical-data/long`     | Retrieve clinical data in long format       |
| `/studies/{study_id}/clinical-data/wide`     | Retrieve clinical data in wide format       |
| `/studies/{study_id}/mutation-panel/summary` | Summarize mutations across multiple genes   |
| `/studies/{study_id}/mutation-panel/matrix`  | Build a sample-by-gene mutation matrix      |
| `/studies/{study_id}/analysis-table`         | Combine clinical and mutation data          |
| `/studies/{study_id}/analysis-table/quality` | Assess missingness and column coverage      |
| `/studies/{study_id}/cohort`                 | Build a filtered clinical-genomic cohort    |
| `/studies/{study_id}/cohort/export.csv`      | Export a filtered cohort                    |

## Testing

Run the automated test suite:

```powershell
python -m pytest -q
```

The tests use synthetic data and do not require a live cBioPortal connection.

They verify:

* unfiltered cohort generation,
* mutated and non-mutated filtering,
* case-insensitive clinical filtering,
* numeric threshold filtering,
* invalid filter combinations,
* and missing-value calculations.

Check syntax separately with:

```powershell
python -m py_compile api.py
```

No output indicates that the file compiled successfully.

## Interpretation Notes

### Missing values

Missing clinical values should not be interpreted as zero, negative, normal, wild type, or otherwise absent.

They indicate that no usable value was available for that field in the selected study.

### Mutation status

A sample classified as `not_mutated` has no retrieved mutation record for the selected gene, molecular profile, and sample list.

This does not establish biological wild-type status under all possible assays or sequencing methods.

### Patient and sample attributes

Clinical attributes may exist at either the sample level or the patient level.

When the same attribute name appears at both levels, the patient-level value is prefixed with:

```text
patient_
```

This prevents patient-level data from overwriting sample-level data.

## Limitations

* All source data are obtained from cBioPortal.
* Available studies and clinical attributes vary in completeness and structure.
* Clinical terminology is not fully standardized across studies.
* The combined analysis table currently focuses on mutation data.
* The cohort builder supports one mutation-gene filter at a time.
* The cohort builder supports one clinical attribute at a time.
* Gene-expression and copy-number data are not yet integrated into the combined table.
* The application does not perform survival analysis, causal inference, or clinical interpretation.
* Results depend on the molecular profile and sample list selected for each study.

## Potential Future Development

Possible future extensions include:

* multiple simultaneous clinical filters,
* Boolean gene-filter logic,
* gene-expression integration,
* copy-number alteration integration,
* survival-analysis endpoints,
* cross-study harmonization,
* caching,
* Docker support,
* and a graphical user interface.

## Data Source

This application uses the public cBioPortal API.

Study contents, service availability, and upstream API behavior are managed externally by cBioPortal.

## Disclaimer

This software is intended for research, educational, and software-development purposes.

It is not intended for clinical diagnosis, treatment selection, or medical decision-making.

## Author

Steven Xie
