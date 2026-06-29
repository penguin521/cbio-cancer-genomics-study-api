"""
Standalone inspection script for viewing how cBioPortal data
moves through the API.

Run with:

    python inspect_cbioportal_data.py

This script does not modify api.py.
It imports and uses functions that already exist there.
"""

from pprint import pprint

from fastapi import HTTPException

from api import (
    cbioportal_get,
    clean_clinical_value,
    clean_clinical_attribute,
    clean_molecular_profile,
    fetch_clinical_values,
    clinical_long_to_wide,
    build_clinical_genomic_analysis_table
)


# ---------------------------------------------------------
# Change these values to inspect a different study or genes
# ---------------------------------------------------------

STUDY_ID = "acbc_mskcc_2015"

GENES = [
    "PALB2",   
    "CHEK2",
    "CDH1",
    "PTEN",
    "TP53"
]


def print_section(title: str):
    print("\n")
    print("=" * 70)
    print(title)
    print("=" * 70)


def inspect_raw_clinical_data():
    """
    Show the original cBioPortal field names and then
    show how clean_clinical_value() renames them.
    """
    print_section("1. RAW SAMPLE-LEVEL CLINICAL DATA")

    raw_records = cbioportal_get(
        f"/api/studies/{STUDY_ID}/clinical-data",
        params={
            "clinicalDataType": "SAMPLE",
            "projection": "SUMMARY",
            "pageSize": 3,
            "pageNumber": 0
        }
    )

    if not raw_records:
        print("No sample-level clinical records were returned.")
        return

    first_raw_record = raw_records[0]

    print("\nOriginal cBioPortal keys:")
    pprint(list(first_raw_record.keys()))

    print("\nOriginal cBioPortal record:")
    pprint(first_raw_record)

    cleaned_record = clean_clinical_value(
        first_raw_record
    )

    print("\nSame record after clean_clinical_value():")
    pprint(cleaned_record)

    print("\nImportant field-name conversion:")

    field_mapping = {
        "studyId": "study_id",
        "sampleId": "sample_id",
        "patientId": "patient_id",
        "clinicalAttributeId": (
            "clinical_attribute_id"
        ),
        "value": "value"
    }

    for raw_name, cleaned_name in field_mapping.items():
        print(f"  {raw_name}  →  {cleaned_name}")


def inspect_clinical_attributes():
    """
    Show clinical-attribute definitions supplied by cBioPortal.
    """
    print_section("2. CLINICAL ATTRIBUTE DEFINITIONS")

    raw_attributes = cbioportal_get(
        f"/api/studies/{STUDY_ID}/clinical-attributes",
        params={
            "projection": "DETAILED"
        }
    )

    print(
        f"\nNumber of clinical attributes returned: "
        f"{len(raw_attributes)}"
    )

    for raw_attribute in raw_attributes[:5]:
        print("\nRaw clinical attribute:")
        pprint(raw_attribute)

        print("\nCleaned clinical attribute:")
        pprint(
            clean_clinical_attribute(raw_attribute)
        )

        print("-" * 50)


def inspect_molecular_profiles():
    """
    Show molecular-profile fields such as datatype and
    molecularAlterationType.
    """
    print_section("3. MOLECULAR PROFILE DEFINITIONS")

    raw_profiles = cbioportal_get(
        f"/api/studies/{STUDY_ID}/molecular-profiles",
        params={
            "projection": "DETAILED"
        }
    )

    print(
        f"\nNumber of molecular profiles returned: "
        f"{len(raw_profiles)}"
    )

    for raw_profile in raw_profiles:
        print("\nRaw molecular profile:")
        pprint(raw_profile)

        print("\nCleaned molecular profile:")
        pprint(
            clean_molecular_profile(raw_profile)
        )

        print("-" * 50)


def inspect_clinical_transformation():
    """
    Reproduce the steps that create clinical_sample_count
    and clinical_patient_count.
    """
    print_section("4. HOW CLINICAL COUNTS ARE CREATED")

    # Retrieve every sample-level clinical value.
    sample_clinical_values = fetch_clinical_values(
        study_id=STUDY_ID,
        clinical_data_type="SAMPLE",
        fetch_all=True
    )

    print(
        "\nNumber of long-form sample clinical records:"
    )
    print(len(sample_clinical_values))

    print("\nExample long-form record:")
    if sample_clinical_values:
        pprint(sample_clinical_values[0])

    # Convert long-form values into one row per sample.
    sample_clinical_rows = clinical_long_to_wide(
        records=sample_clinical_values,
        clinical_data_type="SAMPLE"
    )

    print("\nNumber of wide sample clinical rows:")
    print(len(sample_clinical_rows))

    print("\nExample wide sample row:")
    if sample_clinical_rows:
        pprint(sample_clinical_rows[0])

    # Build the same lookup dictionary used in api.py.
    sample_clinical_by_sample = {
        row.get("sample_id"): row
        for row in sample_clinical_rows
        if row.get("sample_id")
    }

    clinical_sample_count = len(
        sample_clinical_by_sample
    )

    print("\nHow clinical_sample_count is calculated:")
    print(
        "clinical_sample_count = "
        "len(sample_clinical_by_sample)"
    )
    print(
        f"clinical_sample_count = "
        f"{clinical_sample_count}"
    )

    # Repeat the same process for patient-level records.
    patient_clinical_values = fetch_clinical_values(
        study_id=STUDY_ID,
        clinical_data_type="PATIENT",
        fetch_all=True
    )

    patient_clinical_rows = clinical_long_to_wide(
        records=patient_clinical_values,
        clinical_data_type="PATIENT"
    )

    patient_clinical_by_patient = {
        row.get("patient_id"): row
        for row in patient_clinical_rows
        if row.get("patient_id")
    }

    clinical_patient_count = len(
        patient_clinical_by_patient
    )

    print("\nHow clinical_patient_count is calculated:")
    print(
        "clinical_patient_count = "
        "len(patient_clinical_by_patient)"
    )
    print(
        f"clinical_patient_count = "
        f"{clinical_patient_count}"
    )


def inspect_final_analysis_table():
    """
    Build the full analysis table and show where the
    final output values come from.
    """
    print_section("5. FINAL CLINICAL-GENOMIC ANALYSIS TABLE")

    analysis_data = (
        build_clinical_genomic_analysis_table(
            study_id=STUDY_ID,
            gene_symbols=GENES
        )
    )

    print("\nDerived analysis metadata:")

    summary_fields = {
        "molecular_profile_id": analysis_data[
            "molecular_profile_id"
        ],
        "sample_list_id": analysis_data[
            "sample_list_id"
        ],
        "requested_genes": analysis_data[
            "requested_genes"
        ],
        "successful_genes": analysis_data[
            "successful_genes"
        ],
        "profiled_sample_count": analysis_data[
            "profiled_sample_count"
        ],
        "clinical_sample_count": analysis_data[
            "clinical_sample_count"
        ],
        "clinical_patient_count": analysis_data[
            "clinical_patient_count"
        ],
        "analysis_row_count": analysis_data[
            "analysis_row_count"
        ]
    }

    pprint(summary_fields)

    rows = analysis_data["analysis_table"]

    print("\nNumber of final analysis rows:")
    print(len(rows))

    if rows:
        print("\nColumn names in the first analysis row:")
        pprint(list(rows[0].keys()))

        print("\nFirst complete analysis row:")
        pprint(rows[0])

    print("\nMutation retrieval errors:")
    pprint(analysis_data["errors"])


def main():
    try:
        inspect_raw_clinical_data()
        inspect_clinical_attributes()
        inspect_molecular_profiles()
        inspect_clinical_transformation()
        inspect_final_analysis_table()

    except HTTPException as error:
        print("\nThe API returned an error:")
        print(f"Status code: {error.status_code}")
        pprint(error.detail)

    except Exception as error:
        print("\nUnexpected Python error:")
        print(type(error).__name__)
        print(error)


if __name__ == "__main__":
    main()
