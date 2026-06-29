from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import Response
import requests
import csv
import io
from typing import Optional, Literal
from statistics import mean, median
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

app = FastAPI(
    title="Cancer Genomics Study Explorer API",
    description=(
        "A researcher-friendly API that connects to the public cBioPortal API "
        "and returns simplified cancer genomics study and sample metadata."
    ),
    version="1.1.0"
)

CBIOPORTAL_BASE_URL = "https://www.cbioportal.org"

CBIOPORTAL_CONNECT_TIMEOUT = 5
CBIOPORTAL_READ_TIMEOUT = 60
CBIOPORTAL_PAGE_SIZE = 1000

retry_policy = Retry(
    total=3,
    connect=3,
    read=3,
    status=3,
    backoff_factor=0.5,
    status_forcelist=[429,500,502,503,504],
    allowed_methods=frozenset(["GET"]),
    respect_retry_after_header=True,
    raise_on_status=False
)

CBIOPORTAL_SESSION = requests.Session()

CBIOPORTAL_SESSION.mount(
    "https://",
    HTTPAdapter(max_retries=retry_policy)
)

"""
load API in shell, by typing
    python -m uvicorn api:app --reload
        Local example:

        http://localhost:8000/docs

        http://localhost:8000/studies/clean

        http://localhost:8000/studies/search?keyword=breast     

        http://127.0.0.1:8000/studies/search?keyword=lung

        http://127.0.0.1:8000/studies/brca_tcga/summary

        http://127.0.0.1:8000/studies/brca_tcga/samples/clean

        http://127.0.0.1:8000/studies/brca_tcga/clinical-data/long?clinical_data_type=SAMPLE&page_size=100

        http://127.0.0.1:8000/studies/brca_tcga/clinical-data/long?clinical_data_type=SAMPLE&attribute_id=MUTATION_COUNT

        http://127.0.0.1:8000/studies/brca_tcga/clinical-data/export.csv?clinical_data_type=SAMPLE&format=wide
"""

def cbioportal_get(endpoint: str, params: Optional[dict] = None):
    """
    Send a GET request to the public cBioPortal API.
    """
    url = f"{CBIOPORTAL_BASE_URL}{endpoint}"

    try:
        response = CBIOPORTAL_SESSION.get(
            url,
            params=params,
            timeout=(
                CBIOPORTAL_CONNECT_TIMEOUT,
                CBIOPORTAL_READ_TIMEOUT
            )
        )
    except requests.exceptions.RequestException as error:
        raise HTTPException(
            status_code=503,
            detail=f"Could not connect to cBioPortal: {error}"
        )

    if response.status_code != 200:
        raise HTTPException(
            status_code=response.status_code,
            detail=f"cBioPortal API error: {response.text}"
        )

    return response.json()

def cbioportal_get_all_pages(
    endpoint: str,
    params: Optional[dict] = None,
    page_size: int = CBIOPORTAL_PAGE_SIZE,
    max_pages: int =10000
):
    """retrieve every page from paginated endpoint
    """
    all_records = []
    for page_number in range(max_pages):
        page_params = dict(params or {})
        page_params["pageSize"] = page_size
        page_params["pageNumber"]=page_number
        page_records = cbioportal_get(endpoint,params= page_params)
        if not isinstance(page_records, list): 
            raise HTTPException(status_code = 502, detail = ("Expected a list from cbioportal"
                                                             f"for endpoint '{endpoint}'")
            )
        all_records.extend(page_records)

        if len(page_records) < page_size:
            return all_records
        
    raise HTTPException(status_code = 502, detail = "Pagination reached safety limit of "
                        f"{max_pages} pages for endpoint '{endpoint}'"
    )
    

def clean_study(study: dict):
    """
    Cleans up description.
    """
    return {
        "study_id": study.get("studyId"),
        "name": study.get("name"),
        "description": study.get("description"),
        "cancer_type": study.get("cancerTypeId"),
        "sample_count": study.get("allSampleCount"),
        "citation": study.get("citation"),
        "pmid": study.get("pmid")
    }


def clean_sample(sample: dict):
    """
    Cleaner description
    """
    return {
        "sample_id": sample.get("sampleId"),
        "patient_id": sample.get("patientId"),
        "study_id": sample.get("studyId"),
        "sample_type": sample.get("sampleType"),
        "sequenced": sample.get("sequenced"),
        "copy_number_segment_present": sample.get("copyNumberSegmentPresent")
    }

def clean_clinical_value(record: dict):

    return {
        "study_id": record.get("studyId"),
        "sample_id": record.get("sampleId"),
        "patient_id": record.get("patientId"),
        "clinical_attribute_id": record.get("clinicalAttributeId"),
        "value": record.get("value")
    }

def clean_clinical_attribute(attribute: dict):
    """
    clinical attribute more cleaner format
    """
    return {
        "clinical_attribute_id": attribute.get("clinicalAttributeId"),
        "display_name": attribute.get("displayName"),
        "description": attribute.get("description"),
        "datatype": attribute.get("datatype"),
        "patient_attribute": attribute.get("patientAttribute"),
        "priority": attribute.get("priority")
    }


def rows_to_csv_response(rows: list, filename: str):
    """
    Convert a list of dictionaries into a downloadable CSV response.
    """
    output = io.StringIO()

    if rows:
        fieldnames = []

        # Preserve column order while collecting all possible keys
        for row in rows:
            for key in row.keys():
                if key not in fieldnames:
                    fieldnames.append(key)

        writer = csv.DictWriter(output, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    headers = {
        "Content-Disposition": f'attachment; filename="{filename}"'
    }

    return Response(
        content=output.getvalue(),
        media_type="text/csv",
        headers=headers
    )

def clinical_long_to_wide(records: list, clinical_data_type: str):
    """
    convert clinical data from long to wide format (see above)
    """
    wide_rows = {}

    for record in records: 
        if clinical_data_type == "SAMPLE":
            row_id = record.get("sample_id")

        else:
            row_id = record.get("patient_id")

        if row_id is None:
            continue

        if row_id not in wide_rows:
            wide_rows[row_id] = {
                "study_id": record.get("study_id"),
                "sample_id": record.get("sample_id"),
                "patient_id": record.get("patient_id")
            }

        attribute = record.get("clinical_attribute_id")
        value = record.get("value")

        if attribute: 
            wide_rows[row_id][attribute] = value
    return list(wide_rows.values())


def fetch_clinical_values(
    study_id: str,
    clinical_data_type: str = "SAMPLE",
    attribute_id: Optional[str] = None,
    page_size: int = 10000000,
    page_number: int = 0,
    fetch_all: bool = False
    
):
    """
    Pull clinical data values from cBioPortal for a selected study.
    """
    params = {
        "clinicalDataType": clinical_data_type,
        "projection": "SUMMARY"
    }

    if attribute_id:
        params["attributeId"] = attribute_id

    endpoint = f"/api/studies/{study_id}/clinical-data"
    if fetch_all:
        raw_clinical_data = cbioportal_get_all_pages(
            endpoint,
            params = params
        )
    else: 
        params["pageSize"] = page_size
        params["pageNumber"] = page_number
        raw_clinical_data = cbioportal_get(endpoint, params = params)


    return [clean_clinical_value(record) for record in raw_clinical_data]

def clean_copy_number_segment(segment: dict):
    """
    Convert a raw cBioPortal CN segment record into a cleaner format.
    """
    return {
        "study_id": segment.get("studyId"),
        "sample_id": segment.get("sampleId"),
        "chromosome": segment.get("chromosome"),
        "start": segment.get("start"),
        "end": segment.get("end"),
        "number_of_probes": segment.get("numberOfProbes"),
        "segment_mean": segment.get("segmentMean")
    }


@app.get("/")
def home():
    return {
        "message": "Cancer Genomics Study Explorer API is running",
        "purpose": (
            "This API simplifies public cBioPortal cancer genomics metadata "
            "for easier researcher-facing exploration."
        ),
        "try_these_endpoints": [
            "/studies/clean",
            "/studies/search?keyword=breast",
            "/studies/search?keyword=lung",
            "/studies/{study_id}/summary",
            "/studies/{study_id}/samples/clean"
        ]
    }


@app.get("/alive")
def alive_check():
    return {"status": "running"}


@app.get(
    "/studies/clean",
    tags=["Study Metadata"],
    summary="Get a clean paginated list of cancer genomics studies"
)
def get_clean_studies(
    page_size: int = Query(10, description="Number of studies to return"),
    page_number: int = Query(0, description="Page number starting from 0")
):
    """
    Return a clean paginated list of cancer genomics studies.
    """
    params = {
        "projection": "SUMMARY",
        "pageSize": page_size,
        "pageNumber": page_number
    }

    raw_studies = cbioportal_get("/api/studies", params=params)
    clean_studies = [clean_study(study) for study in raw_studies]

    next_page_number = page_number + 1
    previous_page_number = page_number - 1 if page_number > 0 else None

    return {
        "source": "cBioPortal",
        "page_size": page_size,
        "page_number": page_number,
        "count_on_this_page": len(clean_studies),
        "next_page": f"/studies/clean?page_size={page_size}&page_number={next_page_number}",
        "previous_page": (
            f"/studies/clean?page_size={page_size}&page_number={previous_page_number}"
            if previous_page_number is not None
            else None
        ),
        "studies": clean_studies
    }

@app.get("/studies/search")
def search_studies(
    keyword: str = Query(..., description="Search term, such as breast, lung, melanoma, glioma, TCGA"),
    page_size: int = Query(10, description="Number of matching studies to return per page"),
    page_number: int = Query(0, description="Page number starting from 0")
):
    """
    Search studies by keyword and return simplified paginated results.
    """
    params = {
        "projection": "SUMMARY",
        "pageSize": 1000,
        "pageNumber": 0
    }

    raw_studies = cbioportal_get("/api/studies", params=params)
    keyword_lower = keyword.lower()

    matching_studies = []

    for study in raw_studies:
        searchable_text = " ".join([
            str(study.get("studyId", "")),
            str(study.get("name", "")),
            str(study.get("description", "")),
            str(study.get("cancerTypeId", "")),
            str(study.get("citation", ""))
        ]).lower()

        if keyword_lower in searchable_text:
            matching_studies.append(clean_study(study))

    start_index = page_number * page_size
    end_index = start_index + page_size

    paginated_results = matching_studies[start_index:end_index]

    next_page_number = page_number + 1
    previous_page_number = page_number - 1 if page_number > 0 else None

    total_matches = len(matching_studies)
    total_pages = (total_matches + page_size - 1) // page_size

    return {
        "source": "cBioPortal",
        "keyword": keyword,
        "page_size": page_size,
        "page_number": page_number,
        "total_matches": total_matches,
        "total_pages": total_pages,
        "count_on_this_page": len(paginated_results),
        "next_page": (
            f"/studies/search?keyword={keyword}&page_size={page_size}&page_number={next_page_number}"
            if page_number + 1 < total_pages
            else None
        ),
        "previous_page": (
            f"/studies/search?keyword={keyword}&page_size={page_size}&page_number={previous_page_number}"
            if previous_page_number is not None
            else None
        ),
        "matching_studies": paginated_results
    }

@app.get(
    "/study-ids",
    tags=["Study Search"],
    summary="Find clean study IDs by keyword"
)
def find_study_ids(
    keyword: str = Query(..., description="Search term, such as breast, lung, melanoma, TCGA"),
    max_results: int = Query(20, description="Maximum number of study IDs to return")
):
    """
    Search cBioPortal studies and return only the most useful fields:
    study_id, name, cancer_type, and sample_count.
    """

    params = {
        "projection": "SUMMARY",
        "pageSize": 1000,
        "pageNumber": 0
    }

    raw_studies = cbioportal_get("/api/studies", params=params)
    keyword_lower = keyword.lower()

    matches = []

    for study in raw_studies:
        searchable_text = " ".join([
            str(study.get("studyId", "")),
            str(study.get("name", "")),
            str(study.get("description", "")),
            str(study.get("cancerTypeId", "")),
            str(study.get("citation", ""))
        ]).lower()

        if keyword_lower in searchable_text:
            matches.append({
                "study_id": study.get("studyId"),
                "name": study.get("name"),
                "cancer_type": study.get("cancerTypeId"),
                "sample_count": study.get("allSampleCount")
            })

    return {
        "keyword": keyword,
        "count": len(matches[:max_results]),
        "study_ids": matches[:max_results]
    }

@app.get("/studies/{study_id}/summary")
def get_study_summary(study_id: str):
    """
    Return simplified metadata for one study.
    Example: /studies/brca_tcga/summary
    """
    raw_study = cbioportal_get(
        f"/api/studies/{study_id}",
        params={"projection": "DETAILED"}
    )

    return {
        "source": "cBioPortal",
        "study": clean_study(raw_study)
    }


@app.get("/studies/{study_id}/samples/clean")
def get_clean_samples_for_study(
    study_id: str,
    page_size: int = Query(20, description="Number of samples to return"),
    page_number: int = Query(0, description="Page number starting from 0")
):
    """
    Return clean sample metadata for a selected cancer genomics study.
    Example: /studies/brca_tcga/samples/clean
    """
    params = {
        "projection": "SUMMARY",
        "pageSize": page_size,
        "pageNumber": page_number
    }

    raw_samples = cbioportal_get(
        f"/api/studies/{study_id}/samples",
        params=params
    )

    clean_samples = [clean_sample(sample) for sample in raw_samples]

    return {
        "source": "cBioPortal",
        "study_id": study_id,
        "count": len(clean_samples),
        "samples": clean_samples
    }

@app.get(
    "/studies/{study_id}/clinical-data/long",
    tags=["Clinical Data"],
    summary="Get clinical data in long format"
)
def get_clinical_data_long(
    study_id: str,
    clinical_data_type: Literal["SAMPLE", "PATIENT"] = Query(
        "SAMPLE",
        description="Use SAMPLE for sample-level data or PATIENT for patient-level data"
    ),
    attribute_id: Optional[str] = Query(
        None,
        description="Optional clinical attribute ID, such as AGE, SEX, OS_STATUS"
    ),
    page_size: int = Query(10000000, ge=1, le=10000000),
    page_number: int = Query(0, ge=0)
):
    clinical_values = fetch_clinical_values(
        study_id=study_id,
        clinical_data_type=clinical_data_type,
        attribute_id=attribute_id,
        page_size=page_size,
        page_number=page_number
    )

    return {
        "source": "cBioPortal",
        "study_id": study_id,
        "clinical_data_type": clinical_data_type,
        "attribute_id": attribute_id,
        "format": "long",
        "count": len(clinical_values),
        "clinical_data": clinical_values
    }

@app.get(
    "/studies/{study_id}/clinical-data/wide",
    tags=["Clinical Data"],
    summary="Get clinical data in wide format"
)
def get_clinical_data_wide(
    study_id: str,
    clinical_data_type: Literal["SAMPLE", "PATIENT"] = Query(
        "SAMPLE",
        description="Use SAMPLE for sample-level data or PATIENT for patient-level data"
    ),
    attribute_id: Optional[str] = Query(
        None, 
        description="Optional clinical attribute ID, such as AGE, SEX, or OS_STATUS"
    ),
    page_size: int = Query(1000000, ge=1, le=10000000),
    page_number: int = Query(0, ge=0)
):
    clinical_values = fetch_clinical_values(
        study_id=study_id,
        clinical_data_type=clinical_data_type,
        attribute_id=attribute_id,
        page_size=page_size,
        page_number=page_number,
    )
    wide_data = clinical_long_to_wide(
        records=clinical_values,
        clinical_data_type=clinical_data_type
    )
    return {
        "source": "cBioPortal",
        "study_id": study_id,
        "clinical_data_type": clinical_data_type,
        "attribute_id": attribute_id,
        "format": "wide",
        "count": len(wide_data),
        "clinical_data": wide_data
    }

@app.get(
    "/studies/{study_id}/clinical-data/export.csv",
    tags=["Clinical Data"],
    summary="Export clinical data as a CSV file"
)
def export_clinical_data_csv(
    study_id: str,
    clinical_data_type: Literal["SAMPLE", "PATIENT"] = Query(
        "SAMPLE",
        description="Use SAMPLE for sample-level data or PATIENT for patient-level data"
    ),
    table_format: Literal["long", "wide"] = Query(
        "wide",
        alias="format",
        description="Export as long or wide table format"
    ),
    attribute_id: Optional[str] = Query(
        None,
        description="Optional clinical attribute ID, such as AGE, SEX, OS_STATUS"
    ),
    page_size: int = Query(10000000, ge=1, le=10000000),
    page_number: int = Query(0, ge=0)
):
    clinical_values = fetch_clinical_values(
        study_id=study_id,
        clinical_data_type=clinical_data_type,
        attribute_id=attribute_id,
        page_size=page_size,
        page_number=page_number
    )

    if table_format == "wide":
        rows = clinical_long_to_wide(
            records=clinical_values,
            clinical_data_type=clinical_data_type
        )
    else:
        rows = clinical_values

    filename = f"{study_id}_{clinical_data_type.lower()}_clinical_{table_format}.csv"

    return rows_to_csv_response(rows, filename)

@app.get(
    "/studies/{study_id}/samples/{sample_id}/copy-number-segments/clean",
    tags=["Copy Number Segments"],
    summary="Get clean copy-number segments for one sample"
)
def get_clean_copy_number_segments(
    study_id: str,
    sample_id: str,
    chromosome: Optional[str] = Query(
        None,
        description="Optional chromosome filter, such as 1, 8, 17, X"
    ),
    page_size: int = Query(20000, ge=1, le=20000),
    page_number: int = Query(0, ge=0)
):
    params = {
        "projection": "SUMMARY",
        "pageSize": page_size,
        "pageNumber": page_number
    }

    if chromosome:
        params["chromosome"] = chromosome

    raw_segments = cbioportal_get(
        f"/api/studies/{study_id}/samples/{sample_id}/copy-number-segments",
        params=params
    )

    clean_segments = [clean_copy_number_segment(segment) for segment in raw_segments]

    return {
        "source": "cBioPortal",
        "study_id": study_id,
        "sample_id": sample_id,
        "chromosome": chromosome,
        "count": len(clean_segments),
        "copy_number_segments": clean_segments
    }

@app.get(
    "/studies/{study_id}/clinical-attributes/clean",
    tags=["Clinical Data"],
    summary="get clean clinical attribute definition for a study"
)
def get_clean_clinical_attributes(study_id: str):
    """
    return cleaned attribute for a study to help understand what's available 
    """
    raw_attributes = cbioportal_get(
        f"/api/studies/{study_id}/clinical-attributes",
        params={"projection": "DETAILED"}
    )

    clean_attributes = [
        clean_clinical_attribute(attribute)
        for attribute in raw_attributes
    ]

    return {
        "source": "cBioPortal",
        "study_id": study_id,
        "count": len(clean_attributes),
        "clinical_attributes": clean_attributes
    }

def clean_molecular_profile(profile: dict):
    """
    Convert one raw molecular profile into a cleaner format.
    """
    alteration_type = profile.get("molecularAlterationType") or profile.get("geneticAlterationType")

    return {
        "molecular_profile_id": profile.get("molecularProfileId"),
        "study_id": profile.get("studyId"),
        "name": profile.get("name"),
        "description": profile.get("description"),
        "genetic_alteration_type": alteration_type,
        "datatype": profile.get("datatype"),
        "show_profile_in_analysis_tab": profile.get(
            "showProfileInAnalysisTab"
        ),
        "patient_level": profile.get("patientLevel")
    }

@app.get(
    "/studies/{study_id}/molecular-profiles/clean",
    tags=["Molecular Profiles"],
    summary="Get clean molecular profiles for a study"
)
def get_clean_molecular_profiles(study_id: str):
    """
    Return available molecular profiles for a selected study.
    Useful for finding mutation, CNA, expression, or other data profile IDs.
    """
    raw_profiles = cbioportal_get(
        f"/api/studies/{study_id}/molecular-profiles",
        params={"projection": "DETAILED"}
    )

    clean_profiles = [
        clean_molecular_profile(profile)
        for profile in raw_profiles
    ]

    grouped_profiles = {}

    for profile in clean_profiles:
        alteration_type = profile.get("genetic_alteration_type") or "UNKNOWN"

        if alteration_type not in grouped_profiles:
            grouped_profiles[alteration_type] = []

        grouped_profiles[alteration_type].append(profile)

    return {
        "source": "cBioPortal",
        "study_id": study_id,
        "count": len(clean_profiles),
        "profiles_by_alteration_type": grouped_profiles,
        "molecular_profiles": clean_profiles
    }

def clean_sample_list(sample_list: dict):
    """
    Convert one raw sample list into a cleaner format.
    """
    return {
        "sample_list_id": sample_list.get("sampleListId"),
        "study_id": sample_list.get("studyId"),
        "category": sample_list.get("category"),
        "name": sample_list.get("name"),
        "description": sample_list.get("description"),
        "sample_count": sample_list.get("sampleCount")
    }

@app.get(
    "/studies/{study_id}/sample-lists/clean",
    tags=["Sample Lists"],
    summary="Get clean sample lists for a study"
)
def get_clean_sample_lists(study_id: str):
    """
    Return available sample lists for a selected study.
    Sample lists define which samples belong to subsets such as sequenced samples,
    CNA-profiled samples, or all tumor samples.
    """
    raw_sample_lists = cbioportal_get(
        f"/api/studies/{study_id}/sample-lists",
        params={"projection": "DETAILED"}
    )

    clean_lists = [
        clean_sample_list(sample_list)
        for sample_list in raw_sample_lists
    ]

    return {
        "source": "cBioPortal",
        "study_id": study_id,
        "count": len(clean_lists),
        "sample_lists": clean_lists
    }

def group_molecular_profiles_by_type(molecular_profiles: list):
    """
    group molecular profiels by genetic alteration type.
    """
    grouped_profiles = {}

    for profile in molecular_profiles: 
        alteration_type = profile.get("genetic_alteration_type") or "UNKNOWN"

        if alteration_type not in grouped_profiles:
            grouped_profiles[alteration_type] = []

        grouped_profiles[alteration_type].append(profile)
    return grouped_profiles

def group_sample_lists_by_category(sample_lists: list):
    """
    group samples lists by sample list category.
    """
    grouped_lists = {}

    for sample_list in sample_lists:
        category = sample_list.get("category") or "UNKNOWN"
        if category not in grouped_lists:
            grouped_lists[category] = []
        grouped_lists[category].append(sample_list)
    return grouped_lists

def build_study_overview(study_id: str):
    """
    Build a high lvl overview of one study by combining study metadata, clinical attributes, molecular profiles and sample list
    """
    raw_study = cbioportal_get(
        f"/api/studies/{study_id}",
        params={"projection":"DETAILED"}

    )

    raw_clinical_attributes = cbioportal_get(
        f"/api/studies/{study_id}/clinical-attributes",
        params={"projection":"DETAILED"}
    )

    raw_molecular_profiles = cbioportal_get(
        f"/api/studies/{study_id}/molecular-profiles",
        params={"projection":"DETAILED"}
    )
    
    raw_sample_lists = cbioportal_get(
        f"/api/studies/{study_id}/sample-lists",
        params={"projection": "DETAILED"}

    )

    study = clean_study(raw_study)
    clinical_attributes = [
        clean_clinical_attribute(attribute)
        for attribute in raw_clinical_attributes
    ]

    molecular_profiles = [
        clean_molecular_profile(profile)
        for profile in raw_molecular_profiles
    ]

    sample_lists = [
        clean_sample_list(sample_list)
        for sample_list in raw_sample_lists
    ]

    available_data_types = sorted(
        list(
            set(
                profile.get("genetic_alteration_type")
                for profile in molecular_profiles 
                if profile.get("genetic_alteration_type")
            )
        )
    )
    return {
        "study": study,
        "counts": {
            "sample_count": study.get("sample_count"),
            "clinical_attribute_count": len(clinical_attributes),
            "molecular_profile_count": len(molecular_profiles),
            "sample_list_count": len(sample_lists)

        },
        "available_data_types": available_data_types,
        "clinical_atttributes_preview": clinical_attributes[:15],
        "molecular_profiles_by_type": group_molecular_profiles_by_type(molecular_profiles),
        "sample_lists_by_category": group_sample_lists_by_category(sample_lists),
        "recommended_next_steps" :[
            f"/studies/{study_id}/clinical-data/wide?clinical_data_type=SAMPLE",
            f"/studies/{study_id}/clinical-data/export.csv?clinical_data_type=SAMPLE&format=wide",
            f"/studies/{study_id}/molecular-profiles/clean",
            f"/studies/{study_id}/sample-lists/clean"
        ]
    }

@app.get("/studies/{study_id}/overview",
    tags=["Study Overview"],
    summary = "Get a high-level overview of one cancer genomics study"
)
def get_study_overview(study_id: str):
    """
    return a guided overview of one study, profiles, attributes, etc
    """
    overview = build_study_overview(study_id)

    return{
        "source": "cBioPortal",
        "study_id": study_id,
        "overview": overview
    }

@app.get(
    "/studies/overviews",
    tags=["Study Overview"],
    summary="Get high-level overviews for multiple cancer genomics studies"

)
def get_multiple_study_overviews(
    study_ids: str = Query(
        ...,
        description = "Comma-separated study IDs, such as brca_tcga, paad_tcga, cll_broad_2022"

    )
):
    """
    return overview summaries for multiple studies
    """
    study_id_list = [
        study_id.strip()
        for study_id in study_ids.split(",")
        if study_id.strip()

    ]

    overviews = []
    errors = []

    for study_id in study_id_list:
        try: 
            overview= build_study_overview(study_id)
            overviews.append({
                "study_id": study_id,
                "overview": overview
            })
        except HTTPException as error:
            errors.append({
                "study_id":study_id,
                "status_code": error.status_code,
                "detail": error.detail

            })
    return {
        "source": "cBioPortal",
        "requested_study_count": len(study_id_list),
        "successful_study_count": len(overviews),
        "error_count": len(errors),
        "overviews": overviews,
        "errors": errors
    }


def clean_gene(gene: dict):
    """
    convert one raw cBioPortal gene record into cleaner format
    """
    return {
        "entrez_gene_id": gene.get("entrezGeneId"),
        "hugo_gene_symbol": gene.get("hugoGeneSymbol"),
        "gene_type": gene.get("type"),
        "cytoband": gene.get("cytoband"),
        "length": gene.get("length")
    }


def clean_mutation(mutation: dict):
    """
    convert the mutational record into cleaner format
    """
    gene = mutation.get("gene") or {}

    return {
        "molecular_profile_id": mutation.get("molecularProfileId"),
        "sample_id": mutation.get("sampleId"),
        "patient_id": mutation.get("patientId"),
        "entrez_gene_id": mutation.get("entrezGeneId"),
        "gene_symbol": gene.get("hugoGeneSymbol"),
        "mutation_type": mutation.get("mutationType"),
        "protein_change": mutation.get("proteinChange"),
        "chromosome": mutation.get("chromosome"),
        "start_position": mutation.get("startPosition"),
        "end_position": mutation.get("endPosition"),
        "reference_allele": mutation.get("referenceAllele"),
        "variant_allele": mutation.get("variantAllele"),
        "variant_type": mutation.get("variantType"),
        "ncbi_build": mutation.get("ncbiBuild"),
        "tumor_alt_count": mutation.get("tumorAltCount"),
        "tumor_ref_count": mutation.get("tumorRefCount"),
        "keyword": mutation.get("keyword")
    }


def find_mutation_profile_id(study_id: str):
    """
    find a mutation molecular profile for a selected study
    """
    raw_profiles = cbioportal_get(
        f"/api/studies/{study_id}/molecular-profiles",
        params={"projection": "DETAILED"}
    )

    mutation_profiles = [
        profile
        for profile in raw_profiles
        if (
            profile.get("molecularAlterationType")
            or profile.get("geneticAlterationType")
        ) == "MUTATION_EXTENDED"
    ]

    if not mutation_profiles:
        available_profiles = [
            {
                "molecular_profile_id": profile.get("molecularProfileId"),
                "genetic_alteration_type":(
                    profile.get("molecularAlterationType") or profile.get("geneticAlterationType") 
                ), 
                "datatype": profile.get("datatype") 
            }
            for profile in raw_profiles
        ]

        raise HTTPException(
            status_code=404,
            detail={
                "message": (
                    f"No mutation molecular profile was found "
                    f"for study '{study_id}'."
                ),
                "available_profiles": available_profiles,
                "suggestion": (
                    "Check /studies/{study_id}/molecular-profiles/clean "
                    "or choose a study containing MUTATION_EXTENDED data."
                )
            }
        )

    preferred_profiles = [
        profile
        for profile in mutation_profiles
        if profile.get("showProfileInAnalysisTab") is True
    ]

    selected_profile = (
        preferred_profiles[0]
        if preferred_profiles
        else mutation_profiles[0]
    )

    return selected_profile.get("molecularProfileId")


def find_mutation_sample_list_id(study_id: str):
    """
    find most approrpiate sample list for mutation analysis
    """
    raw_sample_lists = cbioportal_get(
        f"/api/studies/{study_id}/sample-lists",
        params={"projection": "DETAILED"}
    )

    sequenced_lists = [
        sample_list
        for sample_list in raw_sample_lists
        if sample_list.get("sampleListId", "").lower().endswith("_sequenced")
    ]

    mutation_lists = [
        sample_list
        for sample_list in raw_sample_lists
        if "mutation" in str(sample_list.get("category", "")).lower()
    ]

    all_sample_lists = [
        sample_list
        for sample_list in raw_sample_lists
        if sample_list.get("sampleListId", "").lower().endswith("_all")
    ]

    if sequenced_lists:
        selected_list = sequenced_lists[0]

    elif mutation_lists:
        selected_list = mutation_lists[0]

    elif all_sample_lists:
        selected_list = all_sample_lists[0]

    elif raw_sample_lists:
        selected_list = raw_sample_lists[0]

    else:
        raise HTTPException(
            status_code=404,
            detail=f"No sample list was found for study '{study_id}'."
        )

    return selected_list.get("sampleListId")


def get_gene_record(gene_symbol: str):
    """
    look up using HUGO gene symbol
    """
    raw_gene = cbioportal_get(
        f"/api/genes/{gene_symbol.upper()}"
    )

    return clean_gene(raw_gene)


def fetch_gene_mutations(
    study_id: str,
    gene_symbol: str,
    molecular_profile_id: Optional[str] = None,
    sample_list_id: Optional[str] = None,
    page_size: int = 1000,
    page_number: int = 0,
    fetch_all: bool = False
):
    """
    pull mutation records for one gene in one study
    """
    gene = get_gene_record(gene_symbol)

    selected_profile_id = (
        molecular_profile_id
        if molecular_profile_id
        else find_mutation_profile_id(study_id)
    )

    selected_sample_list_id = (
        sample_list_id
        if sample_list_id
        else find_mutation_sample_list_id(study_id)
    )

    params = {
        "sampleListId": selected_sample_list_id,
        "entrezGeneId": gene.get("entrez_gene_id"),
        "projection": "DETAILED",
    }

    mutation_endpoint = (f"/api/molecular-profiles/{selected_profile_id}/mutations")
    if fetch_all:
        raw_mutations = cbioportal_get_all_pages(
            mutation_endpoint,
            params = params
        )
    else: 
        params["pageSize"] = page_size
        params["pageNumber"] = page_number

        raw_mutations = cbioportal_get(mutation_endpoint, params=params)

    clean_mutations = [
        clean_mutation(mutation)
        for mutation in raw_mutations
    ]

    return {
        "gene": gene,
        "molecular_profile_id": selected_profile_id,
        "sample_list_id": selected_sample_list_id,
        "mutations": clean_mutations
    }

def build_gene_mutation_summary(
        study_id: str,
        gene_symbol: str,
        molecular_profile_id: str,
        sample_list_id: str,
        profiled_sample_ids: list
):
    """
    build mutation summary for one gene using an already 
    selected molecular profile + sample list
    """

    mutation_data = fetch_gene_mutations(
        study_id=study_id,
        gene_symbol=gene_symbol,
        molecular_profile_id=molecular_profile_id,
        sample_list_id=sample_list_id,
        fetch_all = True
    )

    mutations = mutation_data["mutations"]

    mutated_sample_ids = sorted(
        {
            mutation.get("sample_id")
            for mutation in mutations
            if mutation.get("sample_id")
        }
    )

    mutation_type_counts = {}

    for mutation in mutations: 
        mutation_type = mutation.get("mutation_type") or "UNKNOWN"

        mutation_type_counts[mutation_type] = (
            mutation_type_counts.get(mutation_type, 0) + 1
        )
    profiled_sample_count = len(profiled_sample_ids)
    mutated_sample_count = len(mutated_sample_ids)

    mutation_frequency= (
        mutated_sample_count / profiled_sample_count
        if profiled_sample_count > 0
        else None
    )

    return {
        "gene": mutation_data["gene"],
        "profiled_sample_count": profiled_sample_count,
        "mutated_sample_count": mutated_sample_count,
        "mutation_record_count": len(mutations),
        "mutation_frequency": mutation_frequency,
        "mutation_frequency_percent": round(mutation_frequency *100,2) if mutation_frequency is not None else None,
        "mutation_type_counts": mutation_type_counts,
        "mutated_sample_ids": mutated_sample_ids
    }

def parse_gene_symbols(genes: str):
    """
    convert a comma separated gene into clean list of unique uppercase gene symbols
    """
    gene_symbols = []
    seen_genes = set()

    for gene in genes.split(","):
        gene_symbol = gene.strip().upper()

        if gene_symbol and gene_symbol not in seen_genes:
            gene_symbols.append(gene_symbol)
            seen_genes.add(gene_symbol)

    if not gene_symbols:
        raise HTTPException(
            status_code=400,
            detail="at least one gene symbol must be provided"
        )

    if len(gene_symbols) > 25:
        raise HTTPException(
            status_code=400,
            detail="max of 25 gene symbols may be requested at once"
        )

    return gene_symbols

def build_mutation_panel_summary(
        study_id: str,
        gene_symbols: list,
        molecular_profile_id: Optional[str] = None,
        sample_list_id: Optional[str] = None
):
    """
    build mutation summaries for multiple genes in one study
    """
    selected_profile_id = (
        molecular_profile_id 
        if molecular_profile_id
        else find_mutation_profile_id(study_id)
    )

    selected_sample_list_id = (
        sample_list_id
        if sample_list_id else find_mutation_sample_list_id(study_id)
    )

    profiled_sample_ids = cbioportal_get(
        f"/api/sample-lists/{selected_sample_list_id}/sample-ids"
    )
    
    summaries = []
    errors = []

    for gene_symbol in gene_symbols: 
        try: 
            gene_summary = build_gene_mutation_summary(
                study_id=study_id,
                gene_symbol=gene_symbol,
                molecular_profile_id=selected_profile_id,
                sample_list_id=selected_sample_list_id,
                profiled_sample_ids=profiled_sample_ids
            )

            summaries.append(gene_summary)

        except HTTPException as error:
            errors.append({
                "gene_symbol": gene_symbol,
                "status_code": error.status_code,
                "detail": error.detail
            })

    summaries.sort(
        key=lambda summary: (
            summary.get("mutation_frequency")
            if summary.get("mutation_frequency") is not None
            else -1
    ),
        reverse=True
    )
    return {
        "molecular_profile_id": selected_profile_id,
        "sample_list_id": selected_sample_list_id,
        "profiled_sample_count": len(profiled_sample_ids),
        "summaries": summaries,
        "errors": errors
    }


def build_sample_gene_mutation_matrix(
    study_id: str,
    gene_symbols: list,
    molecular_profile_id: Optional[str] = None,
    sample_list_id: Optional[str] = None
):
    """
    build sample-by-gene mutation matrix for one study
    """
    selected_profile_id = (
        molecular_profile_id
        if molecular_profile_id
        else find_mutation_profile_id(study_id)
    )

    selected_sample_list_id = (
        sample_list_id
        if sample_list_id
        else find_mutation_sample_list_id(study_id)
    )

    profiled_sample_ids = cbioportal_get(
        f"/api/sample-lists/{selected_sample_list_id}/sample-ids"
    )

    matrix_by_sample = {
        sample_id: {
            "study_id": study_id,
            "sample_id": sample_id
        }
        for sample_id in profiled_sample_ids
    }

    successful_genes = []
    errors = []

    for requested_gene_symbol in gene_symbols:
        try:
            mutation_data = fetch_gene_mutations(
                study_id=study_id,
                gene_symbol=requested_gene_symbol,
                molecular_profile_id=selected_profile_id,
                sample_list_id=selected_sample_list_id,
                fetch_all = True
            )

            gene = mutation_data.get("gene") or {}

            gene_symbol = (
                gene.get("hugo_gene_symbol")
                or requested_gene_symbol
            ).upper()

            successful_genes.append(gene_symbol)

            mutation_counts_by_sample = {}
            protein_changes_by_sample = {}

            for mutation in mutation_data["mutations"]:
                sample_id = mutation.get("sample_id")

                # ignore mutation records outside selected samplelist
                if sample_id not in matrix_by_sample:
                    continue

                mutation_counts_by_sample[sample_id] = (
                    mutation_counts_by_sample.get(sample_id, 0) + 1
                )

                protein_change = mutation.get("protein_change")

                if protein_change:
                    if sample_id not in protein_changes_by_sample:
                        protein_changes_by_sample[sample_id] = set()

                    protein_changes_by_sample[sample_id].add(
                        protein_change
                    )

            for sample_id, row in matrix_by_sample.items():
                mutation_count = mutation_counts_by_sample.get(
                    sample_id,
                    0
                )

                protein_changes = sorted(
                    protein_changes_by_sample.get(
                        sample_id,
                        set()
                    )
                )

                row[f"{gene_symbol}_mutated"] = mutation_count > 0

                row[f"{gene_symbol}_mutation_count"] = mutation_count

                row[f"{gene_symbol}_protein_changes"] = (
                    "; ".join(protein_changes)
                    if protein_changes
                    else None
                )

        except HTTPException as error:
            errors.append({
                "gene_symbol": requested_gene_symbol,
                "status_code": error.status_code,
                "detail": error.detail
            })

    matrix_rows = [
        matrix_by_sample[sample_id]
        for sample_id in sorted(matrix_by_sample)
    ]

    return {
        "molecular_profile_id": selected_profile_id,
        "sample_list_id": selected_sample_list_id,
        "profiled_sample_count": len(profiled_sample_ids),
        "requested_genes": gene_symbols,
        "successful_genes": successful_genes,
        "matrix": matrix_rows,
        "errors": errors
    }

def build_clinical_genomic_analysis_table(
    study_id: str,
    gene_symbols: list,
    molecular_profile_id: Optional[str] = None,
    sample_list_id: Optional[str] = None
):
    """
    Combine sample metadata, sample-level clinical data,
    patient-level clinical data, and mutation data.
    """
    mutation_matrix_data = build_sample_gene_mutation_matrix(
        study_id=study_id,
        gene_symbols=gene_symbols,
        molecular_profile_id=molecular_profile_id,
        sample_list_id=sample_list_id
    )

    # Get sample metadata so each sample can be linked
    # to the correct patient.
    raw_samples = cbioportal_get(
        f"/api/studies/{study_id}/samples",
        params={
            "projection": "SUMMARY"
        }
    )

    patient_id_by_sample = {
        sample.get("sampleId"): sample.get("patientId")
        for sample in raw_samples
        if sample.get("sampleId")
    }

    # Sample-level clinical data
    sample_clinical_values = fetch_clinical_values(
        study_id=study_id,
        clinical_data_type="SAMPLE",
        fetch_all = True
    )

    sample_clinical_rows = clinical_long_to_wide(
        records=sample_clinical_values,
        clinical_data_type="SAMPLE"
    )

    sample_clinical_by_sample = {
        row.get("sample_id"): row
        for row in sample_clinical_rows
        if row.get("sample_id")
    }

    # Patient-level clinical data
    patient_clinical_values = fetch_clinical_values(
        study_id=study_id,
        clinical_data_type="PATIENT",
        fetch_all = True
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

    analysis_rows = []
    samples_without_sample_clinical_data = []
    samples_without_patient_clinical_data = []

    for mutation_row in mutation_matrix_data["matrix"]:
        sample_id = mutation_row.get("sample_id")

        sample_clinical_row = sample_clinical_by_sample.get(
            sample_id
        )

        patient_id = patient_id_by_sample.get(sample_id)

        if sample_clinical_row and sample_clinical_row.get(
            "patient_id"
        ):
            patient_id = sample_clinical_row.get("patient_id")

        patient_clinical_row = patient_clinical_by_patient.get(
            patient_id
        )

        analysis_row = {
            "study_id": study_id,
            "sample_id": sample_id,
            "patient_id": patient_id
        }

        # Add sample-level clinical columns
        if sample_clinical_row:
            for column_name, value in sample_clinical_row.items():
                if column_name not in {
                    "study_id",
                    "sample_id",
                    "patient_id"
                }:
                    analysis_row[column_name] = value
        else:
            samples_without_sample_clinical_data.append(
                sample_id
            )

        # Add patient-level clinical columns
        if patient_clinical_row:
            for column_name, value in patient_clinical_row.items():
                if column_name in {
                    "study_id",
                    "sample_id",
                    "patient_id"
                }:
                    continue

                # Avoid overwriting a sample-level column
                # with the same name.
                if column_name in analysis_row:
                    analysis_row[
                        f"patient_{column_name}"
                    ] = value
                else:
                    analysis_row[column_name] = value
        else:
            samples_without_patient_clinical_data.append(
                sample_id
            )

        # Add mutation columns
        for column_name, value in mutation_row.items():
            if column_name not in {
                "study_id",
                "sample_id"
            }:
                analysis_row[column_name] = value

        analysis_rows.append(analysis_row)

    return {
        "molecular_profile_id": mutation_matrix_data[
            "molecular_profile_id"
        ],
        "sample_list_id": mutation_matrix_data[
            "sample_list_id"
        ],
        "requested_genes": mutation_matrix_data[
            "requested_genes"
        ],
        "successful_genes": mutation_matrix_data[
            "successful_genes"
        ],
        "profiled_sample_count": mutation_matrix_data[
            "profiled_sample_count"
        ],
        "clinical_sample_count": len(
            sample_clinical_by_sample
        ),
        "clinical_patient_count": len(
            patient_clinical_by_patient
        ),
        "analysis_row_count": len(analysis_rows),
        "samples_without_clinical_data": (
            samples_without_sample_clinical_data
        ),
        "samples_without_patient_clinical_data": (
            samples_without_patient_clinical_data
        ),
        "analysis_table": analysis_rows,
        "errors": mutation_matrix_data["errors"]
    }


def is_missing_value(value):
    """
    Conservatively identify blank or explicitly unavailable values.
    Values such as 'Unknown' are preserved as real categories.
    """
    if value is None:
        return True

    if isinstance(value, str):
        return value.strip().casefold() in {
            "",
            "na",
            "n/a",
            "nan",
            "none",
            "null",
            "not available",
            "[not available]"
        }

    return False


def build_analysis_quality_summary(
    rows: list,
    sparse_threshold_percent: float = 50.0
):
    """
    Summarize missingness, uniqueness, numeric ranges,
    and mutation frequencies for an analysis table.
    """
    if not rows:
        return {
            "row_count": 0,
            "column_count": 0,
            "complete_row_count": 0,
            "sparse_threshold_percent": sparse_threshold_percent,
            "sparse_columns": [],
            "columns": []
        }

    column_names = []

    for row in rows:
        for column_name in row:
            if column_name not in column_names:
                column_names.append(column_name)

    column_summaries = []
    sparse_columns = []

    for column_name in column_names:
        values = [
            row.get(column_name)
            for row in rows
        ]

        non_missing_values = [
            value
            for value in values
            if not is_missing_value(value)
        ]

        missing_count = (
            len(values) - len(non_missing_values)
        )

        missing_percent = round(
            (missing_count / len(values)) * 100,
            2
        )

        column_summary = {
            "column_name": column_name,
            "non_missing_count": len(
                non_missing_values
            ),
            "missing_count": missing_count,
            "missing_percent": missing_percent,
            "unique_value_count": len({
                str(value)
                for value in non_missing_values
            })
        }

        numeric_values = []

        all_non_missing_values_are_numeric = bool(
            non_missing_values
        )

        for value in non_missing_values:
            # Prevent True and False from being treated
            # as the numbers 1 and 0.
            if isinstance(value, bool):
                all_non_missing_values_are_numeric = False
                break

            try:
                numeric_values.append(float(value))
            except (TypeError, ValueError):
                all_non_missing_values_are_numeric = False
                break

        if all_non_missing_values_are_numeric:
            column_summary["numeric_summary"] = {
                "minimum": min(numeric_values),
                "maximum": max(numeric_values),
                "mean": round(
                    mean(numeric_values),
                    4
                ),
                "median": median(numeric_values)
            }

        if column_name.endswith("_mutated"):
            mutated_count = sum(
                value is True
                for value in non_missing_values
            )

            not_mutated_count = sum(
                value is False
                for value in non_missing_values
            )

            mutation_denominator = (
                mutated_count + not_mutated_count
            )

            column_summary["mutation_summary"] = {
                "mutated_sample_count": mutated_count,
                "not_mutated_sample_count": (
                    not_mutated_count
                ),
                "mutation_frequency": (
                    mutated_count / mutation_denominator
                    if mutation_denominator > 0
                    else None
                )
            }

        if missing_percent >= sparse_threshold_percent:
            sparse_columns.append(column_name)

        column_summaries.append(column_summary)

    complete_row_count = sum(
        all(
            not is_missing_value(
                row.get(column_name)
            )
            for column_name in column_names
        )
        for row in rows
    )

    return {
        "row_count": len(rows),
        "column_count": len(column_names),
        "complete_row_count": complete_row_count,
        "sparse_threshold_percent": (
            sparse_threshold_percent
        ),
        "sparse_columns": sparse_columns,
        "columns": column_summaries
    }


@app.get(
    "/studies/{study_id}/analysis-table/quality",
    tags=["Analysis Tables"],
    summary="Summarize analysis-table completeness and quality"
)
def get_analysis_table_quality(
    study_id: str,
    genes: str = Query(
        ...,
        description=(
            "Comma-separated gene symbols, such as "
            "BRAF,NRAS,NF1,KIT"
        )
    ),
    molecular_profile_id: Optional[str] = Query(None),
    sample_list_id: Optional[str] = Query(None),
    sparse_threshold_percent: float = Query(
        50.0,
        ge=0,
        le=100,
        description=(
            "A column is labeled sparse when its missing "
            "percentage is at least this value."
        )
    )
):
    gene_symbols = parse_gene_symbols(genes)

    analysis_data = build_clinical_genomic_analysis_table(
        study_id=study_id,
        gene_symbols=gene_symbols,
        molecular_profile_id=molecular_profile_id,
        sample_list_id=sample_list_id
    )

    quality_summary = build_analysis_quality_summary(
        rows=analysis_data["analysis_table"],
        sparse_threshold_percent=(
            sparse_threshold_percent
        )
    )

    return {
        "source": "cBioPortal",
        "study_id": study_id,
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
        "quality_summary": quality_summary,
        "errors": analysis_data["errors"]
    }



def build_filtered_clinical_genomic_cohort(
        study_id: str,
        gene_symbols: list,
        molecular_profile_id: Optional[str] = None,
        sample_list_id: Optional[str] = None,
        mutation_gene: Optional[str] = None,
        mutation_status: Optional[str] = None,
        clinical_attribute: Optional[str] = None,
        clinical_value: Optional[str] = None,
        min_numeric_value: Optional[float] = None,
        max_numeric_value: Optional[float] = None
): 
    """
    Build a clinical-genomic analysis table and gilter by mutation status or clinical attribute
    """
    normalized_mutation_gene = mutation_gene.strip().upper() if mutation_gene else None

    if mutation_status and not normalized_mutation_gene:
        raise HTTPException(status_code = 400, detail = ("mutation_gene must be provided when mutation status used"))
    
    if normalized_mutation_gene and not mutation_status: 
        raise HTTPException(status_code=400, detail = ("mutationstatus must be provided when mutation_gene is used"))
    
    if mutation_status not in {None, "mutated", "not_mutated"}:
        raise HTTPException(status_code = 400, detail = ("mutation_status must be either mutated or not_mutated"))
    
    clinical_filter_requested = clinical_value is not None or min_numeric_value is not None or max_numeric_value is not None
    if clinical_filter_requested and not clinical_attribute: 
        raise HTTPException(status_code = 400, detail = ("clinical attribute must be provided when using a clinical filter"))
    
    if clinical_value is not None and (min_numeric_value is not None or max_numeric_value is not None):
        raise HTTPException(status_code = 400, detail = ("use either clinical_value or numeric minimum/maximum filters, not both"))

    if min_numeric_value is not None and max_numeric_value is not None and min_numeric_value > max_numeric_value:
        raise HTTPException(status_code = 400, detail = ("why is your minimum greater than your maximum?"))
    
    analysis_data = build_clinical_genomic_analysis_table(
        study_id = study_id,
        gene_symbols = gene_symbols,
        molecular_profile_id=molecular_profile_id,
        sample_list_id=sample_list_id
    )

    successful_genes = analysis_data["successful_genes"]
    if normalized_mutation_gene and normalized_mutation_gene not in successful_genes:
        raise HTTPException( status_code=400, detail={ "message": 
            (f"Mutation data were not successfully retrieved for gene '{normalized_mutation_gene}'." ), 
            "successful_genes": successful_genes, "suggestion": ( "Include the mutation gene in the genes parameter and confirm it exists in the selected study." ) 
            } )
    
    filtered_rows = []
    for row in analysis_data["analysis_table"]:
        include_row = True
        if normalized_mutation_gene:
            mutation_column = f"{normalized_mutation_gene}_mutated"
            required_mutation_value = mutation_status == "mutated"
            if row.get(mutation_column) != required_mutation_value:
                include_row = False
        
        if include_row and clinical_attribute: 
            row_value = row.get(clinical_attribute)
            if clinical_value is not None: 
                if row_value is None:
                    include_row = False
                else: 
                    actual_value = str(row_value).strip().casefold()
                    required_value = clinical_value.strip().casefold()
                    if actual_value != required_value:
                        include_row = False
            if include_row and (min_numeric_value is not None or max_numeric_value is not None):
                    try: numeric_value = float(row_value) 
                    except (TypeError, ValueError):
                        include_row = False
                    else:
                        if min_numeric_value is not None and numeric_value < min_numeric_value:
                            include_row = False
                        if max_numeric_value is not None and numeric_value > max_numeric_value:
                            include_row = False
        if include_row: 
            filtered_rows.append(row)
    return{
            "molecular_profile_id": analysis_data["molecular_profile_id"],
            "sample_list_id": analysis_data["sample_list_id"],
            "requested_genes": analysis_data["requested_genes"],
            "successful_genes": successful_genes,
            "filters": {
                "mutation_gene": normalized_mutation_gene,
                "mutation_status": mutation_status,
                "clinical_attribute": clinical_attribute,
                "clinical_value": clinical_value,
                "min_numeric_value": min_numeric_value,
                "max_numeric_value": max_numeric_value},
            "unfiltered_row_count": len(analysis_data["analysis_table"]),
            "filtered_row_count": len(filtered_rows),
            "cohort": filtered_rows,
            "errors": analysis_data["errors"]
        }
                    
                        




@app.get(
    "/studies/{study_id}/genes/{gene_symbol}/mutations",
    tags=["Mutations"],
    summary="Get clean mutation records for one gene in one study"
)
def get_gene_mutations(
    study_id: str,
    gene_symbol: str,
    molecular_profile_id: Optional[str] = Query(
        None,
        description="Optional mutation molecular profile ID. Leave blank to detect one automatically"
    ),
    sample_list_id: Optional[str] = Query(
        None,
        description="Optional sample list ID. Leave blank to detect one automatically."
    ),
    page_size: int = Query(100, ge=1, le=10000),
    page_number: int = Query(0, ge=0)
):
    mutation_data = fetch_gene_mutations(
        study_id=study_id,
        gene_symbol=gene_symbol,
        molecular_profile_id=molecular_profile_id,
        sample_list_id=sample_list_id,
        page_size=page_size,
        page_number=page_number

    )

    return {
        "source": "cBioPortal",
        "study_id": study_id,
        "gene": mutation_data["gene"],
        "molecular_profile_id": mutation_data["molecular_profile_id"],
        "sample_list_id": mutation_data["sample_list_id"],
        "page_size": page_size,
        "page_number": page_number,
        "count_on_this_page": len(mutation_data["mutations"]),
        "mutations": mutation_data["mutations"]
    }



@app.get(
    "/studies/{study_id}/genes/{gene_symbol}/mutation-summary",
    tags=["Mutations"],
    summary="Summarize mutations for one gene in one study"
)
def get_gene_mutation_summary(
    study_id: str,
    gene_symbol: str,
    molecular_profile_id: Optional[str] = Query(
        None,
        description="Optional mutation molecular profile ID. Leave blank to detect one automatically."
    ),
    sample_list_id: Optional[str] = Query(
        None,
        description="Optional sample list ID. Leave blank to detect one automatically."
    )
):
    selected_profile_id = (
        molecular_profile_id
        if molecular_profile_id
        else find_mutation_profile_id(study_id)
    )

    selected_sample_list_id = (
        sample_list_id
        if sample_list_id
        else find_mutation_sample_list_id(study_id)
    )

    profiled_sample_ids = cbioportal_get(
        f"/api/sample-lists/{selected_sample_list_id}/sample-ids"
    )

    summary = build_gene_mutation_summary(
        study_id=study_id,
        gene_symbol=gene_symbol.upper(),
        molecular_profile_id=selected_profile_id,
        sample_list_id=selected_sample_list_id,
        profiled_sample_ids=profiled_sample_ids
    )

    return {
        "source": "cBioPortal",
        "study_id": study_id,
        "molecular_profile_id": selected_profile_id,
        "sample_list_id": selected_sample_list_id,
        **summary
    }

@app.get(
    "/studies/{study_id}/mutation-panel/summary",
    tags=["Mutation Panels"],
    summary="Compare mutation summaries for multiple genes in one study"
)
def get_mutation_panel_summary(
    study_id: str,
    genes: str = Query(
        ...,
        description="Comma-separated gene symbols, such as BRAF,NRAS,NF1,KIT"
    ),
    molecular_profile_id: Optional[str] = Query(
        None,
        description="Optional mutation molecular profile ID. Leave blank to detect one automatically."
    ),
    sample_list_id: Optional[str] = Query(
        None,
        description="Optional sample list ID. Leave blank to detect one automatically."
    )
):
    gene_symbols = parse_gene_symbols(genes)

    panel_data = build_mutation_panel_summary(
        study_id=study_id,
        gene_symbols=gene_symbols,
        molecular_profile_id=molecular_profile_id,
        sample_list_id=sample_list_id
    )

    return {
        "source": "cBioPortal",
        "study_id": study_id,
        "requested_gene_count": len(gene_symbols),
        "successful_gene_count": len(panel_data["summaries"]),
        "error_count": len(panel_data["errors"]),
        "molecular_profile_id": panel_data["molecular_profile_id"],
        "sample_list_id": panel_data["sample_list_id"],
        "profiled_sample_count": panel_data["profiled_sample_count"],
        "gene_summaries": panel_data["summaries"],
        "errors": panel_data["errors"]
    }


@app.get(
    "/studies/{study_id}/mutation-panel/export.csv",
    tags=["Mutation Panels"],
    summary="Export a multi-gene mutation summary as CSV"
)
def export_mutation_panel_csv(
    study_id: str,
    genes: str = Query(
        ...,
        description="Comma-separated gene symbols, such as BRAF,NRAS,NF1,KIT"
    ),
    molecular_profile_id: Optional[str] = Query(
        None,
        description="Optional mutation molecular profile ID. Leave blank to detect one automatically."
    ),
    sample_list_id: Optional[str] = Query(
        None,
        description="Optional sample list ID. Leave blank to detect one automatically."
    )
):
    gene_symbols = parse_gene_symbols(genes)

    panel_data = build_mutation_panel_summary(
        study_id=study_id,
        gene_symbols=gene_symbols,
        molecular_profile_id=molecular_profile_id,
        sample_list_id=sample_list_id
    )

    csv_rows = []

    for summary in panel_data["summaries"]:
        gene = summary.get("gene") or {}
        mutation_type_counts = summary.get("mutation_type_counts") or {}

        mutation_type_text = "; ".join(
            f"{mutation_type}: {count}"
            for mutation_type, count in mutation_type_counts.items()
        )

        csv_rows.append({
            "study_id": study_id,
            "gene_symbol": gene.get("hugo_gene_symbol"),
            "entrez_gene_id": gene.get("entrez_gene_id"),
            "molecular_profile_id": panel_data["molecular_profile_id"],
            "sample_list_id": panel_data["sample_list_id"],
            "profiled_sample_count": summary.get("profiled_sample_count"),
            "mutated_sample_count": summary.get("mutated_sample_count"),
            "mutation_record_count": summary.get("mutation_record_count"),
            "mutation_frequency": summary.get("mutation_frequency"),
            "mutation_frequency_percent": summary.get("mutation_frequency_percent"),
            "mutation_type_counts": mutation_type_text
        })

    filename = f"{study_id}_mutation_panel_summary.csv"

    return rows_to_csv_response(csv_rows, filename)


@app.get(
    "/studies/{study_id}/mutation-panel/matrix",
    tags=["Mutation Panels"],
    summary="Build a sample-by-gene mutation matrix"
)
def get_sample_gene_mutation_matrix(
    study_id: str,
    genes: str = Query(
        ...,
        description="Comma-separated gene symbols, such as BRAF, NRAS, NF1, KIT"
    ),
    molecular_profile_id: Optional[str] = Query(
        None,
        description="Optional mutation molecular profile ID. Leave blank to detect automatically"
    ),
    sample_list_id: Optional[str] = Query(
        None,
        description="Optional sample list ID. Leave blank to detect automatically"
    ),
    page_size: int = Query(
        100,
        ge=1,
        le=1000,
        description="Number of sample rows to return per page"
    ),
    page_number: int = Query(
        0,
        ge=0,
        description="Page number starting from 0"
    )
):
    gene_symbols = parse_gene_symbols(genes)

    matrix_data = build_sample_gene_mutation_matrix(
        study_id=study_id,
        gene_symbols=gene_symbols,
        molecular_profile_id=molecular_profile_id,
        sample_list_id=sample_list_id
    )

    all_rows = matrix_data["matrix"]

    start_index = page_number * page_size
    end_index = start_index + page_size

    page_rows = all_rows[start_index:end_index]
    total_rows = len(all_rows)

    total_pages = (
        (total_rows + page_size - 1) // page_size
        if total_rows > 0
        else 0
    )

    return {
        "source": "cBioPortal",
        "study_id": study_id,
        "molecular_profile_id": matrix_data[
            "molecular_profile_id"
        ],
        "sample_list_id": matrix_data["sample_list_id"],
        "requested_genes": matrix_data["requested_genes"],
        "successful_genes": matrix_data["successful_genes"],
        "profiled_sample_count": matrix_data[
            "profiled_sample_count"
        ],
        "page_size": page_size,
        "page_number": page_number,
        "total_pages": total_pages,
        "count_on_this_page": len(page_rows),
        "matrix": page_rows,
        "errors": matrix_data["errors"]
    }



@app.get(
    "/studies/{study_id}/mutation-panel/matrix/export.csv",
    tags=["Mutation Panels"],
    summary = "Export a sample-by-gene mutation matrix as CSV"
)
def export_sample_gene_mutation_matrix_csv(
    study_id: str,
    genes: str = Query(
        ...,
        description = "Comma-separated gene symbols, such as BRAF, NRAS, NF1, KIT"
    ),
    molecular_profile_id: Optional[str] = Query(
        None,
        description = "Optional mutation molecular profile ID. Leave blank to automatically detect"
    ),
    sample_list_id: Optional[str] = Query(
        None,
        description = "Optional sample list ID. Leave blank to automatically detect"
    )
):
    gene_symbols = parse_gene_symbols(genes)
    matrix_data = build_sample_gene_mutation_matrix(
        study_id = study_id,
        gene_symbols = gene_symbols,
        molecular_profile_id=molecular_profile_id,
        sample_list_id = sample_list_id
    )

    filename = f"{study_id}_sample_gene_mutation_matrix.csv"
    return rows_to_csv_response(
        matrix_data["matrix"],
        filename
    )

@app.get(
    "/studies/{study_id}/analysis-table",
    tags = ["Analysis Tables"],
    summary = "Combine clinical data and mutation status by sample"
)
def get_clinical_genomic_analysis_table(
    study_id: str,
    genes: str = Query(
        ..., description =("comma-separated gene symbols such as BRAF, NRAS, NF1, KIT")
    ),
    molecular_profile_id: Optional[str]=Query(
        None,
        description = "Optional mutation profile ID (leave blank to auto detect)"
    ),
    sample_list_id: Optional[str] = Query(None, description = "Optional sample list ID (leave blank to auto detect) "),
    page_size: int = Query(100, ge=1, le=1000, description = "Number of analysis rows per page"),
    page_number: int = Query(0, ge = 0, description = "Page number starting from 0")
):
    gene_symbols = parse_gene_symbols(genes)

    analysis_data = build_clinical_genomic_analysis_table(
        study_id = study_id,
        gene_symbols=gene_symbols,
        molecular_profile_id=molecular_profile_id,
        sample_list_id = sample_list_id,
    )
    all_rows = analysis_data["analysis_table"]
    start_index=page_number*page_size
    end_index = start_index+page_size

    page_rows = all_rows[start_index: end_index]
    total_rows = len(all_rows)

    total_pages= (
        (total_rows+ page_size - 1) // page_size
        if total_rows > 0 
        else 0 
    )
    return {
        "source": "cBioPortal",
        "study_id": study_id,
        "molecular_profile_id": analysis_data["molecular_profile_id"],
        "sample_list_id": analysis_data["sample_list_id"],
        "requested_genes": analysis_data["requested_genes"],
        "successful_genes": analysis_data["successful_genes"],
        "profiled_sample_count": analysis_data["profiled_sample_count"],
        "clinical_sample_count": analysis_data["clinical_sample_count"],
        "clinical_patient_count": analysis_data["clinical_patient_count"],
        "samples_without_patient_clinical_data_count": len(analysis_data["samples_without_patient_clinical_data"]),
        "analysis_row_count": analysis_data["analysis_row_count"],
        "samples_without_clinical_data_count": len(analysis_data["samples_without_clinical_data"]),
        "page_size" : page_size,
        "page_number": page_number,
        "total_pages": total_pages,
        "count_on_this_page": len(page_rows),
        "analysis_table": page_rows,
        "errors": analysis_data["errors"]
    }

@app.get("/studies/{study_id}/analysis-table/export.csv", tags = ["Analysis Tables"], summary = "Export a combined clinical-genomic anaylsis table")
def export_clinical_genomic_analysis_table_csv(
    study_id: str,
    genes: str = Query(..., description= ("Comma-separated gene symbols such as BRAF, NRAS, NF1, KIT")),
    molecular_profile_id: Optional[str] = Query(None, description = ("Optional mutation molecular profile ID (leave blank to auto detect)")),
    sample_list_id: Optional[str] = Query(None, description = ("Optional sample list ID. Leave blank to auto detect")),

):
    gene_symbols = parse_gene_symbols(genes)
    analysis_data = build_clinical_genomic_analysis_table(
        study_id = study_id,
        gene_symbols = gene_symbols,
        molecular_profile_id=molecular_profile_id,
        sample_list_id=sample_list_id,

    )
    filename = (f"{study_id}_clinical_genomic_analysis_table.csv")
    return rows_to_csv_response(analysis_data["analysis_table"], filename)

@app.get("/studies/{study_id}/cohort", tags = ["Cohort Builder"], summary = "Filter a clinical-genomic cohort")
def get_filtered_clinical_genomic_cohort(
    study_id: str,
    genes: str = Query(..., description = "Comma-separated genes used to build the mutation matrix (ie: BRAF, NRAS, Nf1, KIT)"),
    mutation_gene: Optional[str] = Query(None, description = ("Optional gene to filter, such as BRAF. (gene must also appear in genes)")),
    mutation_status: Optional[Literal["mutated","not_mutated"]
                              ] = Query(None, description = "Return samples that are mutated or not mutated for mutation_gene."),
    clinical_attribute: Optional[str] = Query( None, description = "Optional clinical column such as OS_status, age, sex, stage"),
    clinical_value: Optional[str] = Query(None, description = "Optional exact clinical value. (case insensitive)"),
    min_numeric_value: Optional[float] = Query(None, description = "Optional minimum numeric value for selected clinical attribute"),
    max_numeric_value: Optional[float] = Query(None, description = "Optional maximum numeric value for selected clinical attribute"),
    molecular_profile_id: Optional[str] = Query(None, description = "Optional mutation molecular profile ID. leave blank to autodetect"),
    sample_list_id: Optional[str] = Query(None, description = "Optional sample_list_ID. leave blank to autodetect"),
    page_size: int = Query(100, ge = 1, le=1000, description = "Number of cohorts per page"),
    page_number: int = Query(0, ge= 0, description="page number 0")

):
    gene_symbols = parse_gene_symbols(genes)
    cohort_data = build_filtered_clinical_genomic_cohort(
        study_id= study_id,
        gene_symbols=gene_symbols,
        molecular_profile_id=molecular_profile_id,
        sample_list_id=sample_list_id,
        mutation_gene= mutation_gene,
        mutation_status=mutation_status,
        clinical_attribute = clinical_attribute,
        clinical_value = clinical_value,
        min_numeric_value=min_numeric_value,
        max_numeric_value=max_numeric_value
    )
    all_rows = cohort_data["cohort"]

    start_index = page_number * page_size
    end_index = start_index + page_size
    
    page_rows = all_rows[start_index:end_index]
    total_rows  = len(all_rows)
    total_pages= ((total_rows + page_size - 1)//page_size
        if total_rows > 0
        else 0
    )
    return {
        "source": "cBioPortal",
        "study_id": study_id,
        "molecular_profile_id": cohort_data["molecular_profile_id"],
        "sample_list_id": cohort_data["sample_list_id"],
        "requested_genes": cohort_data["requested_genes"],
        "successful_genes": cohort_data["successful_genes"],
        "filters": cohort_data["filters"],
        "unfiltered_row_count": cohort_data["unfiltered_row_count"],
        "filtered_row_count": cohort_data["filtered_row_count"],
        "page_size": page_size,
        "page_number": page_number,
        "total_pages": total_pages,
        "count_on_this_page": len(page_rows),
        "cohort": page_rows,
        "errors": cohort_data["errors"]
    }
                       
@app.get("/studies/{study_id}/cohort/export.csv", tags = ["Cohort Builder"], summary= "Export a filtered clinical-genomic cohort")
def export_filtered_clinical_genomic_cohort_csv(
    study_id: str,
    genes: str = Query(..., description = "CSV gene symbols"),
    mutation_gene: Optional[str] = Query(None),
    mutation_status: Optional[ Literal["mutated","not_mutated"]]= Query(None),
    clinical_attribute : Optional[str] = Query(None),
    clinical_value: Optional[str] = Query(None),
    min_numeric_value: Optional[float] = Query(None),
    max_numeric_value: Optional[float] = Query(None),
    molecular_profile_id: Optional[str] = Query(None),
    sample_list_id: Optional[str] = Query(None),
):
    gene_symbols = parse_gene_symbols(genes)
    cohort_data = build_filtered_clinical_genomic_cohort(
        study_id = study_id,
        gene_symbols = gene_symbols,
        molecular_profile_id=molecular_profile_id,
        sample_list_id = sample_list_id,
        mutation_gene = mutation_gene,
        mutation_status = mutation_status,
        clinical_attribute= clinical_attribute,
        clinical_value = clinical_value,
        min_numeric_value=min_numeric_value,
        max_numeric_value = max_numeric_value

    )
    filename = f"{study_id}_filtered_cohort.csv"
    return rows_to_csv_response(cohort_data["cohort"], filename)
    
