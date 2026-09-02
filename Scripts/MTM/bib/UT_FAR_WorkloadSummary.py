"""Build the publication data used for the UT FAR workload summary."""

import sys
from pathlib import Path

import pandas as pd
from pybliometrics import scopus
from pybliometrics.scopus.exception import Scopus401Error


AUTHOR_ID = "9-s2.0-53984653400"
SCOPUS_API_KEY = "07070f5eff2d979e7792ebcc2096daf6"
# SCOPUS_API_KEY = "2e5e53fcc1fd031f5399bc9c57860d0a"
JOURNAL_EXPORT_PATH = Path(__file__).with_name("UT_FAR_Journals_for_JCR.xlsx")
JIF_DATABASE_PATH = Path(__file__).with_name("JIF.xlsx")
WORKLOAD_EXPORT_PATH = Path(__file__).with_name("UT_FAR_WorkloadSummary.xlsx")


def get_publications(author_id: str, refresh: bool | int = False) -> pd.DataFrame:
    """Return every publication in an author's Scopus profile."""
    scopus.init(keys=[SCOPUS_API_KEY])
    scopus_author_id = author_id.rsplit("-", maxsplit=1)[-1]
    author = scopus.AuthorRetrieval(scopus_author_id)
    try:
        documents = author.get_documents(refresh=refresh) or []
    except Scopus401Error as error:
        raise RuntimeError(
            "Scopus rejected SCOPUS_API_KEY. Replace it with an active key from "
            "https://dev.elsevier.com/apikey/manage."
        ) from error

    columns = {
        "coverDate": "Publication Date",
        "title": "Publication Title",
        "publicationName": "Journal Title",
        "issn": "ISSN",
        "eIssn": "eISSN",
        "aggregationType": "Source Type",
        "doi": "DOI",
        "eid": "Scopus EID",
    }
    publications = pd.DataFrame(documents).reindex(columns=columns)
    publications = publications.rename(columns=columns)
    excluded_titles = publications["Publication Title"].str.contains(
        r"Author Correction|Erratum", case=False, na=False
    )
    publications = publications.loc[~excluded_titles]
    return publications.sort_values(
        by="Publication Date", ascending=False, na_position="last"
    ).reset_index(drop=True)


def export_journals_for_jcr(
    publications: pd.DataFrame,
    output_path: Path = JOURNAL_EXPORT_PATH,
    jif_database_path: Path = JIF_DATABASE_PATH,
) -> pd.DataFrame:
    """Export unique journals, retaining existing JIFs and adding known ones."""
    is_journal = publications["Source Type"].str.casefold().eq("journal")
    journals = publications.loc[
        is_journal, ["Journal Title", "ISSN", "eISSN"]
    ].drop_duplicates()
    journals = journals.sort_values("Journal Title").reset_index(drop=True)

    existing_jifs = _build_jif_lookup(output_path)
    database_jifs = _build_jif_lookup(jif_database_path)
    journals["2025 JIF"] = journals.apply(
        lambda row: _find_jif(row, existing_jifs, database_jifs), axis=1
    )

    temporary_path = output_path.with_name(f"{output_path.stem}.tmp.xlsx")
    journals.to_excel(temporary_path, index=False)
    temporary_path.replace(output_path)
    return journals


def export_workload_summary(
    publications: pd.DataFrame,
    journals: pd.DataFrame,
    output_path: Path = WORKLOAD_EXPORT_PATH,
) -> pd.DataFrame:
    """Export publication details enriched with the available journal JIFs."""
    workload = publications.copy()
    journal_jifs = _build_jif_lookup_from_dataframe(journals)
    workload["2025 JIF"] = workload.apply(
        lambda row: _find_jif(row, journal_jifs), axis=1
    )

    publication_dates = pd.to_datetime(workload["Publication Date"], errors="coerce")
    workload.insert(
        1, "Publication Year", publication_dates.dt.year.astype("Int64")
    )

    temporary_path = output_path.with_name(f"{output_path.stem}.tmp.xlsx")
    with pd.ExcelWriter(temporary_path, engine="openpyxl") as writer:
        workload.to_excel(writer, sheet_name="Publications", index=False)
        journals.to_excel(writer, sheet_name="Journals", index=False)
        for sheet in writer.book.worksheets:
            sheet.freeze_panes = "A2"
            sheet.auto_filter.ref = sheet.dimensions
            for column_cells in sheet.columns:
                values = (
                    str(cell.value) if cell.value is not None else ""
                    for cell in column_cells
                )
                width = min(max(map(len, values), default=0) + 2, 80)
                sheet.column_dimensions[column_cells[0].column_letter].width = width
    temporary_path.replace(output_path)
    return workload


def _normalize_issn(value) -> str | None:
    """Normalize an ISSN for matching while preserving the displayed value."""
    if pd.isna(value):
        return None
    normalized = str(value).replace("-", "").strip().upper()
    return normalized if normalized and normalized != "N/A" else None


def _build_jif_lookup(workbook_path: Path) -> tuple[dict, dict]:
    """Return JIF lookups keyed by normalized ISSN and journal title."""
    by_issn = {}
    by_title = {}
    if not workbook_path.exists():
        return by_issn, by_title

    dataframe = pd.read_excel(workbook_path, dtype={"ISSN": str, "eISSN": str})
    if "2025 JIF" not in dataframe.columns:
        return by_issn, by_title

    return _build_jif_lookup_from_dataframe(dataframe)


def _build_jif_lookup_from_dataframe(
    dataframe: pd.DataFrame,
) -> tuple[dict, dict]:
    """Return JIF lookups from a dataframe containing journal metadata."""
    by_issn = {}
    by_title = {}
    if "2025 JIF" not in dataframe.columns:
        return by_issn, by_title

    for _, row in dataframe.iterrows():
        jif = row["2025 JIF"]
        if pd.isna(jif):
            continue
        title = str(row.get("Journal Title", "")).strip().casefold()
        if title:
            by_title.setdefault(title, jif)
        for column in ("ISSN", "eISSN"):
            issn = _normalize_issn(row.get(column))
            if issn:
                by_issn.setdefault(issn, jif)
    return by_issn, by_title


def _find_jif(row: pd.Series, *lookups: tuple[dict, dict]):
    """Find a JIF, preferring earlier lookups and ISSN matches."""
    identifiers = [_normalize_issn(row[column]) for column in ("ISSN", "eISSN")]
    title = str(row["Journal Title"]).strip().casefold()
    for by_issn, by_title in lookups:
        for identifier in identifiers:
            if identifier in by_issn:
                return by_issn[identifier]
        if title in by_title:
            return by_title[title]
    return pd.NA


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    publications = get_publications(AUTHOR_ID)
    journals = export_journals_for_jcr(publications)
    workload = export_workload_summary(publications, journals)
    print(f"Found {len(publications)} publications for Scopus author {AUTHOR_ID}:\n")
    with pd.option_context("display.max_rows", None, "display.max_colwidth", None):
        print(publications.to_string(index=False))
    print(f"\nExported {len(journals)} unique journals to {JOURNAL_EXPORT_PATH}")
    matched_jifs = workload["2025 JIF"].notna().sum()
    print(
        f"Exported {len(workload)} publications to {WORKLOAD_EXPORT_PATH} "
        f"({matched_jifs} with a 2025 JIF)"
    )
