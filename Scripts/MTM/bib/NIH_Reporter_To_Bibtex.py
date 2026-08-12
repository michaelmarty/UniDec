import csv
import html
import os
import re


def _strip_html_tags(text: str) -> str:
    if not text:
        return ""
    text = re.sub(r"<[^>]+>", "", text)
    return html.unescape(text).strip()


def _extract_href(text: str) -> str:
    if not text:
        return ""
    m = re.search(r"href=['\"]([^'\"]+)['\"]", text, flags=re.IGNORECASE)
    return m.group(1).strip() if m else ""


def _extract_pmcid(pmc_value: str) -> str:
    v = (pmc_value or "").strip()
    if not v:
        return ""
    return v if v.upper().startswith("PMC") else f"PMC{v}"


def _parse_pub_date(pub_date: str):
    parts = (pub_date or "").strip().split()
    if len(parts) >= 2:
        month = parts[1]
        day = parts[2] if len(parts) >= 3 else ""
        return month, day
    return "", ""


def _sanitize_key(text: str) -> str:
    text = re.sub(r"[^A-Za-z0-9]+", "", text or "")
    return text.lower()[:40] if text else "entry"


def _bib_escape(value: str) -> str:
    if value is None:
        return ""
    return (
        str(value)
        .replace("\\", "\\\\")
        .replace("{", "\\{")
        .replace("}", "\\}")
        .replace("&", "\\&")
    )


def _format_field(name: str, value: str) -> str:
    if not value:
        return ""
    return f"  {name:<12} = {{{_bib_escape(value)}}},\n"


def row_to_bibtex_entry(row: dict) -> str:
    core_project = row.get("Core Project Number", "").strip()
    affiliation = row.get("Affilation", "").strip()
    authors_raw = row.get("Authors", "").strip()
    country = row.get("Country", "").strip()
    issn = row.get("ISSN", "").strip()
    issue = row.get("Journal Issue", "").strip()
    journal_html = row.get("Journal (Link to PubMed abstract)", "").strip()
    journal_abbr = row.get("Journal Title ABBR", "").strip()
    volume = row.get("Journal Volume", "").strip()
    language = row.get("Language", "").strip()
    pages = row.get("Page Number", "").strip()
    pmc_id_raw = row.get("PMC ID", "").strip()
    pmid = row.get("PMID", "").strip()
    pub_date = row.get("PUB Date", "").strip()
    title_html = row.get("Title (Link to full-text in PubMed Central)", "").strip()
    pub_year = row.get("PUB Year", "").strip()
    try:
        icite_rcr = row.get("iCite RCR", "").strip()
    except:
        icite_rcr = ""

    journal = _strip_html_tags(journal_html)
    title = _strip_html_tags(title_html)
    pubmed_url = _extract_href(journal_html)
    pmc_url = _extract_href(title_html)

    authors = " and ".join([a.strip() for a in authors_raw.split(";") if a.strip()])
    pages = pages.replace("-", "--")
    month, day = _parse_pub_date(pub_date)
    pmcid = _extract_pmcid(pmc_id_raw)

    first_author_last = authors_raw.split(";")[0].split(",")[0].strip() if authors_raw else "entry"
    key = f"{_sanitize_key(first_author_last)}{pub_year or 'year'}{pmid or _sanitize_key(title)}"

    entry = f"@article{{{key},\n"
    entry += _format_field("author", authors)
    entry += _format_field("title", title)
    entry += _format_field("journal", journal)
    entry += _format_field("journaltitle", journal_abbr)
    entry += _format_field("year", pub_year)
    entry += _format_field("month", month)
    entry += _format_field("day", day)
    entry += _format_field("volume", volume)
    entry += _format_field("number", issue)
    entry += _format_field("pages", pages)
    entry += _format_field("issn", issn)
    entry += _format_field("language", language)
    entry += _format_field("pmid", pmid)
    entry += _format_field("pmcid", pmcid)
    entry += _format_field("url", pubmed_url)
    entry += _format_field("howpublished", pmc_url)
    entry += _format_field("affiliation", affiliation)
    entry += _format_field("country", country)
    entry += _format_field("core_project", core_project)
    if icite_rcr:
        entry += _format_field("note", f"iCite RCR: {icite_rcr}")

    entry = re.sub(r",\n\Z", "\n", entry)
    entry += "}\n"
    return entry


def csv_to_bibtex(csv_path: str, bib_path: str) -> None:
    entries = []
    with open(csv_path, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if not any((v or "").strip() for v in row.values()):
                continue
            entries.append(row_to_bibtex_entry(row))

    with open(bib_path, "w", encoding="utf-8", newline="\n") as out:
        out.write("% Auto-generated BibTeX from CSV\n\n")
        out.write("\n".join(entries))

    print(f"Wrote {len(entries)} entries to: {os.path.abspath(bib_path)}")


if __name__ == "__main__":
    csv_file = r"C:\Users\mm96978\OneDrive - The University of Texas at Austin\Grants\RM1 Renewal\Publications_11193837_20May2026_153311.csv"
    bib_file = r"C:\Users\mm96978\OneDrive - The University of Texas at Austin\Grants\RM1 Renewal\Publications_11193837_20May2026_153311.bib"
    csv_to_bibtex(csv_file, bib_file)
