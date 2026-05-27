import os
import re
import pandas as pd
import xml.etree.ElementTree as ET
import numpy as np
from PublicScripts.Lipids.LipidFunctions import map_ontology_from_lipid_name, apply_sum_comp


def parse_semicolon_array(text):
    """Parse a semicolon-separated string into a numpy float array."""
    if text is None or text.strip() == "":
        return np.array([])
    return np.array([float(x) for x in text.strip().split(";")])


def parse_range_string(text):
    """Parse a range string like '[0.95;1.05]' into (min, max) tuple."""
    if text is None:
        return None, None
    text = text.strip().replace("&#x5d;", "]").replace("[", "").replace("]", "")
    parts = text.split(";")
    if len(parts) == 2:
        return float(parts[0]), float(parts[1])
    return None, None


def parse_spectral_db_entry(entry_el):
    """Parse a <spectraldatabaseentry> element."""
    result = {
        "library_file": entry_el.get("library_file"),
        "mzs": parse_semicolon_array(getattr(entry_el.find("mzs"), "text", None)),
        "intensities": parse_semicolon_array(getattr(entry_el.find("intensities"), "text", None)),
        "fields": {}
    }
    db_fields = entry_el.find("databasefieldslist")
    if db_fields is not None:
        for entry in db_fields.findall("entry"):
            result["fields"][entry.get("name")] = entry.text
    return result


def parse_spectrum_el(el):
    """Parse a spectrum element with <mzs> and <intensities> children."""
    if el is None:
        return None
    return {
        "mzs": parse_semicolon_array(getattr(el.find("mzs"), "text", None)),
        "intensities": parse_semicolon_array(getattr(el.find("intensities"), "text", None)),
    }


def parse_spectral_similarity(sim_el):
    """Parse a <spectralsimilarity> element."""
    if sim_el is None:
        return None
    result = {
        "similarity_function": getattr(sim_el.find("similairtyfunction"), "text", None),
        "overlapping_peaks": int(sim_el.find("overlappingpeaks").text) if sim_el.find("overlappingpeaks") is not None else None,
        "score": float(sim_el.find("score").text) if sim_el.find("score") is not None else None,
        "explained_library_intensity": float(sim_el.find("explainedLibraryIntensity").text) if sim_el.find("explainedLibraryIntensity") is not None else None,
        "library_spectrum": parse_spectrum_el(sim_el.find("libraryspectrum")),
        "query_spectrum": parse_spectrum_el(sim_el.find("queryspectrum")),
        "aligned_spectra": []
    }
    aligned_list = sim_el.find("alignedspectrumlist")
    if aligned_list is not None:
        for aligned in aligned_list.findall("alignedspectrum"):
            result["aligned_spectra"].append(parse_spectrum_el(aligned))
    return result


def parse_feature_annotation(ann_el):
    """Parse a <feature_annotation> element."""
    result = {
        "annotation_type": ann_el.get("annotation_type"),
        "version": ann_el.get("version"),
        "spectral_db_entry": None,
        "spectral_similarity": None,
        "mz_diff": None,
        "rt_absolute_error": None,
        "scan": None,
    }
    db_entry_el = ann_el.find("spectraldatabaseentry")
    if db_entry_el is not None:
        result["spectral_db_entry"] = parse_spectral_db_entry(db_entry_el)
    result["spectral_similarity"] = parse_spectral_similarity(ann_el.find("spectralsimilarity"))
    for dt in ann_el.findall("datatype"):
        dt_type = dt.get("type")
        if dt_type == "mz_diff":
            result["mz_diff"] = float(dt.text) if dt.text else None
        elif dt_type == "rt_absolute_error":
            result["rt_absolute_error"] = float(dt.text) if dt.text else None
    scan_el = ann_el.find("scan")
    if scan_el is not None:
        scan = -1
        try:
            scan = int(scan_el.get("scanindex"))
        except:
            pass
        result["scan"] = {
            "scantype": scan_el.get("scantype"),
            "rawdatafile": scan_el.get("rawdatafile"),
            "scanindex": scan,
        }
    return result


def parse_simple_ion_time_series(el):
    """Parse a <simpleiontimeseries> element."""
    if el is None:
        return None
    return {
        "scans": [int(x) for x in getattr(el.find("scans"), "text", "").split(";") if x],
        "intensities": parse_semicolon_array(getattr(el.find("intensities"), "text", None)),
        "mzs": parse_semicolon_array(getattr(el.find("mzs"), "text", None)),
        "rts": parse_semicolon_array(getattr(el.find("rts"), "text", None)),
    }


def parse_isotope_pattern(el):
    """Parse a <simpleisotopepattern> element."""
    if el is None:
        return None
    return {
        "mzs": parse_semicolon_array(getattr(el.find("mzs"), "text", None)),
        "intensities": parse_semicolon_array(getattr(el.find("intensities"), "text", None)),
        "description": getattr(el.find("description"), "text", None),
        "status": getattr(el.find("status"), "text", None),
        "charge": int(el.find("charge").text) if el.find("charge") is not None else None,
    }


def parse_feature(feature_el):
    """Parse a <feature> element."""
    result = {
        "rawdatafile": feature_el.get("rawdatafile"),
        "rt_range": None,
        "fragment_scans": [],
        "ion_time_series": None,
        "isotopes": None,
    }
    scalar_fields = [
        "rt", "mz", "height", "area", "charge", "fwhm",
        "tailing_factor", "asymmetry_factor", "feature_state",
        "rt_ms2_apex_distance", "datafile",
    ]
    for f in scalar_fields:
        result[f] = None

    for dt in feature_el.findall("datatype"):
        dt_type = dt.get("type")
        if dt_type == "rt_range":
            lo, hi = parse_range_string(dt.text)
            result["rt_range"] = (lo, hi)
        elif dt_type in ("rt", "mz", "height", "area", "fwhm",
                         "tailing_factor", "asymmetry_factor", "rt_ms2_apex_distance"):
            result[dt_type] = float(dt.text) if dt.text else None
        elif dt_type == "charge":
            result["charge"] = int(dt.text) if dt.text else None
        elif dt_type in ("feature_state", "datafile"):
            result[dt_type] = dt.text
        elif dt_type == "fragment_scans":
            for scan_el in dt.findall("scan"):
                result["fragment_scans"].append({
                    "scantype": scan_el.get("scantype"),
                    "rawdatafile": scan_el.get("rawdatafile"),
                    "scanindex": int(scan_el.get("scanindex")),
                })
        elif dt_type == "feature_data":
            fd = dt.find("feature_data")
            if fd is not None:
                its = fd.find("simpleiontimeseries")
                result["ion_time_series"] = parse_simple_ion_time_series(its)
        elif dt_type == "isotopes":
            iso_el = dt.find("simpleisotopepattern")
            result["isotopes"] = parse_isotope_pattern(iso_el)
        elif dt_type == "best_ms1_scan_number":
            scan_el = dt.find("scan")
            if scan_el is not None:
                result["best_ms1_scan"] = {
                    "scantype": scan_el.get("scantype"),
                    "rawdatafile": scan_el.get("rawdatafile"),
                    "scanindex": int(scan_el.get("scanindex")),
                }
        elif dt_type in ("mz_range", "intensity_range"):
            lo, hi = parse_range_string(dt.text)
            result[dt_type] = (lo, hi)
    return result


def parse_row(row_el):
    """Parse a <row> element into a dictionary."""
    row = {"id": int(row_el.get("id"))}
    scalar_types = {
        "height": float, "area": float, "mz": float, "rt": float,
        "charge": int, "feature_group": int,
    }
    bool_types = {"area_box_plot", "height_box_plot", "feature_shape"}
    list_types = {"fragment_scans", "spectral_db_matches"}

    row["spectral_db_matches"] = []
    row["fragment_scans"] = []
    row["features"] = []

    for dt in row_el.findall("datatype"):
        dt_type = dt.get("type")
        if dt_type in scalar_types:
            try:
                row[dt_type] = scalar_types[dt_type](dt.text)
            except (TypeError, ValueError):
                row[dt_type] = None
        elif dt_type in bool_types:
            row[dt_type] = dt.text.strip().lower() == "true" if dt.text else None
        elif dt_type == "rt_range":
            lo, hi = parse_range_string(dt.text)
            row["rt_range"] = (lo, hi)
        elif dt_type == "mz_range":
            lo, hi = parse_range_string(dt.text)
            row["mz_range"] = (lo, hi)
        elif dt_type == "intensity_range":
            lo, hi = parse_range_string(dt.text)
            row["intensity_range"] = (lo, hi)
        elif dt_type == "comment":
            row["comment"] = dt.text
        elif dt_type == "spectral_db_matches":
            for ann_el in dt.findall("feature_annotation"):
                row["spectral_db_matches"].append(parse_feature_annotation(ann_el))
        elif dt_type == "fragment_scans":
            for scan_el in dt.findall("scan"):
                row["fragment_scans"].append({
                    "scantype": scan_el.get("scantype"),
                    "rawdatafile": scan_el.get("rawdatafile"),
                    "scanindex": int(scan_el.get("scanindex")),
                })
        else:
            # Store remaining as raw text or None
            if dt.text and dt.text.strip() not in ("", "NULL_VALUE"):
                row.setdefault("extra", {})[dt_type] = dt.text

    for feature_el in row_el.findall("feature"):
        row["features"].append(parse_feature(feature_el))

    return row


def parse_mzmine_xml(filepath):
    """
    Parse an MZmine feature list XML file.

    Parameters
    ----------
    filepath : str
        Path to the XML file.

    Returns
    -------
    dict with keys:
        - 'metadata'  : dict with featurelistname, numberofrows, date
        - 'rows'      : list of row dicts, each containing all datatypes,
                        spectral_db_matches, fragment_scans, and features
    """
    tree = ET.parse(filepath)
    root = tree.getroot()

    metadata = {
        "featurelistname": root.get("featurelistname"),
        "numberofrows": int(root.get("numberofrows")),
        "date": root.get("date"),
    }

    rows = [parse_row(row_el) for row_el in root.findall("row")]

    return {"metadata": metadata, "rows": rows}


def _format_spectrum(mzs, intensities):
    """Format paired mz/intensity arrays as 'mz:intensity mz:intensity ...' string."""
    if mzs is None or len(mzs) == 0:
        return ""
    return " ".join(f"{m:.5f}:{i:.1f}" for m, i in zip(mzs, intensities))




def rows_to_dataframe(rows):
    """
    Flatten the list of parsed rows into a pandas DataFrame.
    Scalar fields only; complex nested structures (spectra, features) are dropped.

    Parameters
    ----------
    rows : list of dict
        Output from parse_mzmine_xml()['rows'].

    Returns
    -------
    pd.DataFrame
    """
    scalar_keys = ["id", "height", "area", "mz", "rt", "charge",
                   "feature_group", "comment", "area_box_plot",
                   "height_box_plot", "feature_shape"]
    records = []
    for row in rows:
        rec = {k: row.get(k) for k in scalar_keys}
        # Add best spectral match info if present
        if row.get("spectral_db_matches"):
            match = row["spectral_db_matches"][0]
            fields = match.get("spectral_db_entry", {}).get("fields", {})
            rec["match_name"] = fields.get("NAME")
            rec["match_formula"] = fields.get("FORMULA")
            rec["match_inchikey"] = fields.get("INCHIKEY")
            rec["match_precursor_mz"] = fields.get("PRECURSOR_MZ")
            rec["match_ion_type"] = fields.get("ION_TYPE")
            rec["match_score"] = match.get("spectral_similarity", {}).get("score") if match.get("spectral_similarity") else None
            rec["match_rt"] = fields.get("RT")
            rec["mz_diff"] = match.get("mz_diff")
            rec["rt_absolute_error"] = match.get("rt_absolute_error")
        records.append(rec)
    return pd.DataFrame(records)


def rows_to_msdial_dataframe(rows):
    """
    Convert parsed MZmine XML rows into a DataFrame matching MS-DIAL column format.

    Translated columns
    ------------------
    Alignment ID                  <- row id
    Average Rt(min)               <- row rt
    Average Mz                    <- row mz
    Metabolite name               <- spectral_db_entry fields NAME
    Adduct type                   <- spectral_db_entry fields ION_TYPE
    MS/MS assigned                <- True if any spectral_db_matches present
    Reference RT                  <- spectral_db_entry fields RT
    Reference m/z                 <- spectral_db_entry fields PRECURSOR_MZ
    Formula                       <- spectral_db_entry fields FORMULA
    INCHIKEY                      <- spectral_db_entry fields INCHIKEY
    SMILES                        <- spectral_db_entry fields SMILES
    RT matched                    <- rt_absolute_error is not None
    m/z matched                   <- mz_diff is not None
    MS/MS matched                 <- True if spectral_db_matches present
    Comment                       <- spectral_db_entry fields COMMENT (or row comment)
    Weighted dot product          <- spectral_similarity score
    Matched peaks count           <- spectral_similarity overlapping_peaks
    Matched peaks percentage      <- overlapping_peaks / num library peaks * 100
    Spectrum reference file name  <- spectral_db_entry library_file
    MS1 isotopic spectrum         <- feature isotope pattern (mz:intensity pairs)
    MS/MS spectrum                <- spectral_similarity query_spectrum (mz:intensity pairs)
    Ref Spec                      <- spectral_db_entry mzs/intensities (mz:intensity pairs)
    S/N average                   <- row height / intensity_range[0] (approximation)
    Fill %                        <- detected features / total feature elements * 100

    Not available (left blank)
    --------------------------
    Post curation result, Annotation tag (VS1.0),
    Manually modified for quantification, Manually modified for annotation,
    Isotope tracking parent ID, Isotope tracking weight number,
    RT similarity, m/z similarity, Simple dot product, Reverse dot product, Total score

    Parameters
    ----------
    rows : list of dict
        Output from parse_mzmine_xml()['rows'].

    Returns
    -------
    pd.DataFrame
    """
    records = []
    for row in rows:
        rec = {}

        # --- Core feature scalars ---
        rec["Alignment ID"] = row.get("id")
        rec["Average Rt(min)"] = row.get("rt")
        rec["Average Mz"] = row.get("mz")

        # --- Blanks (not available from MZmine XML) ---
        for blank_col in (
             "Ontology", "Annotation tag (VS1.0)",
            "RT similarity", "m/z similarity",
            "Simple dot product", "Reverse dot product", "Total score",
        ):
            rec[blank_col] = ""

        # --- Spectral match fields ---
        matches = row.get("spectral_db_matches", [])
        has_match = len(matches) > 0
        rec["MS/MS assigned"] = has_match
        rec["MS/MS matched"] = has_match

        if has_match:
            match = matches[0]
            db_entry = match.get("spectral_db_entry") or {}
            fields = db_entry.get("fields", {})
            sim = match.get("spectral_similarity") or {}

            rec["Metabolite name"] = fields.get("NAME", "")
            rec["Ontology"] = map_ontology_from_lipid_name(rec["Metabolite name"])
            rec["Adduct type"] = fields.get("ION_TYPE", "")
            rec["Reference RT"] = fields.get("RT", "")
            rec["Reference m/z"] = fields.get("PRECURSOR_MZ", "")
            rec["Formula"] = fields.get("FORMULA", "")
            rec["INCHIKEY"] = fields.get("INCHIKEY", "")
            rec["SMILES"] = fields.get("SMILES", "")
            rec["Comment"] = fields.get("COMMENT", row.get("comment", ""))

            mz_diff = match.get("mz_diff")
            rt_err = match.get("rt_absolute_error")
            rec["m/z matched"] = mz_diff is not None
            rec["RT matched"] = rt_err is not None

            # Similarity / scoring
            rec["Weighted dot product"] = sim.get("score", "")

            n_overlap = sim.get("overlapping_peaks")
            rec["Matched peaks count"] = n_overlap if n_overlap is not None else ""

            lib_mzs = db_entry.get("mzs")
            n_lib = len(lib_mzs) if lib_mzs is not None else 0
            if n_overlap is not None and n_lib > 0:
                rec["Matched peaks percentage"] = round(n_overlap / n_lib * 100, 2)
            else:
                rec["Matched peaks percentage"] = ""

            rec["Spectrum reference file name"] = db_entry.get("library_file", "")

            # Spectra formatted as "mz:intensity ..." strings
            rec["Ref Spec"] = _format_spectrum(
                db_entry.get("mzs"), db_entry.get("intensities")
            )
            qs = sim.get("query_spectrum") or {}
            rec["MS/MS spectrum"] = _format_spectrum(
                qs.get("mzs"), qs.get("intensities")
            )
        else:
            for col in ("Metabolite name", "Adduct type", "Reference RT", "Reference m/z",
                        "Formula", "INCHIKEY", "SMILES", "Comment",
                        "m/z matched", "RT matched", "Weighted dot product",
                        "Matched peaks count", "Matched peaks percentage",
                        "Spectrum reference file name", "Ref Spec", "MS/MS spectrum"):
                rec[col] = ""
            rec["Comment"] = row.get("comment", "")

        # --- MS1 isotopic spectrum from first feature ---
        ms1_iso = ""
        features = row.get("features", [])
        if features:
            iso = features[0].get("isotopes") or {}
            # ms1_iso = _format_spectrum(iso.get("mzs"), iso.get("intensities"))
            ms1_iso = ""
        rec["MS1 isotopic spectrum"] = ms1_iso

        # --- Fill % (detected features / total feature slots) ---
        n_features = len(features)
        n_detected = sum(
            1 for f in features if f.get("feature_state") == "DETECTED"
        )
        rec["Fill %"] = round(n_detected / n_features * 100, 1) if n_features else ""

        # --- S/N approximation: height / min_intensity ---
        height = row.get("height")
        intensity_range = row.get("intensity_range")
        if height and intensity_range and intensity_range[0]:
            rec["S/N average"] = round(height / intensity_range[0], 2)
        else:
            rec["S/N average"] = ""

        # --- Extra pass-through columns not in MS-DIAL format ---
        rec["area"] = row.get("area")
        rec["charge"] = row.get("charge")
        rec["feature_group"] = row.get("feature_group")
        rec["rt_range_min"] = row.get("rt_range", (None, None))[0]
        rec["rt_range_max"] = row.get("rt_range", (None, None))[1]
        rec["mz_range_min"] = row.get("mz_range", (None, None))[0]
        rec["mz_range_max"] = row.get("mz_range", (None, None))[1]
        rec["intensity_range_min"] = row.get("intensity_range", (None, None))[0]
        rec["intensity_range_max"] = row.get("intensity_range", (None, None))[1]
        rec["n_fragment_scans"] = len(row.get("fragment_scans", []))
        rec["n_spectral_matches"] = len(matches)
        if has_match:
            rec["mz_diff"] = matches[0].get("mz_diff", "")
            rec["rt_absolute_error"] = matches[0].get("rt_absolute_error", "")
            rec["explained_library_intensity"] = (
                (matches[0].get("spectral_similarity") or {}).get("explained_library_intensity", "")
            )

        records.append(rec)

    # Enforce MS-DIAL column order first, extras appended after
    msdial_cols = [
        "Alignment ID", "Average Rt(min)", "Average Mz", "Metabolite name",
        "Adduct type", "Fill %", "MS/MS assigned",
        "Reference RT", "Reference m/z", "Formula", "Ontology", "INCHIKEY",
        "SMILES", "Annotation tag (VS1.0)", "RT matched", "m/z matched",
        "MS/MS matched", "Comment",
        "RT similarity", "m/z similarity", "Simple dot product",
        "Weighted dot product", "Reverse dot product",
        "Matched peaks count", "Matched peaks percentage", "Total score",
        "S/N average", "Spectrum reference file name",
        "MS1 isotopic spectrum", "MS/MS spectrum", "Ref Spec",
    ]
    extra_cols = [c for c in records[0].keys() if c not in msdial_cols] if records else []
    df = pd.DataFrame(records, columns=msdial_cols + extra_cols)

    df = apply_sum_comp(df)

    return df

def run_all_xml_files(directory):
    """Process all .xml files in a directory and save as _msdial.csv."""
    for filename in os.listdir(directory):
        if filename.endswith(".xml"):
            filepath = os.path.join(directory, filename)
            print(f"\nProcessing {filepath}...")
            data = parse_mzmine_xml(filepath)
            df = rows_to_msdial_dataframe(data["rows"])
            print("Number of rows:", len(df))
            out_path = os.path.splitext(filepath)[0] + "_msdial.csv"
            df.to_csv(out_path, index=False)
            print("Unique ontologies:", df["Ontology"].nunique())
            # print(f"Saved {out_path}\n")


# ---------------------------------------------------------------------------
# Example usage
# ---------------------------------------------------------------------------
if __name__ == "__main__":

    # file = r"C:\Data\Lipidomics\test.xml"
    # file = r"Z:\Group Share\Annika\DDA_Compare\Stellar_Neg.xml"

    run_all_xml_files(r"Z:\Group Share\Annika\DDA_Compare\mzmine")
    exit()


    data = parse_mzmine_xml(file)
    print("Feature list:", data["metadata"]["featurelistname"])
    print("Number of rows:", data["metadata"]["numberofrows"])
    print("Parsed rows:", len(data["rows"]))

    df = rows_to_msdial_dataframe(data["rows"])
    print(df.head().to_string())
    print(df.columns.tolist())

    print(df["Ontology"].unique())
    # exit()


    out = file.replace(".xml", "_msdial.csv")
    df.to_csv(out, index=False)
    print("Saved to:", out)

