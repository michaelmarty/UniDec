# Import mass-spectrometry data

UniDec delegates acquisition-file reading to UniDecImporter. Use `UniDecImporter.get_importer(path)` to inspect data before running an engine. The engines select the same readers through `UniDecImporter.ImporterFactory.ImporterFactory.create_importer`. Each reader advertises `imms_support`, `cdms_support`, and `chrom_support`. Check the relevant capability and the file's acquisition content. For example, an ordinary mzML file can be readable as LC-MS while containing no ion-mobility coordinates. The formats and requirements below describe the project's declared UniDecImporter dependency, beginning at version 1.0.5. [UniDecImporter documentation](https://pypi.org/project/UniDecImporter/1.0.5/)

## File formats

| Input | Reader and intended use |
| --- | --- |
| `.txt`, `.dat` | Whitespace-delimited numeric text, with an optional leading header. Use two columns for conventional MS, three for IM-MS, or the event columns below for CD-MS. |
| `.csv`, `.tsv`, `.xlsx` | Single-spectrum m/z and intensity from the first two columns. Excel reads the active sheet via openpyxl. Use whitespace text or numeric NPZ to retain extra mobility or event columns. |
| `.npz`, `.bin` | NPZ expects a numeric array under the key `data`. Binary input is interleaved float64 m/z and intensity pairs. |
| `.mzML`, indexed `.mzML.gz` | pymzML reads scan-based files and chromatograms. IM-MS requires an export with a mobility coordinate for each point, such as a combined-ion-mobility mzML. CD-MS requires suitable single-ion event data. |
| `.mzXML` | pyteomics reads scan-based spectra and chromatograms. CD-MS can read event-like values through the generic reader, which fills missing event metadata with defaults. |
| `.dmt`, `.i2ms` | SQLite-backed CD-MS event data. |
| `.h5`, `.hdf5` | A MetaUniDec project is opened through the collection engine or session. These project files have their own HDF5 structure. |

The open acquisition formats work on Windows, macOS, and Linux. A reader context manager calls `close()` to release files, database connections, or vendor handles. The factory raises `UnsupportedFormatError` for an unknown extension. Legacy object-backed NPZ files can trigger pickle deserialization. Use numeric arrays for new files and load legacy object archives only from trusted sources. [UniDecImporter format and resource documentation](https://pypi.org/project/UniDecImporter/1.0.5/)

## Vendor inputs

The vendor readers require 64-bit x86 Windows and the appropriate SDK. Open formats such as mzML can be used on other platforms. Convert with a tool that preserves the scan metadata and, where needed, the ion-mobility or single-ion information.

| Path | Vendor runtime |
| --- | --- |
| Thermo `.raw` **file** | Thermo RawFileReader assemblies and `pythonnet`. UniDecImporter includes the assemblies. For a standalone importer installation, the `thermo` extra installs the bridge. Set `THERMO_RAW_FILE_READER_DIR` to select a different assembly directory if needed. |
| Waters `.raw` **directory** | Licensed Waters MassLynx SDK. Set `MASSLYNX_RAW_DLL` to the local `MassLynxRaw.dll`. |
| Agilent `.d` **directory** | Licensed MassHunter Data Access SDK and `pythonnet`. Set `AGILENT_DA_SDK_DIR` to the local assembly directory. |

A missing SDK or unsupported platform can raise `VendorReaderUnavailableError`. Managed workflow methods surface importer failures as structured errors such as `open_failed`. `UniDecSession.open_spectrum` and the corresponding MCP tool expect a file, so use the direct engine or a converted file for vendor directories. `UniChromSession.open_chromatogram` accepts both files and supported vendor directories. [UniDecImporter vendor requirements](https://pypi.org/project/UniDecImporter/1.0.5/)

## Arrays by mode

| Mode | Reader or project contract |
| --- | --- |
| Conventional MS | `get_single_scan(scan)` and `get_avg_scan(...)` return `(N, 2)` rows of m/z and intensity. Supply finite, ascending m/z for single-spectrum text input. The engine exposes the imported array as `data.rawdata`. |
| IM-MS | `get_imms_scan(scan)` and `get_imms_avg_scan(...)` return `(N, 3)` rows of m/z, drift or arrival time, and intensity. Use milliseconds for the time column of text input configured for the [IM-MS workflow](im-ms.md). The engine densifies sparse input into `data.rawdata3`. |
| CD-MS | `get_cdms_data()` returns `(N, 5)` event rows with m/z, intensity, scan identifier, inverse injection-time metadata, and retention time. Some readers insert defaults or apply their own intensity corrections. For numeric single-scan exports, missing metadata columns default to scan `1`, weight `1`, and time `-1`. The CD-MS engine exposes positive-intensity rows as `darray`. |
| MetaUniDec | The [collection workflow](metaunidec.md) takes an HDF5 project with m/z and intensity rows under `/ms_dataset/<index>/raw_data`. Use the reader to select spectra before adding the resulting arrays to a project, for example with `engine.data.add_data(spectrum, name=...)`. |
| UniChrom | The [chromatogram workflow](unichrom.md) reads TIC and BPC as `(N, 2)` retention-time and intensity arrays. Each time window becomes a separate spectrum with `(M, 2)` m/z and intensity rows in the managed HDF5 project. |

## Select scans or retention times

- For a scan-aware reader, `scan_range=(first, last)` uses inclusive **scan identifiers** from `reader.scans`. Retention times are in minutes. `time_range=(start, end)` resolves to the nearest endpoint scans and takes precedence when both ranges are supplied to `get_avg_scan`. With neither range, it uses the full run. `sum_mode=True` returns summed intensities. Readers can align scans onto a common m/z axis before combining them, so inspect the returned axis.
- Single-spectrum readers return the one spectrum regardless of scan or time selection. For CD-MS, selection depends on the reader. For example, the DMT/I2MS reader provides `get_cdms_data_by_scans((first, last))`. The direct CD-MS engine opens events first and applies `CDScanStart` and `CDScanEnd` during processing when those settings are supplied.
- The direct conventional and IM-MS engines pass `open_file(..., time_range=...)` to their respective averaging readers. Use a fresh working directory or `refresh=True` when changing a selection on a previously opened source to refresh cached data. The single-spectrum session exposes a path-only open. If it needs a subset, import that subset first and write a mode-appropriate input file.
- UniChrom converts each requested time interval to an inclusive scan range before averaging. Neighboring time windows can share a boundary scan. The session's optional time bounds constrain window starts. See the [UniChrom guide](unichrom.md) for the full window contract.

Direct reader calls return arrays. Engine opens can create raw caches, configuration, and other working files. The direct spectrum engine commonly uses a sibling `_unidecfiles` directory. The MetaUniDec session copies its input project, while the UniChrom session creates a derived HDF5 project in the workspace and keeps reading the original acquisition data for new windows.

## Example: select a retention-time range

The example writes a three-scan mzXML to a temporary directory using only Python and NumPy. It verifies the chosen scan identifiers, the difference between averaging and summing, and the conventional engine's imported spectrum. Import alone requires no vendor runtime, GUI, or native UniDec executable.

```python
from base64 import b64encode
from pathlib import Path
from tempfile import TemporaryDirectory
from xml.etree import ElementTree as ET

import numpy as np
from UniDecImporter import get_importer

from unidec.engine import UniDec

mz = np.array([500.0, 501.0, 502.0, 503.0])
root = ET.Element("mzXML")
run = ET.SubElement(root, "msRun", scanCount="3")
for scan_id, scale in enumerate((1.0, 2.0, 4.0), start=1):
    scan = ET.SubElement(
        run, "scan", num=str(scan_id), msLevel="1", peaksCount=str(len(mz)),
        retentionTime=f"PT{60 * (scan_id - 1)}S", polarity="+",
    )
    peaks = ET.SubElement(
        scan, "peaks", precision="32", byteOrder="network",
        pairOrder="m/z-int", compressionType="none",
    )
    pairs = np.column_stack((mz, scale * np.array([1.0, 2.0, 3.0, 1.0])))
    peaks.text = b64encode(pairs.astype(">f4").tobytes()).decode("ascii")

with TemporaryDirectory() as directory:
    source = Path(directory) / "scans.mzXML"
    ET.ElementTree(root).write(source, encoding="utf-8", xml_declaration=True)
    original = source.read_bytes()
    time_range = (0.9, 2.0)

    with get_importer(source) as reader:
        np.testing.assert_allclose(reader.times, [0.0, 1.0, 2.0])
        assert reader.get_scans_from_times(time_range) == [2, 3]
        np.testing.assert_allclose(reader.get_single_scan(2)[:, 1], [2, 4, 6, 2])
        average = reader.get_avg_scan(scan_range=(1, 1), time_range=time_range)
        by_scan = reader.get_avg_scan(scan_range=(2, 3))
        summed = reader.get_avg_scan(scan_range=(2, 3), sum_mode=True)
        np.testing.assert_allclose(average, by_scan)
        np.testing.assert_allclose(summed[:, 0], average[:, 0])
        np.testing.assert_allclose(summed[:, 1], 2 * average[:, 1])

    engine = UniDec(ignore_args=True, silent=True)
    engine.open_file(str(source), time_range=time_range, refresh=True, silent=True)
    np.testing.assert_allclose(engine.data.rawdata, average)
    assert source.read_bytes() == original
```
