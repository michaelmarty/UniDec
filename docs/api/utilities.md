# Scientific utilities and importers

**Tier 2.** Use these areas when you need control over a scientific calculation or input reader. Check the signatures, array conventions, and optional dependencies for the installed version.

The [scientific analysis guide](../guides/scientific-analysis.md) gives runnable examples for preprocessing, scoring, fitting, matching, and mass-defect analysis.

| Area | Import or reference | Requirements and side effects |
| --- | --- | --- |
| Array preprocessing, integration, and peak detection | `unidec.tools`, including `datachop`, `integrate`, and `peakdetect` | These examples operate on in-memory arrays and run without wxPython or native UniDec. Check individual functions for array and configuration mutation. |
| IM-MS calculations | `unidec.modules.IM_functions` | Used by the IM-MS engine. Check calibration settings and units before calling mode-specific helpers directly. |
| Peak and fit scoring | `UniDec.dscore`, `UniDecCD.dscore`, and `unidec.modules.ChromEng.chromatographic_peak_scores` | Engine scoring updates peak objects and uses evidence from the current deconvolution. The chromatographic helper scores an extracted trace in memory. |
| Peak fitting | `unidec.modules.fitting.fit_peak` and `unidec.modules.MassFitter.MassFitter` | Fits a single peak or a mixture. Supply initial guesses and inspect fit parameters in their documented order. |
| Mass matching | `unidec.modules.matchtools`, including `site_match` | Analyze picked peaks against candidate masses. Spreadsheet exporters such as `sitematch_to_excel` write a file. |
| Mass defects | `unidec.tools.simple_mass_defect` and `unidec.tools.kendrick_analysis` | Computes an individual defect or a map in memory. The engine's `kendrick_continuous` method writes text exports. |
| Input readers | `UniDecImporter.get_importer` and the [importer guide](../guides/importing.md) | Selects a reader by extension and path type. The guide covers supported formats, vendor SDKs, array shapes, and scan or time selection. |

Start with the [API map](index.md) if you need an end-to-end workflow. Search the [internal module reference](internals/index.md) for implementation details and additional expert symbols.
