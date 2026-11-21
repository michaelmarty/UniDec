# Scientific utilities and importers

**Tier 2.** Use these areas when you need control over a scientific calculation or input reader. Check the signatures, array conventions, and optional dependencies for the installed version.

| Area | Import or reference | Requirements and side effects |
| --- | --- | --- |
| Array preprocessing, integration, and peak detection | `unidec.tools`, including `datachop`, `integrate`, and `peakdetect` | These examples operate on in-memory arrays and run without wxPython or native UniDec. Check individual functions for array and configuration mutation. |
| IM-MS calculations | `unidec.modules.IM_functions` | Used by the IM-MS engine. Check calibration settings and units before calling mode-specific helpers directly. |
| Mass matching | `unidec.modules.matchtools`, including `site_match` | Analyze picked peaks against candidate masses. Spreadsheet exporters such as `sitematch_to_excel` write a file. |
| Input readers | `UniDecImporter.ImporterFactory.ImporterFactory` and the [file types in the project README](https://github.com/michaelmarty/UniDec#unidec-compatible-file-types) | Reader requirements vary by format and platform. Some vendor readers require external runtimes. Opening inputs through an engine can also create working files. |

Start with the [API map](index.md) if you need an end-to-end workflow. Search the [internal module reference](internals/index.md) for implementation details and additional expert symbols.
