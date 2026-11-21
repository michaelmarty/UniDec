# Engines

The following direct call paths are **Tier 1** for Python scripts and notebooks. Engines keep mutable data and results on the instance. Use a dedicated work location for runs because processing writes intermediate files and native deconvolution writes output that the engine imports. GUI dependencies are unnecessary for these routes. The native executable is required for deconvolution, and vendor formats can have additional runtime requirements.

| Mode | Import | Start the workflow |
| --- | --- | --- |
| Conventional MS | `from unidec.engine import UniDec` | Call `open_file`, `process_data`, `run_unidec`, then `pick_peaks`. |
| IM-MS | `from unidec.modules.IMEng import UniDecIM` | Use the same method sequence with the IM-MS engine and a mobility spectrum. |
| CD-MS | `from unidec.modules.CDEng import UniDecCD` | Open with `open_cdms_file`. Call `run_deconvolution` to process and deconvolve by default, then `pick_peaks`. |
| MetaUniDec | `from unidec.metaunidec.mudeng import MetaUniDec` | Open the HDF5 project with `open`. Call `process_data`, `run_unidec`, and `pick_peaks` to work on the collection. |
| UniChrom | `from unidec.modules.ChromEng import ChromEngine` | Open the raw source with `open_chrom`. Use the [UniChrom session](sessions.md#unichrom) when you want managed construction and analysis of time windows. |

Inspect single-spectrum [arrays and peaks](data-and-peaks.md) after the corresponding stage. Use engine `save_state` and `load_state` for the conventional or IM-MS archive workflow. Use [sessions](sessions.md) when you need path validation, bounded results, or checkpoint policies.
