# Engines

The following direct call paths are **Tier 1** for Python scripts and notebooks. Engines keep mutable data and results on the instance. Use a dedicated work location for runs because processing writes intermediate files and native deconvolution writes output that the engine imports. GUI dependencies are unnecessary for these routes. The native executable is required for deconvolution, and vendor formats can have additional runtime requirements.

| Mode | Import | Start the workflow |
| --- | --- | --- |
| Conventional MS | `from unidec.engine import UniDec` | Follow the [conventional-MS example](../guides/conventional-ms.md) to open, process, deconvolve, and pick peaks. |
| IM-MS | `from unidec.modules.IMEng import UniDecIM` | Open, preprocess, and deconvolve a mobility spectrum. See the [IM-MS workflow](../guides/im-ms.md). |
| CD-MS | `from unidec.modules.CDEng import UniDecCD` | Open with `open_cdms_file`. Call `run_deconvolution` to process and deconvolve by default, then pick and score peaks. See the [CD-MS workflow](../guides/cd-ms.md). |
| MetaUniDec | `from unidec.metaunidec.mudeng import MetaUniDec` | Open the HDF5 project with `open`. Call `process_data`, `run_unidec`, and `pick_peaks`. See the [MetaUniDec workflow](../guides/metaunidec.md). |
| UniChrom | `from unidec.modules.ChromEng import ChromEngine` | Open the raw source with `open_chrom`. See the [UniChrom workflow](../guides/unichrom.md) for managed construction and analysis of time windows. |

Inspect single-spectrum [arrays and peaks](data-and-peaks.md) after the corresponding stage. Use engine `save_state` and `load_state` for the conventional or IM-MS archive workflow. Use [sessions](sessions.md) when you need path validation, bounded results, or checkpoint policies.

## Supported Python reference

IM-MS uses the conventional methods for opening, running the executable, picking peaks, and saving an archive. CD-MS also inherits conventional peak picking. UniChrom inherits MetaUniDec's spectrum preprocessing method. See the mode guides for the input and output differences.

| Entry point | Supported operations |
| --- | --- |
| {any}`UniDec <unidec.engine.UniDec>` | {any}`open_file <unidec.engine.UniDec.open_file>`, {any}`pass_data_in <unidec.engine.UniDec.pass_data_in>`, {any}`process_data <unidec.engine.UniDec.process_data>`, {any}`run_unidec <unidec.engine.UniDec.run_unidec>`, {any}`pick_peaks <unidec.engine.UniDec.pick_peaks>`, {any}`save_state <unidec.engine.UniDec.save_state>`, {any}`load_state <unidec.engine.UniDec.load_state>` |
| {any}`UniDecIM <unidec.modules.IMEng.UniDecIM>` | {any}`process_data <unidec.modules.IMEng.UniDecIM.process_data>`, {any}`load_state <unidec.modules.IMEng.UniDecIM.load_state>` |
| {any}`UniDecCD <unidec.modules.CDEng.UniDecCD>` | {any}`open_cdms_file <unidec.modules.CDEng.UniDecCD.open_cdms_file>`, {any}`process_data <unidec.modules.CDEng.UniDecCD.process_data>`, {any}`run_deconvolution <unidec.modules.CDEng.UniDecCD.run_deconvolution>` |
| {any}`MetaUniDec <unidec.metaunidec.mudeng.MetaUniDec>` | {any}`open <unidec.metaunidec.mudeng.MetaUniDec.open>`, {any}`process_data <unidec.metaunidec.mudeng.MetaUniDec.process_data>`, {any}`run_unidec <unidec.metaunidec.mudeng.MetaUniDec.run_unidec>`, {any}`pick_peaks <unidec.metaunidec.mudeng.MetaUniDec.pick_peaks>` |
| {any}`ChromEngine <unidec.modules.ChromEng.ChromEngine>` | {any}`open_chrom <unidec.modules.ChromEng.ChromEngine.open_chrom>`, {any}`run_unidec <unidec.modules.ChromEng.ChromEngine.run_unidec>`, {any}`pick_peaks <unidec.modules.ChromEng.ChromEngine.pick_peaks>` |
