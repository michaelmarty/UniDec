# Engine state and file lifecycle

An engine holds mutable `config`, `data`, and `pks` objects. `config` owns run settings and filenames, `data` holds the mode's arrays or project spectra, and `pks` tracks picked peaks. The conventional engine and its derived modes share configuration history through `UniDecEngine`.

## Single-spectrum lifecycle

| Operation | In-memory and file effects |
| --- | --- |
| Open | Conventional `open_file` resets data and peaks, sets the managed paths, imports the source, and writes or reuses its raw cache. It can restore a processed spectrum and read the existing configuration. IM-MS starts with the newly loaded mobility grid and explicitly requires preprocessing. |
| Process | Conventional and IM-MS `process_data` use the current settings, write the processed native input, and set `config.procflag` on success. CD-MS rebuilds the event filters and histogram and replaces the previous fit and peak evidence. |
| Deconvolve | The conventional engine exports current settings and clears previous peaks, mass distribution, fit, and charge grids before invoking the executable. It imports outputs on exit status zero. An import can still fail if a required output is absent or malformed. |
| Pick peaks | Peak picking uses the imported mass distribution. Conventional peak picking can rescale its intensity column and export a peak list. Re-run the appropriate stages after changing inputs or settings to refresh downstream results. |

`procflag` indicates whether the conventional engine has prepared data. It can be set from a cache and does not track which settings produced the input. The direct engine's `run_unidec()` attempts processing if the flag is clear. A validated `config.update_parameters()` changes the configuration in memory. Explicitly process again when you change a preprocessing setting. Treat arrays as mutable working state and copy any results you need to retain across operations.

`update_history()`, `undo()`, and `redo()` operate on configuration snapshots. File contents and result arrays belong to the run lifecycle. ZIP state archives capture the working files, and restoring an archive can unpack files and reconstruct an input spectrum. See the [conventional workflow](../guides/conventional-ms.md) for the managed checkpoint route.

## Managed freshness

The single-spectrum session tracks a revision and stage completion for that revision. Successful open, parameter update, or checkpoint restore clears freshness. Running preprocessing invalidates all downstream stages. Running deconvolution requires fresh processing and clears peak availability. Peak picking requires a fresh deconvolution. The session rejects stale results even when arrays remain on the engine. A failed operation has a structured error and may have written working files before it failed. Inspect the returned state and stage availability before using results.

MetaUniDec and UniChrom manage collections instead of individual text work directories. MetaUniDec copies the source project and works in the copy. UniChrom creates a derived project and rebuilds its spectra when windows change. Each session runs its collection stages in order and requires a current peak-picking stage for result summaries. See the [MetaUniDec guide](../guides/metaunidec.md) and the [UniChrom guide](../guides/unichrom.md) for project-specific behavior.
