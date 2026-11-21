# Data and peaks

**Tier 1.** Read these containers to inspect results from an engine. The containers are mutable Python objects that can be inspected without wxPython or the native executable. The stages that populate them have their own [requirements](modes.md).

| Import path | Where to look |
| --- | --- |
| `unidec.modules.unidecstructure.DataContainer` | On a single-spectrum engine, `engine.data` holds arrays such as `rawdata`, `data2` (processed m/z spectrum), and `massdat` (deconvolved mass spectrum). IM-MS adds mobility arrays such as `rawdata3` and `data3`. |
| `unidec.modules.peakstructure.Peaks` and `unidec.modules.peakstructure.Peak` | After peak picking, `engine.pks.peaks` gives individual peaks with fields such as `mass` and `height`. |
| Session `get_state()` and `get_results(max_peaks=...)` | Use these to inspect freshness and retrieve bounded result dictionaries for single spectra, collections, and chromatograms. |

Arrays can be empty before their stage runs. Read the results after a successful run and copy arrays you need to retain independently of later engine operations. Reading a container or a session result summary does not write an export.
