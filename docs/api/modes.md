# Choose a mode

These are **Tier 1** automation routes. Select the mode when you construct the engine or start the MCP server. Each workflow uses native UniDec for deconvolution and can run without wxPython. Some input formats require a vendor importer and runtime.

See [importing data](../guides/importing.md) for supported formats, expected array columns, selection rules, and vendor requirements.

| Input and mode | Direct engine | Managed route | Main file effects |
| --- | --- | --- | --- |
| Conventional m/z and intensity spectrum (`ms`) | `unidec.engine.UniDec` | `unidec.session.UniDecSession` | The direct engine commonly writes configuration, input, and results in a per-spectrum work directory. The session uses its managed work directory. |
| [Ion mobility spectrum with drift time](../guides/im-ms.md) (`imms`) | `unidec.modules.IMEng.UniDecIM` | `UniDecSession` with a `UniDecIM` engine | Processing and deconvolution produce mobility-specific files and grids along with the mass spectrum. |
| [Charge detection ion events](../guides/cd-ms.md) (`cdms`) | `unidec.modules.CDEng.UniDecCD` | `UniDecSession` with `UniDecCD` and a `UniDecCDSessionEngineAdapter` | Processing creates histogram-derived data. The session writes analysis files in its work directory. |
| [HDF5 project containing multiple spectra](../guides/metaunidec.md) (`metaunidec`) | `unidec.metaunidec.mudeng.MetaUniDec` | `unidec.meta_session.MetaUniDecCollectionSession` | The direct engine works on the open HDF5 project. The session copies the project before analysis. |
| [Raw LC-MS chromatogram](../guides/unichrom.md) (`unichrom`) | `unidec.modules.ChromEng.ChromEngine` | `unidec.unichrom_session.UniChromSession` | The session creates a managed HDF5 project and replaces its collection when you regenerate time windows. |

For direct method names, see [engines](engine.md). For Python session setup and MCP mode selection, see [sessions](sessions.md).
