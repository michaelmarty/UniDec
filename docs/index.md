# UniDec documentation

UniDec provides Python and native tools for deconvolving mass and ion-mobility mass spectra.

- **Install or run UniDec.** See the [project README](https://github.com/michaelmarty/UniDec/blob/master/readme.md) for installation, supported file types, GUI launch commands, and release information.
- **Use UniDec from Python.** Start with the [Python API map](api/index.md) to choose an engine or session for your data.
- **Automate UniDec with an agent.** Follow the [UniDec MCP agent interface](guides/agent-interface.md).
- **Explore the implementation.** Start with the [architecture overview](architecture/overview.md) or browse the [internal module reference](api/internals/index.md).

```{toctree}
:caption: Guides
:maxdepth: 2

Agent interface <guides/agent-interface>
Conventional MS from Python <guides/conventional-ms>
Ion mobility MS from Python <guides/im-ms>
Charge detection MS from Python <guides/cd-ms>
MetaUniDec collections from Python <guides/metaunidec>
UniChrom time windows from Python <guides/unichrom>
Scientific analysis from Python <guides/scientific-analysis>
Import mass-spectrometry data <guides/importing>
Extend the UniDec GUI <guides/gui-extensions>
```

```{toctree}
:caption: Architecture
:maxdepth: 2

architecture/overview
```

```{toctree}
:caption: API reference
:maxdepth: 2

Python API map <api/index>
Internal module reference <api/internals/index>
```
