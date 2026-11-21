# Python API map

Choose an entry point by acquisition mode and by how much control you need over files and results. Use the [engines](engine.md) for direct access to NumPy results. Use [sessions](sessions.md) for managed output, validated settings, and bounded results. Use the [MCP interface](../guides/agent-interface.md) when the caller is an MCP client.

| Tier | Audience and stability | Entry points |
| --- | --- | --- |
| 1 · Supported automation | Preferred for scripts, notebooks, and agents. Support covers the workflows and operations described in the curated pages. | [Engines](engine.md), [configuration](configuration.md), [data and peaks](data-and-peaks.md), [sessions and MCP](sessions.md), [modes](modes.md) |
| 2 · Expert | For scientific analysis and advanced integration. Check the callable and its assumptions against the UniDec version you use. | [Scientific utilities and importers](utilities.md), [GUI extension points](../guides/gui-extensions.md) |
| 3 · Internal | For contributors investigating implementation details. Paths and behavior can change as the implementation evolves. | [Internal module reference](internals/index.md) |

The pages below give the current import paths, requirements, and file effects. For older IM-MS code that imports `UniDec` from `unidec.modules.IMEng`, use `UniDecIM` from that module. For conventional MS, use `UniDec` from `unidec.engine`. The curated names are defined in their respective modules.

```{toctree}
:maxdepth: 1

modes
engine
configuration
data-and-peaks
sessions
utilities
```
