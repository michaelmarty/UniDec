# Configuration

**Tier 1.** Python callers can inspect configuration through the engine's `config` attribute. Import `UniDecConfig` and `ParameterSpec` from `unidec.modules.unidecstructure`. This metadata and validation layer is headless and needs no native executable for schema lookup or validation.

- `ParameterSpec` describes registered settings, including their types, bounds or choices, units when applicable, and invalidation stages.
- `UniDecConfig.parameter_schema()` returns the public settings for the active mode. Read it from the selected engine's configuration so that it reflects the workflow's mode.
- `UniDecConfig.update_parameters(values)` validates the proposed mapping and updates the in-memory configuration. Use session `get_parameters` and `set_parameters` when you want the session to track revisions and invalidate downstream results.
- Engine and workflow operations handle persistence. Applying a configuration update alone does not run processing or produce new results.

The public schema is an allowlist for automation. The internal module reference also covers serialization and settings used by the GUI and native engine.

See the [Python-to-native boundary](../architecture/native-boundary.md) for text keys, HDF5 attributes, and native dispatch.

For conventional MS, `parameter_values()` returns a snapshot of the public values. Pass a name or sequence of names to select a subset. Unset optional values are represented as `None`. `validate_parameters(values)` checks and normalizes a proposed mapping without changing the config. `update_parameters(values)` applies the normalized mapping only after validation succeeds and returns the values it applied. Invalid names, types, values, or ranges raise `ValueError`. When updating paired limits, pass both in the same mapping if one limit on its own would make the range invalid.

Conventional MS requires `startz <= endz` and `masslb < massub`. When both m/z limits are set, `minmz < maxmz`. Native execution also checks run-specific settings and input requirements.

This example checks that validation alone leaves the current settings intact, then applies the updates together.

```{code-block} python
:name: validate-parameters

from unidec.modules.unidecstructure import UniDecConfig

config = UniDecConfig()
initial = config.parameter_values(("startz", "endz"))
selection = {"startz": 8, "endz": 24}
assert config.validate_parameters(selection) == selection
assert config.parameter_values(("startz", "endz")) == initial
assert config.update_parameters(selection) == selection
assert config.parameter_values(("startz", "endz")) == selection
assert config.numz == 17
```

```{toctree}
:maxdepth: 1

generated/parameters
```

## Supported Python reference

| Entry point | Supported operations |
| --- | --- |
| {any}`ParameterSpec <unidec.modules.unidecstructure.ParameterSpec>` | Defines parameter metadata for the schema. |
| {any}`UniDecConfig <unidec.modules.unidecstructure.UniDecConfig>` | {any}`parameter_schema <unidec.modules.unidecstructure.UniDecConfig.parameter_schema>`, {any}`parameter_values <unidec.modules.unidecstructure.UniDecConfig.parameter_values>`, {any}`validate_parameters <unidec.modules.unidecstructure.UniDecConfig.validate_parameters>`, {any}`update_parameters <unidec.modules.unidecstructure.UniDecConfig.update_parameters>` |
