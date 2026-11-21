# Configuration

**Tier 1.** Python callers can inspect configuration through the engine's `config` attribute. Import `UniDecConfig` and `ParameterSpec` from `unidec.modules.unidecstructure`. This metadata and validation layer is headless and needs no native executable for schema lookup or validation.

- `ParameterSpec` describes registered settings, including their types, bounds or choices, units when applicable, and invalidation stages.
- `UniDecConfig.parameter_schema()` returns the public settings for the active mode. Read it from the selected engine's configuration so that it reflects the workflow's mode.
- `UniDecConfig.update_parameters(values)` validates the proposed mapping and updates the in-memory configuration. Use session `get_parameters` and `set_parameters` when you want the session to track revisions and invalidate downstream results.
- Engine and workflow operations handle persistence. Applying a configuration update alone does not run processing or produce new results.

The public schema is an allowlist for automation. The internal module reference also covers serialization and settings used by the GUI and native engine.
