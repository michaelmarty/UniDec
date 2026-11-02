# UniDec MCP Agent Interface

The optional MCP adapter exposes one headless, conventional-MS UniDec session through a bounded set of tools. It does not provide shell access, arbitrary Python or filesystem operations, native command selection, GUI workflows, MetaUniDec, UniChrom, IM-MS, CD-MS, DoubleDec, or SMAS.

## Install and start

From `public/UniDec`, install UniDec and the optional MCP dependency into the Python environment that Codex will use:

```shell
python -m pip install -e ".[agent]"
```

The adapter currently targets the MCP 1.x `FastMCP` API, so the optional dependency is constrained to `mcp>=1.9,<2`. If MCP 2.x was installed previously, rerun the command above after updating the checkout, or downgrade it directly:

```shell
python -m pip install --upgrade "mcp>=1.9,<2"
```

Verify the selected API before registering the server:

```shell
python -c "from mcp.server.fastmcp import FastMCP; print('MCP 1.x FastMCP import OK')"
```

The server requires at least one readable input root and a separate working directory. Repeat `--input-root` when spectra or checkpoints live under multiple directories.

```shell
python -m unidec.mcp_server --input-root C:/Data/UniDec --working-directory C:/Data/UniDecAgent
```

The installed `unidec-mcp` entry point accepts the same arguments.

The server uses STDIO for MCP messages. Engine progress is sent to stderr so it does not corrupt the protocol stream. Do not use the command as an interactive console program.

## Configure Codex

Codex supports local STDIO MCP servers, and its CLI and IDE extension share MCP configuration. The default user configuration is `~/.codex/config.toml`; a trusted project can instead use `.codex/config.toml`. See the [official Codex MCP documentation](https://learn.chatgpt.com/docs/extend/mcp) for the current configuration and approval options.

Add the server from a terminal:

```shell
codex mcp add unidec -- python -m unidec.mcp_server --input-root C:/Data/UniDec --working-directory C:/Data/UniDecAgent
codex mcp list
```

Alternatively, configure it directly. Use paths valid on the host; forward slashes make Windows TOML examples easier to read.

```toml
[mcp_servers.unidec]
command = "python"
args = [
  "-m", "unidec.mcp_server",
  "--input-root", "C:/Data/UniDec",
  "--working-directory", "C:/Data/UniDecAgent",
]
cwd = "<path-to-UniDecDev>/public/UniDec"
startup_timeout_sec = 30
tool_timeout_sec = 600
default_tools_approval_mode = "prompt"
enabled = true
```

Restart the IDE extension after changing its MCP configuration. Keep approval mode at `prompt` until the local roots and expected tool behavior have been reviewed. A host timeout can stop waiting for a tool response, but it cannot safely interrupt a native UniDec calculation already in progress.

## Tool workflow

Use the tools in this order:

1. `open_spectrum(path)` opens a file below a configured input root.
2. `get_state()` and `get_parameters(names)` inspect the current revision and supported parameter schema.
3. `set_parameters(values)` validates and applies a parameter update atomically.
4. `run(stage)` executes `process`, `deconvolve`, `pick_peaks`, or `all`.
5. `get_results(max_peaks)` returns bounded structured results; `get_plot(name)` writes a managed PNG artifact.
6. `save_state(path)` writes a new ZIP below the working directory; `load_state(path)` restores a ZIP from an input root or the working directory.

Every mutation has an operation ID and invalidates stale downstream results where appropriate. Mutations are serialized and can return `busy`. Check `ok`, `error`, `state.revision`, and `state.available` before continuing. `save_state` refuses to overwrite an existing checkpoint.

## Security boundary

- Input paths are resolved and must remain below an explicitly configured root.
- Plots, operation files, and saved checkpoints remain below the working directory.
- Checkpoint archives with absolute or parent-traversal members are rejected before extraction.
- The adapter declares exactly nine domain tools; it has no generic command, import, or file tool.
- Run the server with a dedicated working directory and only the minimum input roots needed for the task.
