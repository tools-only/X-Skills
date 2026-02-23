# mcp-server-motherduck — Cross-Platform Path Configuration

Verified against v1.0.1 source code read from `motherduckdb/mcp-server-motherduck` via GitHub API (2026-02-22).

---

## The Core Problem

MCP client `mcpServers` JSON configs pass `args` as a string array. The array is passed to the subprocess directly — no shell interprets it. Therefore:

- `~` is NOT expanded to the home directory
- `$HOME` is NOT substituted
- `%USERPROFILE%` is NOT substituted

The `mcp-server-motherduck` source code also performs no expansion. Evidence: zero matches for `expanduser`, `expandvars`, `os.path.expand*` in `__init__.py` and `database.py` (grep verified 2026-02-22).

---

## Solution Patterns by Client

### Pattern 1: MCP_DB_PATH Environment Variable (Recommended)

The `--db-path` option reads from `MCP_DB_PATH` envvar (source: `__init__.py` line 51, `envvar="MCP_DB_PATH"`).

The MCP client host process inherits the OS environment. Setting `MCP_DB_PATH` in the system or user environment means the args array never needs a path.

```json
{
  "mcpServers": {
    "duckdb": {
      "command": "uvx",
      "args": ["mcp-server-motherduck", "--read-write"],
      "env": {
        "MCP_DB_PATH": "/Users/jane/data/my.duckdb"
      }
    }
  }
}
```

Each user/machine sets `MCP_DB_PATH` in their profile. The checked-in JSON stays static.

### Pattern 2: Absolute Path Hardcoded Per Machine

Acceptable when the config is per-machine (not shared). Use the absolute path directly in args:

```json
{
  "args": ["mcp-server-motherduck", "--db-path", "/Users/jane/data/my.duckdb", "--read-write"]
}
```

### Pattern 3: In-Memory (No Path Required)

When persistent storage is not needed:

```json
{
  "args": ["mcp-server-motherduck", "--db-path", ":memory:", "--read-write", "--allow-switch-databases"]
}
```

Note: `--read-write` is mandatory with `:memory:`. Omitting it causes a `UsageError` at startup.

### Pattern 4: MotherDuck (No Local Path)

When using MotherDuck cloud storage, no local path is needed:

```json
{
  "args": ["mcp-server-motherduck", "--db-path", "md:", "--read-write"],
  "env": {
    "motherduck_token": "<TOKEN>"
  }
}
```

---

## Windows-Specific: HOME and DuckDB Extensions

DuckDB uses the `HOME` environment variable to locate its extensions directory. On Windows, `HOME` is not set by default (Windows uses `USERPROFILE`).

The `--home-dir` flag sets `os.environ["HOME"]` for the server process (source: `database.py` lines 63-64):

```python
if home_dir:
    os.environ["HOME"] = home_dir
```

**Limitation**: `--home-dir` must be an absolute path in the args array (same problem as `--db-path`). Since it also cannot use tilde, the value must be literal.

Workaround options:

1. Set `HOME` in the system environment before launching the MCP client. DuckDB will find it without `--home-dir`.
2. Pass `--home-dir C:\Users\username` as a literal absolute path in the args array.
3. Set `HOME` in the `env` section of the MCP config:

```json
{
  "env": {
    "HOME": "C:\\Users\\username",
    "MCP_DB_PATH": "C:\\Users\\username\\data\\my.duckdb"
  }
}
```

---

## Environment Variable Precedence

For `--db-path` / `MCP_DB_PATH`:

1. CLI flag `--db-path` takes precedence over `MCP_DB_PATH` envvar (Click behavior: envvar is fallback when CLI flag not provided)
2. If neither is provided, default is `:memory:`

For MotherDuck token:

1. `--motherduck-token` CLI flag
2. `motherduck_token` environment variable
3. `MOTHERDUCK_TOKEN` environment variable

Source: `__init__.py` line 58: `envvar=["motherduck_token", "MOTHERDUCK_TOKEN"]` (Click checks in order).

---

## Client-Specific Notes

### Claude Code

Supports `env` section in `mcpServers` config. `MCP_DB_PATH` in `env` works reliably. Claude Code issue #1254 reports that env vars from `env` section were not always passed in earlier versions — verify on current Claude Code version.

Source: <https://github.com/anthropics/claude-code/issues/1254> (referenced in web search, 2026-02-22)

### Claude Desktop

JSON config at `~/Library/Application Support/Claude/claude_desktop_config.json` (macOS) or `%APPDATA%\Claude\claude_desktop_config.json` (Windows). The `env` section is supported. Set `MCP_DB_PATH` there.

### Cursor

Settings -> MCP -> Add new global MCP server. Uses JSON with `command`, `args`, `env` structure. `env` section is supported.

### OpenCode / Codex

Both support `env` section alongside `args`. `MCP_DB_PATH` in `env` is the correct approach.

---

## What Does NOT Work

| Approach | Why it fails |
|---|---|
| `"--db-path", "~/mydb.duckdb"` in args | Tilde is not expanded by subprocess or by server |
| `"--db-path", "$HOME/mydb.duckdb"` in args | Shell variable not interpolated in JSON args array |
| `"--db-path", "%USERPROFILE%\\mydb.duckdb"` in args | Windows variable not interpolated in JSON args array |
| `"--home-dir", "~"` in args | Same — no expansion |

---

## Sources

- `src/mcp_server_motherduck/__init__.py` read via `gh api repos/motherduckdb/mcp-server-motherduck/contents/...` (2026-02-22)
- `src/mcp_server_motherduck/database.py` read via `gh api` (2026-02-22)
- `README.md` read via `gh api` (2026-02-22)
- <https://github.com/motherduckdb/mcp-server-motherduck> (official repository, accessed 2026-02-22)
