---
summary: |
  LSP diagnostics provide automatic feedback when files are written or updated.
  Language servers (ruff, tsserver, etc.) check code and report errors/warnings
  as part of the tool result. The agent receives this context and decides how
  to respond.
when_to_read: |
  - When you want to enable or configure LSP diagnostics
  - When adding support for a new language
  - When debugging why diagnostics are not appearing
  - When understanding the architecture of LSP integration
---

# LSP Diagnostics

TunaCode provides automatic LSP diagnostics feedback when files are written or updated.

## File Map

```
tools/
├── lsp/                          # LSP implementation (tools layer)
│   ├── __init__.py               # get_diagnostics(), format_diagnostics()
│   ├── client.py                 # LSPClient - JSON-RPC over stdio
│   └── servers.py                # Server command mapping by extension
├── write_file.py                 # Calls tools.lsp after writing
└── update_file.py                # Calls tools.lsp after updating
```

## How It Works

When `write_file` or `update_file` tools are called:
1. File operation completes
2. Language server spawns (ruff, tsserver, gopls, etc.)
3. Diagnostics fetched for the file
4. Results appended to tool output

## Configuration

Enable/disable in `~/.tunacode/config.yaml`:

```yaml
settings:
  lsp:
    enabled: true
    timeout: 5.0
```

## Requirements

Install language servers for your languages:

- **Python**: `pip install ruff`
- **TypeScript/JavaScript**: `npm install -g typescript-language-server`
- **Go**: `go install golang.org/x/tools/gopls@latest`
- **Rust**: `rustup component add rust-analyzer`

## Output Format

```
<file_diagnostics>
ACTION REQUIRED: 2 error(s) found - fix before continuing
Error (line 15): undefined name 'x'
Error (line 23): missing return statement
</file_diagnostics>
```

The agent receives this as contextual feedback and decides how to respond.

## Architecture

LSP lives in `tools/lsp/` as a tool-layer concern:
- `tools.lsp.get_diagnostics()` - fetch diagnostics
- `tools.lsp.format_diagnostics()` - format for output
- `tools.write_file` and `tools.update_file` call these automatically
