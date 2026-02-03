---
title: LSP Module
path: src/tunacode/lsp
type: directory
depth: 1
description: Language Server Protocol client for diagnostics
exports: [LSPClient, get_diagnostics, start_server]
seams: [M, D]
---

# LSP Module

## Purpose
Implements a minimal Language Server Protocol client to fetch diagnostic feedback (errors, warnings) from language servers.

## Key Components

### client.py
**LSPClient Class**
Main LSP client implementation:
- **start_server()** - Spawn language server process
- **initialize()** - Initialize LSP handshake
- **get_diagnostics()** - Fetch diagnostics for file
- **shutdown()** - Gracefully terminate server
- **did_open()** - Notify file open
- **did_change()** - Notify file change
- **did_close()** - Notify file close

**Features:**
- Async subprocess management
- JSON-RPC communication
- Server lifecycle management
- Error handling and recovery

### diagnostics.py
**Diagnostic Functions:**
- **get_diagnostics()** - Fetch diagnostics for file
- **format_diagnostics()** - Format for display
- **parse_diagnostics()** - Parse LSP response

**Diagnostic Types:**
- Error - Critical issues
- Warning - Potential problems
- Information - Suggestions
- Hint - Style recommendations

### servers.py
**Language Server Registry**
Supported language servers:
- **Python** - pylsp, python-lsp-server
- **JavaScript/TypeScript** - typescript-language-server
- **Go** - gopls
- **Rust** - rust-analyzer
- **etc.**

**Server Configuration:**
```python
{
  "language": "python",
  "command": "pylsp",
  "args": [],
  "workspace": Path
}
```

## LSP Workflow

```
1. Start server process
2. Initialize LSP handshake
3. Notify file open (did_open)
4. Get diagnostics
5. Notify file change (did_change)
6. Get updated diagnostics
7. Notify file close (did_close)
8. Shutdown server
```

## Integration Points

- **tools/write_file.py** - Trigger diagnostics after write
- **tools/update_file.py** - Trigger diagnostics after update
- **tools/decorators.py** - LSP integration in file_tool decorator
- **ui/renderers/tools/diagnostics.py** - Display diagnostics

## Usage

Tool decorator integration:
```python
@file_tool(writes=True)
async def write_file(filepath: str, content: str):
    # Write file
    ...
    # Fetch diagnostics
    diagnostics = await lsp_client.get_diagnostics(filepath)
    # Append to output
    return result + diagnostics
```

## Seams (M, D)

**Modification Points:**
- Add new language servers
- Customize LSP client behavior
- Extend diagnostic formatting
- Add LSP features (completion, hover, etc.)

**Extension Points:**
- Implement additional LSP capabilities
- Add server auto-detection
- Create custom diagnostic formatters
- Add LSP server management
