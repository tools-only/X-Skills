# Research - Issue #314: Layer 2 Lateral Coupling (tools, indexing, lsp)

**Date:** 2026-01-27
**Owner:** claude-agent
**Phase:** Research
**Git Commit:** 1de46c7bba8a7a322fc00b878eb918cf224bd9ae

## Goal

Research the lateral coupling violations between Layer 2 modules (`tools`, `indexing`, `lsp`) to understand:
1. Which specific imports violate the architecture
2. How the imported code is being used
3. What refactoring strategies from prior PRs apply

## Architecture Context

From `docs/architecture/layers_html.html`:

```
Application Layers: ui → core → tools | indexing | lsp → utils
```

The `|` syntax means `tools`, `indexing`, and `lsp` are **peers** at Layer 2. They should NOT import from each other - only downward to `utils` and foundation (`types`, `configuration`).

**Current violations (from architecture diagram):**
- `tools → indexing` (1 import)
- `tools → lsp` (2 imports)

---

## Findings

### Violation 1: tools → indexing (CodeIndex)

**Location:** `src/tunacode/tools/glob.py:12`

```python
from tunacode.indexing import CodeIndex
```

**Usage (4 call sites):**
- Line 131: `index = CodeIndex.get_instance()`
- Line 132: `index.build_index()`
- Line 185: `all_files = code_index.get_all_files()`
- Line 189: `abs_path = code_index.root_dir / file_path`

**Purpose:** Performance optimization. The glob tool has two strategies:
1. **Filesystem scan** - Direct `os.scandir()` traversal (slow)
2. **Index lookup** - Use pre-built `CodeIndex` (fast)

When directory is project root and conditions are met, glob uses the index for faster pattern matching.

**GitHub permalink:** https://github.com/alchemiststudiosDOTai/tunacode/blob/1de46c7bba8a7a322fc00b878eb918cf224bd9ae/src/tunacode/tools/glob.py#L12

---

### Violation 2: tools → lsp (diagnostics)

**Location:** `src/tunacode/tools/decorators.py:73`

```python
from tunacode.lsp import format_diagnostics, get_diagnostics
```

**Usage (2 call sites):**
- Line 77: `get_diagnostics(Path(filepath), timeout=timeout)`
- Line 80: `return format_diagnostics(diagnostics)`

**Purpose:** Type-checking feedback for file modifications. When tools decorated with `@file_tool(writes=True)` write to files, the decorator fetches LSP diagnostics to provide immediate feedback to the LLM about type errors.

**Note:** Import is lazy (inside function, not module-level) to avoid circular imports and allow tools to function when LSP is unavailable.

**GitHub permalink:** https://github.com/alchemiststudiosDOTai/tunacode/blob/1de46c7bba8a7a322fc00b878eb918cf224bd9ae/src/tunacode/tools/decorators.py#L73

---

### Violation 3: tools → lsp (server status)

**Location:** `src/tunacode/tools/lsp_status.py:5`

```python
from tunacode.lsp.servers import get_server_command
```

**Usage (1 call site):**
- Line 34: `command = get_server_command(LSP_STATUS_CHECK_PATH)`

**Purpose:** UI status indicator. The resource bar needs to display which LSP server is active. `get_server_command()` detects available language servers by file extension.

**Data flow:**
1. UI's `ResourceBar` calls `core.lsp_status.get_lsp_status()`
2. Core facade delegates to `tools.lsp_status.get_lsp_status()`
3. Tools function calls `lsp.servers.get_server_command()`
4. Returns `(enabled: bool, server_name: str | None)`

**GitHub permalink:** https://github.com/alchemiststudiosDOTai/tunacode/blob/1de46c7bba8a7a322fc00b878eb918cf224bd9ae/src/tunacode/tools/lsp_status.py#L5

---

## Key Patterns / Solutions Found

Based on PR #317 (commit `076cbf3c`) which fixed issue #313:

### Strategy 1: Move to Foundation Layer (Preferred)

If the imported code is **read-only configuration or data**, move it to `configuration/` or `types/`.

**Applies to:** Violation 3 (`get_server_command`)
- `get_server_command()` is a pure function: file extension → command lookup
- `SERVER_CONFIG` is static configuration data
- Could move to `configuration/lsp_servers.py`

### Strategy 2: Protocol Abstraction

Define interface in `types/`, implement in each module, inject from core.

**Applies to:** Violation 1 (CodeIndex)
```python
# types/protocols.py
class FileRegistry(Protocol):
    def get_all_files(self) -> list[Path]: ...
    @property
    def root_dir(self) -> Path: ...
```

Then `glob()` accepts optional `registry: FileRegistry | None` injected by core.

### Strategy 3: Push Orchestration to Core

If the cross-layer call is **orchestration logic**, move it up to core.

**Applies to:** Violation 2 (LSP diagnostics)
- The decorator shouldn't orchestrate LSP calls - that's business logic
- Core should wrap file tools, check if write occurred, fetch diagnostics
- `@file_tool(writes=True)` loses `on_write` callback; core handles it

### Strategy 4: Merge Modules

If a tools module is only consumed by core, merge it into core.

**Applies to:** Violation 3 (lsp_status.py)
- `tools/lsp_status.py` is already wrapped by `core/lsp_status.py`
- Could merge the tools file into core and eliminate the indirection

---

## Recommended Refactoring Plan

| Violation | Strategy | Effort | Impact |
|-----------|----------|--------|--------|
| tools/glob.py → indexing.CodeIndex | Protocol abstraction | Medium | High (breaks coupling) |
| tools/decorators.py → lsp.get_diagnostics | Push to core | Medium | High (cleaner layers) |
| tools/lsp_status.py → lsp.servers | Move to configuration OR merge into core | Low | Medium |

### Order of Operations

1. **Violation 3 first** (lowest effort, self-contained)
   - Move `get_server_command()` and `SERVER_CONFIG` to `configuration/lsp_servers.py`
   - OR merge `tools/lsp_status.py` into `core/lsp_status.py`

2. **Violation 2 second** (medium effort)
   - Remove LSP orchestration from `@file_tool` decorator
   - Create hook point in core agent loop for post-write diagnostics

3. **Violation 1 last** (highest effort, architectural)
   - Define `FileRegistry` protocol in `types/protocols.py`
   - Have `glob()` accept optional registry parameter
   - Core injects `CodeIndex` when appropriate

---

## Knowledge Gaps

- **Testing:** Need to verify each module can be tested in isolation after refactoring
- **import-linter:** Need to update configuration after moves
- **Performance:** Need to verify CodeIndex optimization still works after protocol indirection

---

## References

- Issue: https://github.com/alchemiststudiosDOTai/tunacode/issues/314
- Architecture: `docs/architecture/layers_html.html`
- Prior fix PR #317: commit `076cbf3c`
- Prior fix PR #316: commit `68e16f12`
- Related research: `memory-bank/research/2026-01-27_11-53-37_issue-313-core-utils-layer-violation.md`

---

## Files Requiring Changes

| File | Change Type | Notes |
|------|-------------|-------|
| `src/tunacode/tools/glob.py` | Modify | Remove direct CodeIndex import, accept protocol |
| `src/tunacode/tools/decorators.py` | Modify | Remove LSP import, delegate to core |
| `src/tunacode/tools/lsp_status.py` | Delete or Merge | Move to core or configuration |
| `src/tunacode/lsp/servers.py` | May move | `get_server_command()` to configuration |
| `src/tunacode/core/agents/main.py` | Modify | Add post-write diagnostic hook |
| `src/tunacode/types/protocols.py` | Create | Add `FileRegistry` protocol |
