---
title: Tools Layer
summary: LLM-callable tool implementations and supporting subsystems (grep engine, LSP).
read_when: Adding a new tool, modifying tool error handling, or changing how tool schemas are generated.
depends_on: [types, infrastructure, configuration]
feeds_into: [core]
---

# Tools Layer

**Package:** `src/tunacode/tools/`

## What

Every capability the agent can invoke during a conversation. Each tool is an async function decorated with `@base_tool` or `@file_tool`, then converted to a tinyagent `AgentTool` via `to_tinyagent_tool()`.

## Key Files

### Tool Implementations

| File | Tool Name | Purpose |
|------|-----------|---------|
| `bash.py` | `bash` | Execute shell commands with timeout and output capture. |
| `glob.py` | `glob` | Find files by glob pattern, respecting ignore rules. |
| `grep.py` | `grep` | Regex search across files using ripgrep. |
| `read_file.py` | `read_file` | Read file contents with hash-tagged line numbers for validation. |
| `write_file.py` | `write_file` | Create or overwrite a file. |
| `hashline_edit.py` | `hashline_edit` | Apply validated edits using content-hash line references. |
| `list_dir.py` | `list_dir` | List directory contents with file metadata. |
| `web_fetch.py` | `web_fetch` | Fetch and convert web pages with URL security validation. |
| `discover.py` | `discover` | Find and map code related to concepts via natural language. |

### Framework

| File | Purpose |
|------|---------|
| `decorators.py` | `@base_tool` -- wraps tools with consistent error handling (`ToolRetryError` passthrough, catch-all to `ToolExecutionError`). `@file_tool` -- adds path-specific error mapping (`FileNotFoundError` to `ToolRetryError`, `PermissionError` to `FileOperationError`). `to_tinyagent_tool()` -- converts a decorated async function to an `AgentTool` with auto-generated OpenAI-function JSON schema. |
| `xml_helper.py` | Loads tool descriptions from XML prompt files. If a tool has a matching XML file, its docstring is replaced with the XML content at decoration time. |
| `ignore.py` | Core ignore-pattern matching logic. |
| `ignore_manager.py` | Manages the full ignore stack (built-in + `.gitignore` + user overrides). |

### Hashline Edit System

| File | Purpose |
|------|---------|
| `hashline.py` | Content-hash line tagging and validation. Provides `HashedLine`, `format_hashline()`, and `parse_line_ref()` for read/write validation. |
| `line_cache.py` | In-memory cache for edit validation. Stores `{path: {line_number: HashedLine}}` to detect stale references. |

### Discover Engine (`utils/`)

| File | Purpose |
|------|---------|
| `discover_pipeline.py` | Core discovery pipeline logic. Implements term extraction, glob generation, candidate scoring, and clustering. |
| `discover_terms.py` | Lexical vocabularies for search heuristics. Defines `SOURCE_EXTENSIONS`, `CONCEPT_EXPANSIONS`, and noise filtering. |
| `discover_types.py` | Data structures for discovery reports. Provides `DiscoveryReport`, `ConceptCluster`, `FileEntry`, and `Relevance` enum. |

### Grep Engine (`grep_components/`)

| File | Purpose |
|------|---------|
| `file_filter.py` | Decides which files to search based on ignore rules and include patterns. |
| `pattern_matcher.py` | Regex compilation and matching with timeout protection. |
| `result_formatter.py` | Formats grep results for the LLM (line numbers, context lines, truncation). |
| `search_result.py` | `SearchResult` dataclass for a single match. |

### LSP (`lsp/`)

| File | Purpose |
|------|---------|
| `client.py` | LSP client for language-server communication. |
| `diagnostics.py` | Fetch and format diagnostics from the LSP server. |
| `servers.py` | Server configuration and lifecycle management. |

### Utilities (`utils/`)

| File | Purpose |
|------|---------|
| `formatting.py` | Text formatting helpers (truncation, line numbering). |
| `ripgrep.py` | Ripgrep binary detection and invocation. |

### Cache Accessors (`cache_accessors/`)

| File | Purpose |
|------|---------|
| `ignore_manager_cache.py` | Cached ignore-manager instance. |
| `ripgrep_cache.py` | Cached ripgrep binary path. |
| `xml_prompts_cache.py` | Cached XML prompt content. |

## How

Tool registration flow:
1. `agent_config.py::_build_tools()` calls `to_tinyagent_tool()` on each decorated tool function.
2. `to_tinyagent_tool()` introspects the function signature to build an OpenAI-function JSON schema.
3. The tool's docstring (possibly replaced by XML prompt content) becomes the tool description the model sees.
4. At runtime, tinyagent calls `AgentTool.execute(tool_call_id, args, signal, on_update)`.
5. The adapter binds `args` to the function signature, checks for abort signal, calls the tool, and wraps the result in `AgentToolResult`.

Error contract:
- `ToolRetryError` -- model should try again with corrected arguments (surfaces as tool error to model).
- `ToolExecutionError` -- hard failure, reported to user.
- `FileOperationError` -- file-specific hard failure.

### UI Renderers

Each tool with visual output has a renderer in `src/tunacode/ui/renderers/`. The `hashline_edit` renderer displays diffs with syntax highlighting and line change indicators.

## Why

The decorator pattern means tool authors only write the business logic. Error handling, schema generation, and abort-signal checking are handled uniformly. The XML prompt system lets tool descriptions be edited without touching Python code.

The hashline edit system replaces the previous fuzzy matching approach with cryptographic validation, preventing edits to stale file content. Each line read by `read_file` is tagged with a content hash; `hashline_edit` validates these hashes before applying any changes.
