---
title: Tools Module
path: src/tunacode/tools
type: directory
depth: 1
description: Tool implementations and decorators
exports: [base_tool, file_tool, bash, glob, grep, read_file, write_file, update_file]
seams: [M, D]
---

# Tools Module

## Purpose
Provides AI agent capabilities through a robust tool system with decorators and specialized implementations.

## Tool Decorator System (decorators.py)

### base_tool
Foundation decorator for all tools:
- Standardizes error handling (Exception → ToolExecutionError)
- Loads tool descriptions from XML prompt files
- Ensures consistent error propagation
- Preserves ModelRetry and ToolExecutionError exceptions

### file_tool
Specialized decorator for file system tools:
- Extends base_tool functionality
- Handles file-specific errors (FileNotFoundError, PermissionError, UnicodeDecodeError)
- Maps errors to appropriate exceptions (ModelRetry for recoverable, FileOperationError for fatal)

**XML Prompt Files:**
- `bash_prompt.xml` - Bash tool documentation
- `read_file_prompt.xml` - Read file documentation
- `submit_prompt.xml` - Submit tool documentation
- etc.

## Tool Implementations

### bash.py
Executes shell commands securely:
- Destructive pattern validation (blocks `rm -rf`, etc.)
- Configurable timeouts
- Stdout/stderr capture with exit codes
- Working directory management

### glob.py
Fast file pattern matching:
- Supports `**/*.py` style patterns
- Recursive searching
- .gitignore awareness
- Sorted results
- Shared ignore manager in `src/tunacode/tools/ignore.py`

### grep.py
Advanced content search:
**Search Strategies:**
- `smart` - Auto-selects best method
- `ripgrep` - Uses ripgrep binary
- `python` - Pure Python implementation
- `hybrid` - Combined approach

**Features:**
- Regex support
- Case sensitivity control
- File type filtering
- Context lines
- Output modes: content, files-only, count, JSON
- Broad pattern prevention
- Shared ignore manager for fast-glob prefiltering and `.gitignore` rules

### read_file.py
Reads file contents safely:
- Line range support (offset/limit)
- Line number preservation
- Long line truncation
- Streaming line reads to avoid loading entire files into memory
- File size limits
- Encoding handling

### write_file.py
Creates new files:
- Parent directory auto-creation
- Overwrite prevention (use update_file instead)
- LSP diagnostic triggers
- Atomic writes

### update_file.py
Modifies existing files:
**Fuzzy Matching Algorithms:**
- `line-trimmed` - Whitespace-tolerant
- `indentation-flexible` - Indent-tolerant
- `block-anchor` - Multi-line block matching

**Features:**
- Diff generation
- LSP diagnostic updates
- Change validation

### list_dir.py
Lists directory contents:
- Recursive tree view
- Ignore pattern support
- Hidden file visibility
- Output size limits
- Shared ignore manager in `src/tunacode/tools/ignore.py`

### ignore.py
Shared ignore manager for discovery tools:
- Loads default ignore patterns plus root `.gitignore`
- Root-only `.gitignore` support (no nested ignore files)
- Pathspec `gitwildmatch` matching for ignore rules
- Provides `should_ignore`/`should_ignore_dir` helpers and cached manager instances

### web_fetch.py
Fetches web content:
- HTTP/HTTPS only
- HTML to plain text conversion
- Security validation (blocks private IPs)
- Content size limits
- Timeout handling

### lsp/ (sub-package)
LSP integration for file tools:
- `client.py` - Minimal LSP client (JSON-RPC over stdin/stdout)
- `servers.py` - Server config mapping (extension → language server)
- `diagnostics.py` - `maybe_prepend_lsp_diagnostics()` helper for write/update tools

### submit.py
Completion signaling tool:
- `submit` marks the task as complete for the orchestrator
- Optional summary text for the final response

## Tool Utilities (tools/utils/)

### formatting.py
Diagnostic message formatting helpers:
- **truncate_diagnostic_message()** - Shortens verbose diagnostic text for LSP output

### ripgrep.py
Ripgrep integration wrapper:
- Binary detection
- Output parsing
- Error handling
- Async subprocess execution with timeouts

### text_match.py
Fuzzy matching algorithms for update_file:
- Line-trimmed matching
- Indentation-flexible matching
- Block-anchor matching

## Integration Points

- **core/agents/** - Tool registration and execution
- **core/agents/tool_executor.py** - Parallel execution
- **ui/renderers/** - Specialized output rendering

## Seams (M, D)

**Modification Points:**
- Add new tool implementations
- Extend decorator system with new features
- Add new fuzzy matching algorithms

**Extension Points:**
- Create custom tool types
- Add tool-specific validation
- Extend prompt XML system
