# Context Engineering Internals

This page documents the implementation of the ZERG context engineering system. For a conceptual overview, see [[Context Engineering]].

---

## Architecture

The context engineering system is implemented across four modules:

| Module | File | Responsibility |
|--------|------|----------------|
| Context Plugin | `zerg/context_plugin.py` | Orchestrates context building for each task |
| Command Splitter | `zerg/command_splitter.py` | Splits large command files into core + details |
| Security Rules | `zerg/security_rules.py` | Filters and summarizes security rules by file type |
| Spec Loader | `zerg/spec_loader.py` | Loads and scopes feature specs to task keywords |
| Context Tracker | `zerg/context_tracker.py` | Heuristic token counting and threshold monitoring |

The `ContextEngineeringPlugin` class in `context_plugin.py` is the main entry point. It implements the `ContextPlugin` abstract base class from the [[Plugin API Reference]] and is registered with the `PluginRegistry` at startup.

```
PluginRegistry
  |
  v
ContextEngineeringPlugin
  |-- CommandSplitter         (command splitting)
  |-- filter_rules_for_files  (security rule filtering)
  |-- summarize_rules         (security rule summarization)
  |-- SpecLoader              (spec scoping)
  |
  v
Worker receives scoped context string
```

---

## Command Splitting

### Splitting Algorithm

The `CommandSplitter` class in `zerg/command_splitter.py` handles file analysis and splitting.

**Eligibility.** Only command files with 300 or more lines are candidates for splitting. Files below this threshold are left as-is because the token savings would not justify the complexity.

**Split point selection.** The splitter:

1. Parses all `#` and `##` section headers to identify structural boundaries.
2. Calculates a target split point at approximately 30% of total lines.
3. Advances to the nearest section boundary at or after the target.
4. Ensures the split point is at least 10 lines from either end of the file.

```python
# Pseudocode for split point selection
target_line = total_lines * 0.30
for section in sections:
    cumulative_lines += section.lines
    if cumulative_lines >= target_line:
        split_at = section.end_line + 1
        break
split_at = clamp(split_at, 10, total_lines - 10)
```

**Output files.** Splitting produces three files:

| File | Content |
|------|---------|
| `{name}.core.md` | First ~30% of the original (essential instructions) |
| `{name}.details.md` | Remaining ~70% (examples, templates, edge cases) with a header |
| `{name}.md` | Updated to contain core content plus a `<!-- SPLIT: ... -->` reference comment |

The original `.md` file retains the core content so that existing symlinks from `.claude/commands/` continue to work.

### Token Estimation

The splitter estimates tokens using a simple character-based heuristic:

```python
estimated_tokens = len(text) // 4  # ~4 characters per token
```

This is intentionally conservative. Actual token counts will be lower for English prose and higher for code with special characters.

### Loading Split Files

`CommandSplitter.load_command()` is the entry point for loading commands at runtime:

1. Check if `{name}.core.md` exists. If so, load it.
2. If the caller requests details (`include_details=True`), append `{name}.details.md`.
3. If no `.core.md` exists, fall back to the full `{name}.md` file.

---

## Task-Scoped Context Building

### Entry Point

`ContextEngineeringPlugin.build_task_context()` is called for each task before the worker starts. It delegates to `_build_context_inner()`, which assembles the context from two sources:

1. **Security rules** (30% of token budget)
2. **Spec excerpts** (50% of token budget)

The remaining 20% is reserved as buffer for formatting overhead.

### Security Rule Filtering

**Module:** `zerg/security_rules.py`, function `filter_rules_for_files()`

Given a list of file paths from the task, the filter:

1. Always includes `_core/owasp-2025.md` (core OWASP rules).
2. Maps file extensions to language-specific rule files:
   - `.py`, `.pyx`, `.pyi` maps to `languages/python/CLAUDE.md`
   - `.js`, `.mjs`, `.cjs`, `.ts`, `.tsx`, `.jsx` maps to `languages/javascript/CLAUDE.md`
3. Checks for Docker-related files (`Dockerfile`, `docker-compose.*`) and maps to `containers/docker/CLAUDE.md`.
4. Returns the list of matching `Path` objects.

### Security Rule Summarization

**Module:** `zerg/security_rules.py`, function `summarize_rules()`

Full security rule files are too large for task context. The summarizer:

1. Reads each rule file.
2. Strips code blocks (content between triple-backtick fences).
3. Extracts only lines that match:
   - `## Rule:` or `### Rule:` (rule headers)
   - `**Level**:` (severity indicators)
   - `**When**:` (applicability conditions)
4. Assembles the extracted lines into a compact summary.
5. Enforces a character budget: `max_tokens * 4` characters.
6. Stops adding rule files when the budget is exhausted.

The result is a markdown summary that captures rule names, severity levels, and applicability without code examples or detailed explanations.

### Spec Scoping

**Module:** `zerg/spec_loader.py`, class `SpecLoader`

The `format_task_context()` method extracts spec sections relevant to a specific task:

**Step 1: Keyword extraction** (`_extract_task_keywords`)

- Splits the task title and description into words, keeping those longer than 3 characters.
- Extracts file path stems, splits on `_` and `-`, keeps parts longer than 2 characters.
- Returns a set of lowercase keywords.

**Step 2: Section matching** (`_extract_relevant_sections`)

- Splits the full spec text on double-newline boundaries (paragraphs).
- Scores each paragraph by counting keyword matches.
- Sorts by score (descending) and selects the top 5 paragraphs.

**Step 3: Token budgeting**

- Requirements get up to half the available token budget.
- Design gets the remainder.
- Text is truncated at paragraph or sentence boundaries when possible.
- Truncated text ends with `[... truncated for context limits ...]`.

### Fallback Behavior

If `fallback_to_full` is `true` (the default) and `build_task_context()` raises an exception, it returns an empty string. The caller interprets an empty string as "use full/global context instead."

If `fallback_to_full` is `false`, the exception propagates to the caller.

---

## Context Tracker

**Module:** `zerg/context_tracker.py`

The `ContextTracker` class monitors estimated token usage within a worker session and signals when to checkpoint.

### Token Estimation Heuristics

| Source | Estimation Formula | Constant |
|--------|-------------------|----------|
| File content | `file_size_bytes * 0.25` | ~4 chars per token |
| File read overhead | `100` per file read | Fixed overhead |
| Task context | `500` per task executed | Estimated task prompt size |
| Tool calls | `50` per invocation | Tool call framing overhead |
| Conversation growth | `100` per elapsed minute | Proxy for response accumulation |

```python
total_tokens = (
    sum(size * 0.25 + 100 for size in file_sizes)
    + tasks_executed * 500
    + tool_calls * 50
    + elapsed_minutes * 100
)
```

The maximum context window is set to 200,000 tokens.

### Checkpoint Decision

`should_checkpoint()` returns `True` when the estimated token usage exceeds the configured `context_threshold_percent` of the maximum context window:

```python
usage_percent = (estimated_tokens / 200_000) * 100
should_checkpoint = usage_percent >= threshold_percent
```

With the default threshold of 70%, this triggers at approximately 140,000 estimated tokens.

### Task Budget Allocation

`budget_for_task()` distributes a total token budget across tasks based on file count:

```python
budget_per_task = total_budget // max(file_count, 1)
clamped = max(500, min(budget_per_task, total_budget))
```

Tasks with fewer files receive a larger share of the budget. The minimum per-task budget is 500 tokens.

### Context Usage Snapshot

`get_usage()` returns a `ContextUsage` dataclass with:

| Field | Type | Description |
|-------|------|-------------|
| `estimated_tokens` | int | Current estimated token count |
| `threshold_percent` | float | Configured threshold |
| `files_read` | int | Number of files read in this session |
| `tasks_executed` | int | Number of tasks completed |
| `tool_calls` | int | Number of tool invocations |
| `usage_percent` | float | Current usage as percentage of max context |
| `is_over_threshold` | bool | Whether checkpoint should be triggered |

---

## File Collection

`ContextEngineeringPlugin._collect_task_files()` extracts the list of files a task will operate on from the task dict:

```python
def _collect_task_files(self, task: dict) -> list[str]:
    files_section = task.get("files", {})
    result = []
    for key in ("create", "modify"):
        entries = files_section.get(key)
        if isinstance(entries, list):
            result.extend(entries)
    return result
```

Only `create` and `modify` keys are considered. Read-only file references are excluded because they do not determine which security rules apply to the task's output.

---

## Token Constants Summary

| Constant | Value | Defined In |
|----------|------:|-----------|
| `CHARS_PER_TOKEN` (splitter) | 4 | `command_splitter.py` |
| `MIN_LINES_TO_SPLIT` | 300 | `command_splitter.py` |
| `TOKENS_PER_CHAR` | 0.25 | `context_tracker.py` |
| `TOKENS_PER_LINE` | 15 | `context_tracker.py` |
| `TOKENS_PER_FILE_READ` | 100 | `context_tracker.py` |
| `TOKENS_PER_TASK` | 500 | `context_tracker.py` |
| `TOKENS_PER_TOOL_CALL` | 50 | `context_tracker.py` |
| `MAX_CONTEXT_TOKENS` | 200,000 | `context_tracker.py` |
| `MAX_SPEC_TOKENS` | 2,000 | `spec_loader.py` |
| `CHARS_PER_TOKEN` (spec loader) | 4 | `spec_loader.py` |
| Default `task_context_budget_tokens` | 4,000 | `plugin_config.py` |

---

## Error Handling

Every external operation in the context engineering pipeline is wrapped in try/except to prevent context failures from blocking worker execution.

| Component | On Error | Behavior |
|-----------|----------|----------|
| `build_task_context()` | Any exception | Returns empty string (fallback to full context) if `fallback_to_full=True`; re-raises otherwise |
| `_build_security_section()` | `filter_rules_for_files` fails | Logs at debug level, returns empty string |
| `_build_spec_section()` | `SpecLoader` fails | Logs at debug level, returns empty string |
| `PluginRegistry.build_task_context()` | Individual plugin fails | Logs warning, continues with other plugins |

This design ensures that a missing spec file, a malformed task dict, or a filesystem permission error never prevents a worker from starting.

---

## See Also

- [[Context Engineering]] -- Conceptual overview and configuration
- [[Plugin API Reference]] -- ContextPlugin abstract base class
- [[Configuration]] -- Full YAML reference
- [[Tuning Guide]] -- Token budget tuning recommendations
