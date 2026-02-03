# NeXTSTEP Panel Uniformity & Architecture

**Date:** 2025-12-06
**Scope:** UI / Tooling
**Status:** Canonical

## Core Philosophy

The TunaCode UI follows the **NeXTSTEP** design philosophy: **Information Hierarchy**, **Direct Manipulation**, and **Visual Feedback**.

To maintain a professional, predictable, and high-quality user experience, all panels—especially tool outputs—must strictly adhere to a uniform architectural pattern. We do not create "special" snowflake panels for individual tools unless they are extensions of the standard base.

## The Panel Standard

Every tool execution in the system, regardless of its function (search, file edit, shell command), must be rendered using the **Standard Tool Panel** layout.

### Anatomy of a Tool Panel

```
┌──────────────────────────────────────────────────────────────────┐
│  TOOL_NAME [STATUS]                                              │ ← Header
├──────────────────────────────────────────────────────────────────┤
│                                                                  │
│  [Arguments Section]                                             │ ← Context
│  filepath: src/main.py                                           │
│  pattern:  "def hello"                                           │
│                                                                  │
│  [Content Section]                                               │ ← The "Hero"
│  (The primary output of the tool. Text, Diff, Table, etc.)       │
│                                                                  │
│  [Footer Section]                                                │ ← Metadata
│  145ms • 12 lines                                                │
│                                                                  │
└──────────────────────────────────────────────────────────────────┘
   09:41 AM                                                          ← Subtitle
```

### 1. The Header
*   **Left:** Tool Name (e.g., `update_file`, `grep`).
*   **Right:** Status Indicator (`[running]`, `[done]`, `[fail]`).
*   **Style:** Primary color for tool name; status color-coded (Green=Success, Red=Error).

### 2. The Context (Arguments)
*   **Must always be present.**
*   Displays the *inputs* that produced the result.
*   **Format:** Key-Value table, truncated if values are too long (>60 chars).
*   **Why:** Without context, the result is meaningless. The user must see *what* was asked to understand *what* is returned.

### 3. The Content (Hero)
*   This is the *only* variable region.
*   **Text Tools:** Standard text output.
*   **Diff Tools:** Syntax-highlighted unified diff (using `rich.syntax.Syntax`).
*   **Search Tools:** Structured list of matches.
*   **Constraint:** Must fit within the panel boundaries. Large content is truncated/paginated.

### 4. The Footer (Metadata)
*   **Duration:** Execution time in milliseconds (e.g., `145ms`).
*   **Metrics:** Line counts, token usage, or other quantitative data.
*   **Style:** Dimmed/Muted text.

### 5. The Subtitle
*   **Timestamp:** The wall-clock time the action occurred.
*   **Placement:** Bottom border of the panel.

## Implementation Guidelines

### Python Architecture
We enforce uniformity via the `RichPanelRenderer` class in `src/tunacode/ui/renderers/panels.py`.

*   **Do not** create ad-hoc `Panel()` objects in random parts of the code.
*   **Do not** manually build strings for standard metadata.
*   **Always** use `render_tool`, `render_error`, or `render_diff_tool`.

### The `render_diff_tool` Spec
When rendering a diff (e.g., for `update_file` or `write_file`), you must pass the standard context:

```python
RichPanelRenderer.render_diff_tool(
    tool_name="update_file",
    message="File updated successfully",
    diff="--- a/file.py\n+++ b/file.py...",
    args={"filepath": "...", "target": "..."},  # <--- CRITICAL: Must include args
    duration_ms=120,                            # <--- CRITICAL: Must include metrics
    timestamp=datetime.now()                    # <--- CRITICAL: Must include timestamp
)
```

**Anti-Pattern (Forbidden):**
```python
# WRONG: Missing context and metadata
Panel(
    Syntax(diff, "diff"),
    title="update_file"
)
```

## Tool Conformance

All tools must return data that allows the UI to populate this structure.
*   **Return Type:** Tools return `str` (enforced by `test_tool_conformance.py`).
*   **Rich Content:** If a tool needs to return rich data (like a diff), it must serialise it into the string (e.g., appending the diff after a double newline), which the UI renderer then parses.

## Verification

Any changes to panel rendering must be verified against:
1.  **Visual Consistency:** Does it look like the other panels?
2.  **Information Completeness:** Are args, duration, and timestamp visible?
3.  **NeXTSTEP Alignment:** Does it provide immediate visual feedback of the action?
