# Research – Permission System Removal Map

**Date:** 2026-01-26 17:34:12
**Owner:** Claude (research agent)
**Phase:** Research
**Git commit:** 48071e33c62516b5a9018bfb2b20c42a48215640
**Git branch:** master
**Repo:** alchemiststudiosDOTai/tunacode

---

## Goal

Map the entire permission system boundaries to enable clean removal. The system will be removed entirely because users are power users who use git and proper backups—permission confirmation is unnecessary overhead.

---

## Findings

### Module: `tools/authorization/` (PRIMARY REMOVAL TARGET)

The entire `src/tunacode/tools/authorization/` directory is the core of the permission system and should be **completely removed**.

| File | Lines | Purpose | Removal Action |
|------|-------|---------|----------------|
| `policy.py` | ~40 | `AuthorizationPolicy` class with `get_authorization()` and `should_confirm()` methods | DELETE |
| `handler.py` | ~65 | `ToolHandler` class orchestrating authorization checks and confirmation requests | DELETE |
| `types.py` | ~20 | `AuthorizationResult` enum (ALLOW, CONFIRM, DENY) | DELETE |
| `rules.py` | ~75 | All authorization rules (ReadOnly, YoloMode, ToolIgnore, Template) | DELETE |
| `context.py` | ~35 | `AuthContext` frozen dataclass | DELETE |
| `factory.py` | ~25 | `create_default_authorization_policy()` function | DELETE |
| `notifier.py` | ~50 | `ToolRejectionNotifier` for agent notification on rejection | DELETE |
| `requests.py` | ~130 | `ConfirmationRequestFactory` with diff preview generation | DELETE |
| `__init__.py` | ~20 | Exports all authorization components | DELETE |

**Total:** ~460 lines to remove from this directory alone.

---

### Module: `types/state.py` (Protocol Cleanup)

**File:** [`src/tunacode/types/state.py`](https://github.com/alchemiststudiosDOTai/tunacode/blob/48071e33/src/tunacode/types/state.py)

| Protocol/Type | Lines | Purpose | Removal Action |
|---------------|-------|---------|----------------|
| `TemplateProtocol` | ~10 | Template metadata for authorization (`allowed_tools` field) | REMOVE `allowed_tools` field or DELETE protocol |
| `AuthorizationProtocol` | ~10 | Authorization flow state access (includes `yolo: bool`) | DELETE |
| `ToolHandlerProtocol` | ~15 | Tool handler interface with `should_confirm()` and `process_confirmation()` | DELETE |
| `MinimalStateForAuth` | ~10 | Minimal state for authorization tools | DELETE |

---

### Module: `types/dataclasses.py` (Dataclass Cleanup)

**File:** [`src/tunacode/types/dataclasses.py`](https://github.com/alchemiststudiosDOTai/tunacode/blob/48071e33/src/tunacode/types/dataclasses.py)

| Dataclass | Lines | Purpose | Removal Action |
|-----------|-------|---------|----------------|
| `ToolConfirmationRequest` | ~15 | Request with tool_name, args, filepath, diff_content | DELETE |
| `ToolConfirmationResponse` | ~15 | Response with approved, skip_future, abort, instructions | DELETE |

---

### Module: `core/state.py` (Session State Cleanup)

**File:** [`src/tunacode/core/state.py`](https://github.com/alchemiststudiosDOTai/tunacode/blob/48071e33/src/tunacode/core/state.py)

| Field | Line | Purpose | Removal Action |
|-------|------|---------|----------------|
| `yolo: bool` | ~48 | Global auto-approve toggle | DELETE |
| `tool_ignore: list[ToolName]` | ~49 | User-selected tools to auto-approve | DELETE |
| `tool_handler: ToolHandler` | ~97 | Tool handler instance | DELETE |

---

### Module: `ui/app.py` (UI Cleanup)

**File:** [`src/tunacode/ui/app.py`](https://github.com/alchemiststudiosDOTai/tunacode/blob/48071e33/src/tunacode/ui/app.py)

| Component | Lines | Purpose | Removal Action |
|-----------|-------|---------|----------------|
| `pending_confirmation: PendingConfirmationState` | ~93 | Tracks active confirmation | DELETE |
| `request_tool_confirmation()` | 299-308 | Async method to request user confirmation | DELETE |
| `_show_inline_confirmation()` | 472-519 | Displays inline confirmation prompt with diff preview | DELETE |
| `on_key()` confirmation handling | 521-542 | Key event handler for 1/2/3/Esc keys | REMOVE confirmation branches |

---

### Module: `ui/repl_support.py` (Callback Cleanup)

**File:** [`src/tunacode/ui/repl_support.py`](https://github.com/alchemiststudiosDOTai/tunacode/blob/48071e33/src/tunacode/ui/repl_support.py)

| Component | Lines | Purpose | Removal Action |
|-----------|-------|---------|----------------|
| `build_textual_tool_callback()` | 129-151 | Creates tool callback with permission checks | SIMPLIFY to direct execution |
| `PendingConfirmationState` | 102-108 | Dataclass tracking future and request | DELETE |
| `ConfirmationRequester` protocol | 109-121 | Protocol for confirmation requests | DELETE |
| `AppForCallbacks` protocol | 123-126 | Protocol combining ConfirmationRequester | DELETE |

**Key change:** The `build_textual_tool_callback()` function needs to be simplified to execute tools directly without any confirmation logic.

---

### Module: `ui/commands/__init__.py` (YOLO Command Removal)

**File:** [`src/tunacode/ui/commands/__init__.py`](https://github.com/alchemiststudiosDOTai/tunacode/blob/48071e33/src/tunacode/ui/commands/__init__.py)

| Component | Lines | Purpose | Removal Action |
|-----------|-------|---------|----------------|
| `YoloCommand` class | 113-119 | Toggles yolo mode | DELETE |

---

### Module: `core/agents/agent_components/orchestrator/tool_dispatcher.py` (Dispatcher Cleanup)

**File:** [`src/tunacode/core/agents/agent_components/orchestrator/tool_dispatcher.py`](https://github.com/alchemiststudiosDOTai/tunacode/blob/48071e33/src/tunacode/core/agents/agent_components/orchestrator/tool_dispatcher.py)

| Component | Lines | Purpose | Removal Action |
|-----------|-------|---------|----------------|
| `dispatch_tools()` categorization | 278-337 | Separates read-only from write/execute tools | SIMPLIFY to single batch |
| `READ_ONLY_TOOLS` import | ~15 | Tool categorization | REMOVE usage |
| `UserAbortError` handling | 332-334 | Catches rejected tools | REMOVE try/except |

**Key change:** All tools can execute in parallel without categorization. No need for separate read-only vs write/execute batches.

---

### Module: `constants.py` (Constant Removal)

**File:** [`src/tunacode/constants.py`](https://github.com/alchemiststudiosDOTai/tunacode/blob/48071e33/src/tunacode/constants.py)

| Constant | Lines | Purpose | Removal Action |
|----------|-------|---------|----------------|
| `READ_ONLY_TOOLS` | 76-83 | Auto-approved tools list | DELETE |
| `WRITE_TOOLS` | 86-88 | Write tool list | DELETE (if unused) |
| `EXECUTE_TOOLS` | 91-93 | Execute tool list | DELETE (if unused) |

---

### Module: `configuration/defaults.py` (Config Cleanup)

**File:** [`src/tunacode/configuration/defaults.py`](https://github.com/alchemiststudiosDOTai/tunacode/blob/48071e33/src/tunacode/configuration/defaults.py)

| Config | Lines | Purpose | Removal Action |
|--------|-------|---------|----------------|
| `"tool_ignore": []` | ~25 | Default empty ignore list | DELETE |

---

### Module: `tests/unit/ui/test_confirmation_preview.py` (Test Removal)

**File:** [`tests/unit/ui/test_confirmation_preview.py`](https://github.com/alchemiststudiosDOTai/tunacode/blob/48071e33/tests/unit/ui/test_confirmation_preview.py)

**Action:** DELETE entire file (tests for confirmation preview truncation).

---

### Module: `ui/styles/modals.tcss` (Style Cleanup)

**File:** [`src/tunacode/ui/styles/modals.tcss`](https://github.com/alchemiststudiosDOTai/tunacode/blob/48071e33/src/tunacode/ui/styles/modals.tcss)

| Style | Purpose | Removal Action |
|-------|---------|----------------|
| `ToolConfirmationModal` | Modal styling (if present) | DELETE (unused with inline confirmation) |

---

## Data Flow Summary

```
CURRENT (WITH PERMISSIONS):
1. Agent generates tool call
2. tool_dispatcher categorizes (read-only vs write/execute)
3. tool_callback invoked → tool_handler.should_confirm()
4. AuthorizationPolicy evaluates rules (yolo, tool_ignore, read-only, template)
5. If CONFIRM: create request → show UI → await user response
6. process_response() → add to tool_ignore if skip_future
7. If rejected: notify agent via UserPromptPart
8. Tool executes or aborts

TARGET (AFTER REMOVAL):
1. Agent generates tool call
2. tool_dispatcher executes all tools in parallel
3. Tools execute directly
```

---

## Removal Strategy (Recommended Order)

### Phase 1: Core Authorization Module (DELETE)
1. Delete entire `src/tunacode/tools/authorization/` directory
2. Delete `src/tunacode/types/dataclasses.py` confirmation types
3. Delete `tests/unit/ui/test_confirmation_preview.py`

### Phase 2: Protocol Cleanup (DELETE/EDIT)
4. Edit `src/tunacode/types/state.py` - remove authorization protocols
5. Edit `src/tunacode/core/state.py` - remove yolo, tool_ignore, tool_handler fields
6. Edit `src/tunacode/constants.py` - remove tool categorization constants
7. Edit `src/tunacode/configuration/defaults.py` - remove tool_ignore config

### Phase 3: UI Cleanup (DELETE/EDIT)
8. Edit `src/tuna/tunacode/ui/app.py` - remove confirmation UI methods and key handlers
9. Edit `src/tunacode/ui/repl_support.py` - simplify build_textual_tool_callback()
10. Edit `src/tunacode/ui/commands/__init__.py` - remove YoloCommand
11. Edit `src/tunacode/ui/styles/modals.tcss` - remove confirmation styles (if present)

### Phase 4: Dispatcher Simplification (EDIT)
12. Edit `src/tunacode/core/agents/agent_components/orchestrator/tool_dispatcher.py` - remove categorization, execute all tools in parallel

### Phase 5: Documentation Update (DELETE/EDIT)
13. Update `docs/codebase-map/modules/tools-overview.md` - remove permission references
14. Update `docs/codebase-map/structure/03-tools-directory.md` - remove authorization directory

---

## Key Dependencies to Break

| From | To | Break Action |
|------|-----|--------------|
| `ui/app.py` | `types.dataclasses.ToolConfirmationRequest` | Remove import and usage |
| `ui/app.py` | `types.dataclasses.ToolConfirmationResponse` | Remove import and usage |
| `ui/repl_support.py` | `types.state.ConfirmationRequester` | Remove protocol |
| `ui/repl_support.py` | `types.state.ToolHandlerProtocol` | Remove protocol |
| `core/state.py` | `tools.authorization.ToolHandler` | Remove field |
| `core/agents/agent_components/orchestrator/tool_dispatcher.py` | `constants.READ_ONLY_TOOLS` | Remove usage |
| `ui/commands/__init__.py` | `core.state.AuthorizationProtocol` | Remove YoloCommand |

---

## Knowledge Gaps

- **Template system integration:** The `TemplateAllowedToolsRule` references templates. Need to verify if templates have other uses beyond authorization.
- **UserAbortError propagation:** Need to trace where else `UserAbortError` might be used outside permission flow.
- **Tool retry logic:** Verify if tool retry interacts with permission system.

---

## References

### Core Files (Primary Removal Targets)
- [`src/tunacode/tools/authorization/`](https://github.com/alchemiststudiosDOTai/tunacode/tree/48071e33/src/tunacode/tools/authorization/) - Entire directory (~9 files, ~460 LOC)
- [`src/tunacode/types/state.py`](https://github.com/alchemiststudiosDOTai/tunacode/blob/48071e33/src/tunacode/types/state.py) - Authorization protocols
- [`src/tunacode/types/dataclasses.py`](https://github.com/alchemiststudiosDOTai/tunacode/blob/48071e33/src/tunacode/types/dataclasses.py) - Confirmation types

### UI Files (Secondary Removal Targets)
- [`src/tunacode/ui/app.py`](https://github.com/alchemiststudiosDOTai/tunacode/blob/48071e33/src/tunacode/ui/app.py) - Confirmation UI
- [`src/tunacode/ui/repl_support.py`](https://github.com/alchemiststudiosDOTai/tunacode/blob/48071e33/src/tunacode/ui/repl_support.py) - Tool callback
- [`src/tunacode/ui/commands/__init__.py`](https://github.com/alchemiststudiosDOTai/tunacode/blob/48071e33/src/tunacode/ui/commands/__init__.py) - YoloCommand

### Core Files (Edits Required)
- [`src/tunacode/core/state.py`](https://github.com/alchemiststudiosDOTai/tunacode/blob/48071e33/src/tunacode/core/state.py) - Session state fields
- [`src/tunacode/core/agents/agent_components/orchestrator/tool_dispatcher.py`](https://github.com/alchemiststudiosDOTai/tunacode/blob/48071e33/src/tunacode/core/agents/agent_components/orchestrator/tool_dispatcher.py) - Tool dispatcher
- [`src/tunacode/constants.py`](https://github.com/alchemiststudiosDOTai/tunacode/blob/48071e33/src/tunacode/constants.py) - Tool categorization

### Test Files (Delete)
- [`tests/unit/ui/test_confirmation_preview.py`](https://github.com/alchemiststudiosDOTai/tunacode/blob/48071e33/tests/unit/ui/test_confirmation_preview.py) - Confirmation tests

### Documentation (Update)
- [`docs/codebase-map/modules/tools-overview.md`](https://github.com/alchemiststudiosDOTai/tunacode/blob/48071e33/docs/codebase-map/modules/tools-overview.md)
- [`docs/codebase-map/structure/03-tools-directory.md`](https://github.com/alchemiststudiosDOTai/tunacode/blob/48071e33/docs/codebase-map/structure/03-tools-directory.md)
