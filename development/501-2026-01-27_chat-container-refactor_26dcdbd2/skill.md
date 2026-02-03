---
title: "Chat Container Refactor – Plan"
phase: Plan
date: "2026-01-27"
owner: "agent"
parent_research: "inline (user-provided research in command args)"
git_commit_at_plan: "5efc5423"
tags: [plan, ui, streaming, chat, refactor]
---

## Goal

- **ONE outcome:** Replace dual-display streaming architecture (RichLog + streaming_output Static) with a unified ChatContainer using individual MessageWidget components that stream in-place.

**Non-goals:**
- Deployment/release process
- Performance benchmarking beyond visual verification
- Tool result rendering changes (keep current panel rendering – deferred decision)

## Scope & Assumptions

**In scope:**
- New `ChatContainer` and `MessageWidget` widgets
- Replace `RichLog` with `ChatContainer` in `app.py`
- Simplify streaming callback to direct widget updates
- Update all `rich_log.write()` call sites to use `chat.add_message()`
- Update CSS for new widget classes

**Out of scope:**
- Tool result widget changes (keep `ToolResultDisplay` → `tool_panel_smart` flow)
- Session replay refactoring (will use new API)
- Logger TUI callback (will use new API)

**Assumptions:**
- Textual's `VerticalScroll` + `Static.update()` handles batched refresh efficiently
- Plain text streaming is acceptable (Markdown only on finalize) – matches Claude.ai, ChatGPT, Cursor
- Clean break – no `rich_log` alias during transition

## Deliverables

1. `src/tunacode/ui/widgets/chat.py` – New `MessageWidget` and `ChatContainer` classes
2. Updated `src/tunacode/ui/app.py` – Replace RichLog with ChatContainer
3. Updated `src/tunacode/ui/commands/__init__.py` – Use `chat.add_message()` API
4. Updated `src/tunacode/ui/welcome.py` – Use `chat.add_message()` API
5. Updated `src/tunacode/ui/styles/layout.tcss` – CSS for new widget classes

## Readiness

**Preconditions:**
- [x] Research complete (user-provided analysis)
- [x] Git state captured (5efc5423)
- [x] Touch points identified (26 `rich_log.write()` calls across 3 files)
- [x] User decisions: clean break, NeXTSTEP styling, tool results unchanged

## Milestones

- **M1:** Widget skeleton and CSS foundation
- **M2:** Core streaming integration (app.py refactor)
- **M3:** Call site migration (commands, welcome, logging)
- **M4:** Visual verification and edge cases

## Work Breakdown (Tasks)

### M1: Widget Skeleton and CSS Foundation

| ID | Task | Files | Acceptance Test |
|----|------|-------|-----------------|
| T1 | Create `MessageWidget` class with `append()` and `finalize()` methods | `widgets/chat.py` | Widget can be mounted, content appended, then finalized to Markdown panel |
| T2 | Create `ChatContainer` class with `add_message()`, `start_stream()`, `stream()`, `end_stream()` methods | `widgets/chat.py` | Container scrolls, can add static messages and stream new ones |
| T3 | Add CSS rules for `MessageWidget` and `ChatContainer` | `styles/layout.tcss` | Widgets render with proper height/margin/scrolling |
| T4 | Export new widgets from `widgets/__init__.py` | `widgets/__init__.py` | Import works from `tunacode.ui.widgets` |

**Dependencies:** None (foundation)

### M2: Core Streaming Integration

| ID | Task | Files | Acceptance Test |
|----|------|-------|-----------------|
| T5 | Add `chat: ChatContainer` to `TextualReplApp` | `app.py` | App composes with ChatContainer instead of RichLog |
| T6 | Remove streaming-related state (`streaming_output`, `current_stream_text`, `_last_display_update`, `STREAM_THROTTLE_MS`, `_update_streaming_panel`) | `app.py` | All removed state attributes and methods deleted |
| T7 | Simplify `streaming_callback` to call `chat.stream()` | `app.py` | Tokens appear immediately (no throttle) |
| T8 | Refactor `_process_request` to use `chat.start_stream()` / `chat.end_stream()` | `app.py` | Stream starts before request, finalizes in finally block |
| T9 | Update user message rendering to use `chat.add_message()` | `app.py:on_editor_submit_requested` | User messages appear as MessageWidgets |
| T10 | Update tool result rendering to mount into ChatContainer | `app.py:on_tool_result_display` | Tool panels appear in chat flow |

**Dependencies:** T1-T4

### M3: Call Site Migration

| ID | Task | Files | Acceptance Test |
|----|------|-------|-----------------|
| T11 | Update `/help` command to use `chat.add_message()` | `commands/__init__.py:69` | Help table renders in chat |
| T12 | Update `/clear` command to call `chat.clear()` | `commands/__init__.py:79` | Chat clears properly |
| T13 | Update `/debug` command output | `commands/__init__.py:132` | Debug messages appear in chat |
| T14 | Update `/model` command error messages | `commands/__init__.py:163` | API key warnings in chat |
| T15 | Update `/resume` session loaded message | `commands/__init__.py:385` | Load confirmation in chat |
| T16 | Update `/update` command outputs | `commands/__init__.py:411-460` | Version info in chat |
| T17 | Refactor `show_welcome()` to accept ChatContainer | `welcome.py` | Welcome renders in ChatContainer |
| T18 | Update logger TUI callback to use `chat.add_message()` | `app.py:165-168` | Debug logs appear in chat |
| T19 | Update `_replay_session_messages()` to use chat API | `app.py:310-331` | Session restore works |
| T20 | Update `write_shell_output()` to use chat | `app.py:403-404` | Shell output in chat |
| T21 | Update error rendering to use chat | `app.py:182,218` | Errors appear in chat |
| T22 | Update system notices to use chat | `app.py:299-301` | Notices in chat |

**Dependencies:** T5-T10

### M4: Visual Verification and Edge Cases

| ID | Task | Files | Acceptance Test |
|----|------|-------|-----------------|
| T23 | Handle pause/resume streaming with ChatContainer | `app.py` | Ctrl+P pauses/resumes, buffer flush works |
| T24 | Handle cancel streaming (Escape) | `app.py` | Cancel clears in-progress message |
| T25 | Remove `#streaming-output` CSS and `#viewport` streaming classes | `styles/layout.tcss` | No dead CSS rules |
| T26 | Manual visual verification | N/A | Streaming immediate, finalize styled, scroll works |

**Dependencies:** T11-T22

## Risks & Mitigations

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| Scroll jank during rapid token updates | Low | Medium | Textual batches refresh; if needed, add minimal throttle (16ms) |
| MessageWidget height calculation issues | Medium | Low | Use `height: auto` CSS; test with long/short messages |
| Tool panels break in new container | Low | High | Keep existing `tool_panel_smart` flow unchanged; just mount to ChatContainer |

## Test Strategy

- **T26:** Manual visual verification (streaming, finalize, scroll, pause, cancel)
- No unit tests required – this is pure UI refactoring with no business logic changes

## References

- Research: User-provided inline analysis of dual-display architecture
- Current streaming: `app.py:333-356` (`streaming_callback`, `_update_streaming_panel`)
- Current finalization: `app.py:228-243` (finally block in `_process_request`)
- CSS: `styles/layout.tcss:66-81` (`#streaming-output` rules)
- Touch points: 26 `rich_log.write()` calls (Grep results above)

## Final Gate

| Item | Value |
|------|-------|
| Plan path | `memory-bank/plan/2026-01-27_chat-container-refactor.md` |
| Milestone count | 4 |
| Task count | 26 |
| Ready for coding | Yes |

**Next command:** `/context-engineer:execute "memory-bank/plan/2026-01-27_chat-container-refactor.md"`
