---
title: "Claude Code Release History"
description: "Condensed changelog of official Claude Code releases with highlights and breaking changes"
tags: [reference, release]
---

# Claude Code Release History

> Condensed changelog of Claude Code official releases.
> **Full details**: [github.com/anthropics/claude-code/CHANGELOG.md](https://github.com/anthropics/claude-code/blob/main/CHANGELOG.md)
> **Machine-readable**: [claude-code-releases.yaml](../machine-readable/claude-code-releases.yaml)

**Latest**: v2.1.59 | **Updated**: 2026-02-26

---

## Quick Jump

- [2.1.x Series (January-February 2026)](#21x-series-january-february-2026) — Worktree isolation, background agents, ConfigChange hook, Fast mode Opus 4.6, 1M context, claude.ai MCP connectors, remote-control, auto-memory, /copy command
- [2.0.x Series (Nov 2025 - Jan 2026)](#20x-series-november-2025---january-2026) — Opus 4.5, Claude in Chrome, Background agents
- [Breaking Changes Summary](#breaking-changes-summary)
- [Milestone Features](#milestone-features)

---

## 2.1.x Series (January-February 2026)

### v2.1.59 (2026-02-26)

- **New**: Auto-memory — Claude automatically saves useful context to memory; manage with `/memory`
- **New**: `/copy` command — interactive picker when code blocks are present, select individual blocks or full response
- **Improved**: Smarter "always allow" prefix suggestions for compound bash commands (per-subcommand prefixes instead of treating whole command as one)
- **Improved**: Memory usage in multi-agent sessions (releases completed subagent task state)
- **Improved**: Ordering of short task lists
- **Fixed**: MCP OAuth token refresh race condition when running multiple Claude Code instances simultaneously
- **Fixed**: Shell commands not showing clear error message when working directory has been deleted
- **Fixed**: Config file corruption that could wipe authentication when multiple Claude Code instances ran simultaneously

### v2.1.58 (2026-02-26)

- **Expanded**: Remote Control available to more users

### v2.1.56 (2026-02-25)

- **Fixed**: VSCode: Another cause of "command 'claude-vscode.editor.openLast' not found" crashes

### v2.1.55 (2026-02-25)

- **Fixed**: BashTool failing on Windows with EINVAL error

### v2.1.53 (2026-02-25)

- **Fixed**: UI flicker where user input briefly disappeared after submission before rendering
- **Fixed**: Bulk agent kill (ctrl+f) now sends single aggregate notification instead of one per agent, and properly clears command queue
- **Fixed**: Graceful shutdown sometimes leaving stale sessions when using Remote Control (parallelized teardown)
- **Fixed**: `--worktree` flag sometimes being ignored on first launch
- **Fixed**: Panic ("switch on corrupted value") on Windows
- **Fixed**: Crash when spawning many processes on Windows
- **Fixed**: Crash in WebAssembly interpreter on Linux x64 & Windows x64
- **Fixed**: Crash that sometimes occurred after 2 minutes on Windows ARM64

### v2.1.52 (2026-02-24)

- **Fixed**: VSCode extension crash on Windows ("command 'claude-vscode.editor.openLast' not found")

### v2.1.51 (2026-02-24)

- **New**: `claude remote-control` subcommand for external builds — enables local environment serving for all users
- **New**: Custom npm registries and specific version pinning when installing plugins from npm sources
- **New**: SDK: `CLAUDE_CODE_ACCOUNT_UUID`, `CLAUDE_CODE_USER_EMAIL`, `CLAUDE_CODE_ORGANIZATION_UUID` env vars to provide account info synchronously (eliminates race conditions in early telemetry)
- **Changed**: BashTool now skips login shell (`-l` flag) by default when shell snapshot is available — performance improvement (previously required `CLAUDE_BASH_NO_LOGIN=true`)
- **Changed**: Tool results larger than 50K characters now persisted to disk (previously 100K threshold)
- **Improved**: `/model` picker now shows human-readable labels (e.g., "Sonnet 4.5") instead of raw model IDs for pinned versions, with upgrade hint when newer version available
- **Fixed**: Security issue where `statusLine` and `fileSuggestion` hook commands could execute without workspace trust acceptance in interactive mode
- **Fixed**: Duplicate `control_response` messages from WebSocket reconnects causing API 400 errors
- **Fixed**: Slash command autocomplete crashing when a plugin's SKILL.md description is a YAML array or other non-string type

### v2.1.50 (2026-02-21)

- **New**: `WorktreeCreate` and `WorktreeRemove` hook events — custom VCS setup/teardown when agent worktree isolation creates or removes worktrees
- **New**: `isolation: worktree` in agent definitions for declarative worktree isolation (no longer requires setting in each call)
- **New**: `claude agents` CLI command to list all configured agents
- **New**: `startupTimeout` configuration for LSP servers
- **New**: `CLAUDE_CODE_DISABLE_1M_CONTEXT` env var to disable 1M context window support
- **New**: Pre-configured OAuth client credentials for MCP servers that don't support Dynamic Client Registration (Slack); use `--client-id` and `--client-secret` with `claude mcp add`
- **New**: VSCode `/extra-usage` command support
- **Changed**: Opus 4.6 (fast mode) now includes full 1M context window
- **Changed**: `CLAUDE_CODE_SIMPLE` mode now also disables MCP tools, attachments, hooks, and CLAUDE.md loading for fully minimal experience
- **Fixed**: Bug where resumed sessions could be invisible when working directory involved symlinks
- **Fixed**: `disableAllHooks` setting to respect managed settings hierarchy (non-managed settings can no longer disable managed hooks)
- **Fixed**: Linux: native modules not loading on systems with glibc older than 2.30 (RHEL 8)
- **Fixed**: Memory leak in agent teams where completed teammate tasks were never garbage collected
- **Fixed**: Memory leak where completed task state objects were never removed from AppState
- **Fixed**: Memory leak where LSP diagnostic data was never cleaned up after delivery
- **Fixed**: Unbounded memory growth in long sessions (file history snapshots capped; circular buffer fix; stream buffers released after use)
- **Fixed**: MCP tools not discovered when tool search is enabled and prompt passed as launch argument
- **Fixed**: Prompt suggestion cache regression that reduced cache hit rates
- **Improved**: Startup performance for headless mode (`-p`) by deferring Yoga WASM and UI component imports
- **Improved**: Memory usage during long sessions by clearing internal caches after compaction and clearing large tool results after processing

### v2.1.49 (2026-02-20)

- **New**: `--worktree` / `-w` CLI flag to start Claude in an isolated git worktree
- **New**: Subagents support `isolation: "worktree"` for working in a temporary git worktree
- **New**: `background: true` field in agent definitions to always run as a background task
- **New**: `ConfigChange` hook event — fires when configuration files change during a session (enterprise security auditing + blocking)
- **New**: Plugins can ship `settings.json` for default configuration
- **New**: `--from-pr` flag to resume sessions linked to a specific GitHub PR (+ sessions auto-linked when created via `gh pr create`)
- **New**: `PreToolUse` hooks can return `additionalContext` to the model
- **New**: `plansDirectory` setting to customize where plan files are stored
- **New**: `auto:N` syntax for configuring MCP tool search auto-enable threshold
- **New**: `Setup` hook event triggered via `--init`, `--init-only`, or `--maintenance` CLI flags
- **Changed**: Sonnet 4.5 1M context removed from Max plan — Sonnet 4.6 now has 1M context (switch in `/model`)
- **Changed**: Simple mode now includes file edit tool (not just Bash)
- **Fixed**: File-not-found errors now suggest corrected paths when model drops repo folder
- **Fixed**: Ctrl+C and ESC silently ignored when background agents running + main thread idle (double-press within 3s now kills all agents)
- **Fixed**: Plugin `enable`/`disable` auto-detects correct scope (no longer defaults to user scope)
- **Fixed**: Context window blocking limit calculated too aggressively (~65% instead of ~98%)
- **Fixed**: Memory issues causing crashes with parallel subagents
- **Fixed**: Memory leak in long sessions where stream resources not cleaned up
- **Fixed**: `@` symbol incorrectly triggering file autocomplete in bash mode
- **Fixed**: Background agent results returning raw transcript data instead of final answer
- **Fixed**: Slash command autocomplete selecting wrong command (e.g. `/context` vs `/compact`)
- **Improved**: `@` mention file suggestion speed (~3× faster in git repos)
- **Improved**: MCP connection: `list_changed` notification support for dynamic tool updates without reconnection
- **Improved**: Skills invoke progress display; skill suggestions prioritize recently/frequently used
- **Improved**: Incremental output for async agents; token count includes background agent tokens

### v2.1.47 (2026-02-19)

- **Improved**: VS Code plan preview auto-updates as Claude iterates; commenting enabled only when plan is ready for review; preview stays open when rejected for revision
- **New**: `ctrl+f` kills all background agents simultaneously (replaces double-ESC); ESC now cancels main thread only, background agents keep running
- **New**: `last_assistant_message` field added to Stop and SubagentStop hook inputs (access final response without parsing transcript files)
- **New**: `chat:newline` keybinding action; `added_dirs` in statusline JSON workspace section
- **Fixed**: Compaction failing when conversation contains many PDF documents (strips document blocks alongside images)
- **Fixed**: Edit tool corrupting Unicode curly quotes (`"` `"` `'` `'`) by replacing with straight quotes
- **Fixed**: Parallel file write/edit — single file failure no longer aborts sibling operations
- **Fixed**: OSC 8 hyperlinks only clickable on first line when link text wraps across multiple terminal lines
- **Fixed**: Bash permission classifier now validates match descriptions against actual input rules (prevents hallucinated permissions)
- **Fixed**: Config backups timestamped and rotated (5 most recent kept) instead of overwriting
- **Fixed**: Session name lost after context compaction; plan mode lost after compaction
- **Fixed**: Hooks (PreToolUse, PostToolUse) silently failing on Windows (now uses Git Bash)
- **Fixed**: Custom agents/skills not discovered in git worktrees (main repo `.claude/` now included)
- **Fixed**: 70+ additional rendering, session, permission, and platform fixes

### v2.1.46 (2026-02-19)

- **Fixed**: Orphaned Claude Code processes after terminal disconnect on macOS
- **New**: Support for using claude.ai MCP connectors in Claude Code

### v2.1.45 (2026-02-17)

- **New**: Claude Sonnet 4.6 model support
- **New**: `spinnerTipsOverride` setting — customize spinner tips via `tips` array, opt out of built-in tips with `excludeDefault: true`
- **New**: SDK `SDKRateLimitInfo` and `SDKRateLimitEvent` types for rate limit status tracking (utilization, reset times, overage)
- **Fixed**: Agent Teams teammates failing on Bedrock, Vertex, and Foundry (env vars now propagated to tmux-spawned processes)
- **Fixed**: Sandbox "operation not permitted" errors on macOS temp file writes
- **Fixed**: Task tool (backgrounded agents) crashing with `ReferenceError` on completion
- **Improved**: Memory usage for large shell command outputs (RSS no longer grows unboundedly)
- **Improved**: Startup performance (removed eager session history loading)
- **Improved**: Plugin-provided commands, agents, and hooks available immediately after install (no restart needed)

### v2.1.44 (2026-02-17)

- Fixed: Auth refresh errors

### v2.1.43 (2026-02-17)

- Fixed: AWS auth refresh hanging indefinitely (added 3-minute timeout)
- Fixed: Structured-outputs beta header being sent unconditionally on Vertex/Bedrock
- Fixed: Spurious warnings for non-agent markdown files in `.claude/agents/` directory

### v2.1.42 (2026-02-14)

- **Improved**: Startup performance via deferred Zod schema construction (faster on large projects)
- **Improved**: Prompt cache hit rate by moving date outside the system prompt (avoids daily cache invalidation)
- **New**: Opus 4.6 effort callout for eligible users (one-time onboarding)
- Fixed: `/resume` showing interrupt messages as session titles
- Fixed: Image dimension limit errors now suggest using `/compact` instead of opaque failure

### v2.1.41 (2026-02-13)

- **New**: Guard against launching Claude Code inside another Claude Code session
- **New**: `claude auth login`, `claude auth status`, `claude auth logout` CLI subcommands
- **New**: Windows ARM64 (win32-arm64) native binary support
- Added `speed` attribute to OTel events and trace spans for fast mode visibility
- **Improved**: `/rename` auto-generates session name from conversation context when called without arguments
- Improved narrow terminal layout for prompt footer
- Fixed: Agent Teams using wrong model identifier for Bedrock, Vertex, and Foundry customers
- Fixed: Crash when MCP tools return image content during streaming
- Fixed: `/resume` session previews showing raw XML tags instead of readable command names
- Fixed: Opus 4.6 launch announcement showing for Bedrock/Vertex/Foundry users
- Fixed: Hook blocking errors (exit code 2) not showing stderr to the user
- Fixed: Structured-outputs beta header sent unconditionally on Vertex/Bedrock
- Fixed: File resolution for @-mentions with anchor fragments (e.g., `@README.md#installation`)
- Fixed: FileReadTool blocking on FIFOs, `/dev/stdin`, and large files
- Fixed: Background task notifications not delivered in streaming Agent SDK mode
- Fixed: Auto-compact failure error notifications shown to users
- Fixed: Stale permission rules not clearing when settings change on disk
- Fixed: Permission wait time included in subagent elapsed time display
- Fixed: Proactive ticks firing while in plan mode
- Improved: Model error messages for Bedrock/Vertex/Foundry with fallback suggestions

### v2.1.39 (2026-02-10)

- Improved: Terminal rendering performance
- Fixed: Fatal errors being swallowed instead of displayed
- Fixed: Process hanging after session close
- Fixed: Character loss at terminal screen boundary
- Fixed: Blank lines in verbose transcript view

### v2.1.38 (2026-02-10)

- Fixed: VS Code terminal scroll-to-top regression introduced in 2.1.37
- Fixed: Tab key queueing slash commands instead of autocompleting
- Fixed: Bash permission matching for commands using environment variable wrappers
- Fixed: Text between tool uses disappearing when not using streaming
- **Security**: Improved heredoc delimiter parsing to prevent command smuggling
- **Security**: Blocked writes to `.claude/skills` directory in sandbox mode

### v2.1.37 (2026-02-08)

- Fixed `/fast` not immediately available after enabling `/extra-usage`

### v2.1.36 (2026-02-08) ⭐

- ⭐ **Fast mode now available for Opus 4.6** — Same model, faster output. Toggle with `/fast`. [Learn more](https://code.claude.com/docs/en/fast-mode)

### v2.1.34 (2026-02-07)

- Fixed a crash when agent teams setting changed between renders
- **Security fix**: Commands excluded from sandboxing (via `sandbox.excludedCommands` or `dangerouslyDisableSandbox`) could bypass the Bash ask permission rule when `autoAllowBashIfSandboxed` was enabled

### v2.1.33 (2026-02-06)

**Highlights**:
- **Agent teams fixes** — Improved tmux session handling and availability warnings
- **New hook events** — `TeammateIdle` and `TaskCompleted` for multi-agent workflows
- **Agent frontmatter enhancements**:
  - `memory` field for user/project/local scope memory selection
  - `Task(agent_type)` syntax to restrict sub-agent spawning in agent definitions
- **Plugin identification** — Plugin name now shown in skill descriptions and `/skills` menu
- **VSCode improvements** — Remote sessions support, branch/message count in session picker
- Fixed: Thinking interruption, streaming abort, proxy settings, `/resume` XML markup
- Improved: API connection errors show specific cause instead of generic message
- Improved: Invalid managed settings errors now surfaced properly
- Multiple stability fixes across agent workflows and tool interactions

### v2.1.32 (2026-02-05) ⭐ MAJOR

**Highlights**:
- ⭐ **Claude Opus 4.6 is now available!**
- ⭐ **Agent teams research preview** — Multi-agent collaboration for complex tasks (token-intensive, requires `CLAUDE_CODE_EXPERIMENTAL_AGENT_TEAMS=1`)
- ⭐ **Automatic memory recording and recall** — Claude now automatically records and recalls memories as it works
- **"Summarize from here"** — Message selector now allows partial conversation summarization
- Skills from `.claude/skills/` in `--add-dir` directories auto-load
- Fixed: `@` file completion showing incorrect relative paths from subdirectories
- Fixed: Bash tool no longer throws "Bad substitution" errors with JavaScript template literals (e.g., `${index + 1}`)
- Improved: Skill character budget now scales with context window (2% of context)
- Improved: `--resume` re-uses `--agent` value from previous conversation by default
- Fixed: Thai/Lao spacing vowels rendering issues
- [VSCode] Fixed slash commands incorrectly executing when pressing Enter with preceding text
- [VSCode] Added spinner when loading past conversations list

### v2.1.31 (2026-02-03)

- **Session resume hint** — Exit message now shows how to continue your conversation later
- **Full-width (zenkaku) space support** — Added Japanese IME checkbox selection support
- Fixed: PDF too large errors permanently locking sessions (now recoverable without starting new conversation)
- Fixed: Bash commands incorrectly reporting "Read-only file system" when sandbox enabled
- Fixed: Plan mode crash when project config missing default fields
- Fixed: `temperatureOverride` being silently ignored in streaming API path
- Fixed: LSP shutdown/exit compatibility with strict language servers
- Improved: System prompts now guide model toward Read/Edit/Glob/Grep tools instead of bash equivalents
- Improved: PDF and request size error messages show actual limits (100 pages, 20MB)
- Reduced: Layout jitter when spinner appears/disappears during streaming

### v2.1.30 (2026-02-02)

- **⭐ PDF page range support** — `pages` parameter in Read tool for PDFs (e.g., `pages: "1-5"`) with lightweight references for large PDFs (>10 pages)
- **⭐ Pre-configured OAuth for MCP servers** — Built-in client credentials for servers without Dynamic Client Registration (Slack support via `--client-id` and `--client-secret`)
- **⭐ New `/debug` command** — Claude can help troubleshoot current session issues
- **Additional git flags** — Support for `git log` and `git show` read-only flags (`--topo-order`, `--cherry-pick`, `--format`, `--raw`)
- **Task tool metrics** — Results now include token count, tool uses, and duration
- **Reduced motion mode** — New config option for accessibility
- Fixed: Phantom "(no content)" text blocks in API history (reduces token waste)
- Fixed: Prompt cache not invalidating when tool schemas changed
- Fixed: 400 errors after `/login` with thinking blocks
- Fixed: Session resume hang with corrupted `parentUuid` cycles
- Fixed: Rate limit showing wrong "/upgrade" for Max 20x users
- Fixed: Permission dialogs stealing focus while typing
- Fixed: Subagents unable to access SDK MCP tools
- Fixed: Windows users with `.bashrc` unable to run bash
- Improved: Memory usage for `--resume` (68% reduction for many sessions)
- Improved: TaskStop displays stopped command description instead of generic message
- Changed: `/model` executes immediately instead of queuing
- [VSCode] Added multiline input in "Other" text fields (Shift+Enter for new lines)
- [VSCode] Fixed duplicate sessions in session list

### v2.1.29 (2026-01-31)

- **Performance**: Fixed startup performance issues when resuming sessions with saved hook context
- Significantly improved session recovery speed for long-duration sessions

### v2.1.27 (2026-01-29)

- **New**: `--from-pr` flag to resume sessions linked to a specific GitHub PR number or URL
- **New**: Sessions automatically linked to PRs when created via `gh pr create`
- Added tool call failures and denials to debug logs
- Fixed context management validation error for Bedrock/Vertex gateway users
- Fixed `/context` command not displaying colored output
- Fixed status bar duplicating background task indicator when PR status was shown
- [Windows] Fixed bash command execution failing for users with `.bashrc` files
- [Windows] Fixed console windows flashing when spawning child processes
- [VSCode] Fixed OAuth token expiration causing 401 errors after extended sessions

### v2.1.25 (2026-01-30)

- Fixed beta header validation for Bedrock and Vertex gateway users — Ensures `CLAUDE_CODE_DISABLE_EXPERIMENTAL_BETAS=1` environment variable works correctly

### v2.1.23 (2026-01-29)

- **Customizable spinner verbs** — New `spinnerVerbs` setting allows personalization of spinner action words
- mTLS and corporate proxy connectivity fixes — Improved support for users behind corporate proxies with client certificates
- Per-user temp directory isolation — Prevents permission conflicts on shared systems
- Improved terminal rendering performance — Optimized screen data layout for faster updates
- Fixed: Prompt caching race condition causing 400 errors
- Fixed: Async hooks not canceling when headless streaming ends
- Fixed: Tab completion not updating input field
- Fixed: Ripgrep search timeouts returning empty results instead of errors
- Changed: Bash commands show timeout duration alongside elapsed time
- Changed: Merged PRs show purple status indicator in prompt footer
- [IDE] Fixed: Model options displaying incorrect region strings for Bedrock users in headless mode

### v2.1.22 (2026-01-28)

- Improved task UI performance with virtualization — Task list now uses virtual scrolling for better responsiveness with many tasks
- Vim selection and deletion fixes — Fixed visual mode selections and `dw` command behavior
- LSP improvements: Kotlin support, UTF-16 range handling, better error recovery
- Tasks now consistently use `task-N` IDs instead of internal UUIDs
- Fixed: `#` keyboard shortcut not working in task creation fields
- Fixed: Compact tool use rendering in chat history
- Fixed: Session URL escaping in git commit messages
- Fixed: Command output handling improvements

### v2.1.21 (2026-01-28)

- **Skills/commands can specify required/recommended Claude Code version** — Use `minClaudeCodeVersion` and `recommendedClaudeCodeVersion` in frontmatter
- **New TaskCreate fields**: `category` (testing, implementation, documentation, etc.), `checklist` (subtasks as markdown list), `parentId` (task hierarchy)
- **Automatic Claude Code update checking** at session start (respects auto-update settings)
- Tasks appear in `/context` output with 'Disable tasks' shortcut for quick toggling
- Improved task UI: Delete button added, better empty state messaging
- Fixed: Task deletion now properly removes all related task data
- Fixed: Shell environment variables expanded correctly in hook commands
- Fixed: Pasted URLs with parentheses properly formatted in markdown
- Fixed: Bash output capture for commands with large output

### v2.1.20 (2026-01-27)

- **New**: TaskUpdate tool can delete tasks via `status="deleted"`
- **New**: PR review status indicator in prompt footer — Shows PR state (approved, changes requested, pending, draft) as colored dot with clickable link
- Arrow key history navigation in vim normal mode when cursor cannot move further
- External editor shortcut (Ctrl+G) added to help menu
- Support for loading CLAUDE.md from `--add-dir` directories (requires `CLAUDE_CODE_ADDITIONAL_DIRECTORIES_CLAUDE_MD=1`)
- Fixed: Session compaction issues causing full history load instead of compact summary
- Fixed: Agents ignoring user messages while actively working
- Fixed: Wide character (emoji, CJK) rendering artifacts
- Improved: Task list dynamically adjusts to terminal height
- Changed: Background agents prompt for tool permissions before launching
- Changed: Config backups timestamped and rotated (keeps 5 most recent)

### v2.1.19 (2026-01-25)

- **New**: `CLAUDE_CODE_ENABLE_TASKS` environment variable — Set to `false` to temporarily revert to old task system
- **New**: Argument shorthand in custom commands — Use `$0`, `$1`, etc. instead of verbose syntax
- [VSCode] Session forking and rewind functionality enabled for all users
- Fixed: Crashes on processors without AVX instruction support
- Fixed: Dangling Claude Code processes when terminal closed (SIGKILL fallback)
- Fixed: `/rename` and `/tag` not updating correct session when resuming from different directory (git worktrees)
- Fixed: Resuming sessions by custom title from different directory
- Fixed: Pasted text lost when using prompt stash (Ctrl+S) and restore
- Fixed: Agent list displaying "Sonnet (default)" instead of "Inherit (default)" for agents without explicit model
- Fixed: Backgrounded hook commands blocking session instead of returning early
- Fixed: File write preview omitting empty lines
- Changed: Skills without additional permissions/hooks allowed without approval
- [SDK] Added replay of queued_command attachment messages when `replayUserMessages` enabled

**⚠️ Breaking**:
- Indexed argument syntax changed: `$ARGUMENTS.0` → `$ARGUMENTS[0]` (bracket syntax)

### v2.1.18 (2026-01-24) ⭐

- ⭐ **Customizable keyboard shortcuts** — Configure keybindings per context, create chord sequences, personalize workflow
- Run `/keybindings` to get started
- Learn more: [code.claude.com/docs/en/keybindings](https://code.claude.com/docs/en/keybindings)

### v2.1.17 (2026-01-23)

- Fix: Crashes on processors without AVX instruction support

### v2.1.16 (2026-01-22) ⭐

- ⭐ **New task management system** with dependency tracking
- [VSCode] Native plugin management support
- [VSCode] OAuth users can browse and resume remote sessions from Sessions dialog
- Fixed: Out-of-memory crashes when resuming sessions with heavy subagent usage
- Fixed: "Context remaining" warning not hidden after `/compact`
- [IDE] Fixed race condition on Windows where sidebar view container wouldn't appear

### v2.1.15 (2026-01-22)

- **⚠️ Deprecation notice for npm installations** — Run `claude install` or see [docs](https://docs.anthropic.com/en/docs/claude-code/getting-started)
- Improved UI rendering performance with React Compiler
- Fixed: MCP stdio server timeout not killing child process, which could cause UI freezes

### v2.1.14 (2026-01-21)

- **History-based autocomplete in bash mode** — Type `!` followed by a partial command and press Tab to complete from bash history
- Search functionality in installed plugins list
- Support for pinning plugins to specific git commit SHAs for exact version control
- Fixed: Context window blocking limit calculated too aggressively (~65% instead of ~98%)
- Fixed: Memory issues and leaks in long-running sessions with parallel subagents
- Fixed: `@` symbol incorrectly triggering file autocomplete in bash mode
- Fixed: Slash command autocomplete selecting wrong command for similar names
- Improved: Backspace deletes pasted text as single token

### v2.1.12 (2026-01-18)

- Bug fix: Message rendering

### v2.1.11 (2026-01-17)

- Fix: Excessive MCP connection requests for HTTP/SSE transports

### v2.1.10 (2026-01-17)

- New `Setup` hook event (--init, --init-only, --maintenance flags)
- Keyboard shortcut 'c' to copy OAuth URL
- File suggestions show as removable attachments
- [VSCode] Plugin install count + trust warnings

### v2.1.9 (2026-01-16)

- **`auto:N` syntax for MCP tool search threshold** — Configure when Tool Search activates: `ENABLE_TOOL_SEARCH=auto:5` (5% context), `auto:10` (default), `auto:20` (conservative). See [architecture.md](./architecture.md#mcp-tool-search-lazy-loading) for details.
- `plansDirectory` setting for custom plan file locations
- Session URL attribution to commits/PRs from web sessions
- PreToolUse hooks can return `additionalContext`
- `${CLAUDE_SESSION_ID}` string substitution for skills

### v2.1.7 (2026-01-15)

- `showTurnDuration` setting to hide turn duration messages
- **MCP Tool Search auto mode enabled by default** — Lazy loading for MCP tools when definitions exceed 10% of context. Based on Anthropic's [Advanced Tool Use](https://www.anthropic.com/engineering/advanced-tool-use) API feature. Result: **85% token reduction** on tool definitions, improved tool selection accuracy (Opus 4: 49%→74%, Opus 4.5: 79.5%→88.1%)
- Inline display of agent final response in task notifications

**⚠️ Breaking**:
- OAuth/API Console URLs changed: `console.anthropic.com` → `platform.claude.com`
- Security fix: Wildcard permission rules could match compound commands

### v2.1.6 (2026-01-14)

- Search functionality in `/config` command
- Date range filtering in `/stats` (press `r` to cycle)
- Auto-discovery of skills from nested `.claude/skills` directories
- Updates section in `/doctor` showing auto-update channel

**⚠️ Security Fix**: Permission bypass via shell line continuation

### v2.1.5 (2026-01-13)

- `CLAUDE_CODE_TMPDIR` environment variable for custom temp directory

### v2.1.4 (2026-01-12)

- `CLAUDE_CODE_DISABLE_BACKGROUND_TASKS` environment variable

### v2.1.3 (2026-01-11)

- Merged slash commands and skills (simplified mental model)
- Release channel toggle (stable/latest) in `/config`
- `/doctor` warnings for unreachable permission rules

### v2.1.2 (2026-01-10)

- Windows Package Manager (winget) support
- Clickable hyperlinks for file paths (OSC 8 terminals)
- Shift+Tab shortcut in plan mode for auto-accept edits
- Large bash outputs saved to disk instead of truncated

**⚠️ Breaking**:
- Security fix: Command injection in bash command processing
- Deprecated: `C:\ProgramData\ClaudeCode` managed settings path

### v2.1.0 (2026-01-08) ⭐ MAJOR

**Highlights**:
- ⭐ **Automatic skill hot-reload** — Skills modified in `~/.claude/skills` or `.claude/skills` immediately available
- ⭐ **Shift+Enter works OOTB** in iTerm2, WezTerm, Ghostty, Kitty
- ⭐ **New Vim motions**: `;` `,` `y` `yy` `Y` `p` `P` text objects (`iw` `aw` `i"` etc.) `>>` `<<` `J`
- **Unified Ctrl+B** for backgrounding all running tasks
- `/plan` command shortcut to enable plan mode
- Slash command autocomplete anywhere in input
- `language` setting for response language (e.g., `language: "japanese"`)
- Skills `context: fork` support for forked sub-agent context
- Hooks support in agent/skill/command frontmatter
- MCP `list_changed` notifications support
- `/teleport` and `/remote-env` commands for web sessions
- Disable specific agents with `Task(AgentName)` syntax
- `--tools` flag in interactive mode
- YAML-style lists in frontmatter `allowed-tools`

**⚠️ Breaking**:
- OAuth URLs: `console.anthropic.com` → `platform.claude.com`
- Removed permission prompt for entering plan mode
- [SDK] Minimum zod peer dependency: `^4.0.0`

---

## 2.0.x Series (November 2025 - January 2026)

### v2.0.76 (2026-01-05)

- Fix: macOS code-sign warning with Claude in Chrome

### v2.0.74 (2026-01-04) ⭐

- ⭐ **LSP (Language Server Protocol) tool** for code intelligence (go-to-definition, find references, hover)
- `/terminal-setup` for Kitty, Alacritty, Zed, Warp
- Ctrl+T in `/theme` to toggle syntax highlighting
- Grouped skills/agents by source in `/context`

### v2.0.72 (2026-01-02) ⭐

- ⭐ **Claude in Chrome (Beta)** — Control browser directly from Claude Code
- Reduced terminal flickering
- QR code for mobile app download
- Thinking toggle changed: Tab → Alt+T

### v2.0.70 (2025-12-30)

- Enter key accepts/submits prompt suggestions immediately
- Wildcard syntax `mcp__server__*` for MCP tool permissions
- Auto-update toggle for plugin marketplaces
- 3x memory usage improvement for large conversations

**⚠️ Breaking**: Removed `#` shortcut for quick memory entry

### v2.0.67 (2025-12-26) ⭐

- ⭐ **Thinking mode enabled by default for Opus 4.5**
- Thinking config moved to `/config`
- Search in `/permissions` with `/` shortcut

### v2.0.64 (2025-12-22) ⭐

- ⭐ **Instant auto-compacting**
- ⭐ **Async agents and bash commands** with wake-up messages
- `/stats` with usage graphs, streaks, favorite model
- Named sessions: `/rename`, `/resume <name>`
- Support for `.claude/rules/` directory
- Image dimension metadata for coordinate mappings

### v2.0.60 (2025-12-18) ⭐

- ⭐ **Background agents** — Agents run while you work
- `--disable-slash-commands` CLI flag
- Model name in Co-Authored-By commits
- `/mcp enable|disable [server-name]`

### v2.0.51 (2025-12-10) ⭐ MAJOR

- ⭐ **Opus 4.5 released**
- ⭐ **Claude Code for Desktop**
- Updated usage limits for Opus 4.5
- Plan Mode builds more precise plans

### v2.0.45 (2025-12-05) ⭐

- ⭐ **Microsoft Foundry support**
- `PermissionRequest` hook for auto-approve/deny
- `&` prefix for background tasks to web

### v2.0.28 (2025-11-18) ⭐

- ⭐ **Plan mode: introduced Plan subagent**
- Subagents: resume capability
- Subagents: dynamic model selection
- `--max-budget-usd` flag (SDK)
- Git-based plugins branch/tag support (`#branch`)

### v2.0.24 (2025-11-10)

- Claude Code Web: Web → CLI teleport
- Sandbox mode for BashTool (Linux & Mac)
- Bedrock: `awsAuthRefresh` output display

---

## Breaking Changes Summary

### URLs

| Version | Change |
|---------|--------|
| v2.1.0, v2.1.7 | OAuth/API Console: `console.anthropic.com` → `platform.claude.com` |

### Windows

| Version | Change |
|---------|--------|
| v2.0.58 | Managed settings prefer `C:\Program Files\ClaudeCode` |
| v2.1.2 | Deprecated `C:\ProgramData\ClaudeCode` path |

### SDK

| Version | Change |
|---------|--------|
| v2.0.25 | Removed legacy SDK entrypoint → `@anthropic-ai/claude-agent-sdk` |
| v2.1.0 | Minimum zod peer dependency: `^4.0.0` |

### API Ecosystem

| Date | Feature |
|------|---------|
| 2026-01-29 | **Structured Outputs GA**: `output_config.format` remplace `output_format`. [Docs](https://platform.claude.com/docs/en/build-with-claude/structured-outputs) |

### Shortcuts

| Version | Change |
|---------|--------|
| v2.0.70 | Removed `#` shortcut for quick memory entry |

### Security Fixes

| Version | Issue |
|---------|-------|
| v2.1.2 | Command injection in bash command processing |
| v2.1.6 | Shell line continuation permission bypass |
| v2.1.7 | Wildcard permission rules compound commands |
| v2.1.38 | Heredoc delimiter command smuggling prevention |

### Syntax

| Version | Change |
|---------|--------|
| v2.1.19 | Indexed argument syntax changed: `$ARGUMENTS.0` → `$ARGUMENTS[0]` (bracket syntax) |

---

## Milestone Features

| Version | Key Features |
|---------|--------------|
| **v2.1.32** | Opus 4.6, Agent teams preview, Automatic memory |
| **v2.1.18** | Customizable keyboard shortcuts with /keybindings |
| **v2.1.16** | New task management system with dependency tracking |
| **v2.1.0** | Skill hot-reload, Shift+Enter OOTB, Vim motions, /plan command |
| **v2.0.74** | LSP tool for code intelligence |
| **v2.0.72** | Claude in Chrome (browser control) |
| **v2.0.67** | Thinking mode default for Opus 4.5 |
| **v2.0.64** | Instant auto-compact, async agents, named sessions |
| **v2.0.60** | Background agents |
| **v2.0.51** | Opus 4.5, Claude Code for Desktop |
| **v2.0.45** | Microsoft Foundry, PermissionRequest hook |
| **v2.0.28** | Plan subagent, subagent resume/model selection |
| **v2.0.24** | Web teleport, Sandbox mode |

---

## Updating This Document

1. **Watch**: [github.com/anthropics/claude-code/releases](https://github.com/anthropics/claude-code/releases)
2. **Update**: `machine-readable/claude-code-releases.yaml` (source of truth)
3. **Regenerate**: Update this markdown accordingly
4. **Sync landing**: Run `./scripts/check-landing-sync.sh`

---

*Last updated: 2026-02-13 | [Back to main guide](./ultimate-guide.md)*
