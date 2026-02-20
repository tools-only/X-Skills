# Kaizen Transcript Analysis: Data Schema & Feasibility Report

**Date:** 2026-02-18
**Analyst:** kaizen-data-analyst subagent
**Dataset:** `-home-ubuntulinuxqa2-repos-claude-skills` project (richest project)
**Sessions sampled:** 10 sessions (tiny/small/medium/large) + corpus-wide scans

---

## 1. Data Schema Reference

### 1.1 Corpus Overview

```
Location:    ~/.claude/projects/-home-ubuntulinuxqa2-repos-claude-skills/
Total files: 720 JSONL files in the main project directory
  - 614 main session JSONL files (named {uuid}.jsonl)
  - 106 top-level agent JSONL files (named agent-{id}.jsonl, orphan agents)
  - 83 session subdirectories (named {uuid}/)
    - 79 of those contain a `subagents/` sub-directory
    - Total subagent JSONL files inside subdirs: 687
    - Also contains `tool-results/` dir with individual .txt files (one per async task)
Versions seen: 2.0.50, 2.1.14–2.1.45 (progressive)
Date range: 2026-01-09 to 2026-02-18 (most activity on 2026-01-09, 2026-01-13, 2026-01-15)
Sessions with meaningful tool use: ~88 (rest are <5 lines or billing-error sessions)
```

### 1.2 Record Types

Every line in a JSONL file is one JSON object. The `type` field is the primary discriminator.

**Corpus-wide type counts:**

```
progress              21,084  (hook progress + subagent streaming)
assistant             15,682  (LLM response turns)
user                  11,772  (human input + tool results)
queue-operation        4,609  (conversation summary queue)
file-history-snapshot  1,944  (file edit tracking)
system.stop_hook_summary 920  (per-turn hook feedback)
system.turn_duration     600  (timing metadata)
summary                  524  (session title/summary)
system.compact_boundary   84  (context compaction events)
system.local_command      20  (user slash commands like /clear, /skills)
system.api_error          14  (connection failures + retries)
system.microcompact_boundary 3
saved_hook_context         1  (hook-injected context persisted to transcript)
```

### 1.3 Common Fields (All Record Types)

Every record with session context carries:

```json
{
  "type": "...",
  "uuid": "<message-uuid>",
  "parentUuid": "<parent-uuid | null>",
  "timestamp": "2026-01-30T16:39:40.005Z",
  "sessionId": "<session-uuid>",
  "cwd": "~/my/current/dir",
  "version": "2.1.25",
  "gitBranch": "main",
  "isSidechain": false,
  "userType": "external",
  "slug": "squishy-discovering-breeze"
}
```

**Subagent-specific additions:**

```json
{
  "agentId": "a0ec468",
  "isSidechain": true,
  "agentName": "external-researcher",  // only when named agent in team
  "teamName": "kaizen-scoping"         // only when part of a team
}
```

### 1.4 Record Type Schemas

#### `user` — Human or Tool Result

```json
{
  "type": "user",
  "message": {
    "role": "user",
    "content": "<string> | [{type, text}] | [{type: tool_result, tool_use_id, content, is_error}]"
  },
  "toolUseResult": true,          // present when this is a tool result
  "sourceToolAssistantUUID": "...", // UUID of the assistant turn that called the tool
  "uuid": "...",
  "parentUuid": "..."
}
```

**Tool result content structure:**

```json
{
  "type": "tool_result",
  "tool_use_id": "toolu_01...",
  "content": [{"type": "text", "text": "..."}] | "<string>",
  "is_error": false  // true for errors
}
```

#### `assistant` — LLM Response Turn

```json
{
  "type": "assistant",
  "requestId": "req_011CYEpNE8U6bRtNVbuJBjaW",
  "message": {
    "id": "msg_018oRhgn7XnKrXg9DKBiKsvS",
    "model": "claude-sonnet-4-6",
    "role": "assistant",
    "type": "message",
    "stop_reason": "tool_use | end_turn | stop_sequence | null",
    "stop_sequence": null,
    "content": [
      {"type": "text", "text": "..."},
      {
        "type": "tool_use",
        "id": "toolu_01...",
        "name": "Bash",
        "input": {"command": "..."},
        "caller": {"type": "direct"}
      }
    ],
    "usage": {
      "input_tokens": 1234,
      "output_tokens": 567,
      "cache_creation_input_tokens": 49394,
      "cache_read_input_tokens": 0,
      "cache_creation": {
        "ephemeral_5m_input_tokens": 0,
        "ephemeral_1h_input_tokens": 49394
      },
      "service_tier": "standard"
    }
  },
  "error": "billing_error"  // present only on error turns
}
```

#### `progress` — Background Processing Events

Two subtypes within `data.type`:

**`hook_progress`** — Hook execution notification:

```json
{
  "type": "progress",
  "data": {
    "type": "hook_progress",
    "hookEvent": "SessionStart",
    "hookName": "SessionStart:startup",
    "command": "node \"~/.claude/hooks/gsd-check-update.js\""
  },
  "parentToolUseID": "...",
  "toolUseID": "..."
}
```

**`agent_progress`** — Subagent streaming (full inner turn):

```json
{
  "type": "progress",
  "data": {
    "type": "agent_progress",
    "message": {
      "type": "user | assistant",
      "timestamp": "...",
      "message": { /* full user or assistant message object */ }
    }
  }
}
```

#### `system` — Metadata Events

**`stop_hook_summary`** — After every turn:

```json
{
  "type": "system",
  "subtype": "stop_hook_summary",
  "hookCount": 3,
  "hookInfos": [{"command": "node ..."}],
  "hookErrors": [],
  "preventedContinuation": false,
  "stopReason": "",
  "hasOutput": true,
  "level": "suggestion",
  "toolUseID": "..."
}
```

**`turn_duration`** — Turn timing:

```json
{
  "type": "system",
  "subtype": "turn_duration",
  "durationMs": 12500,
  "isMeta": true
}
```

**`api_error`** — API failures:

```json
{
  "type": "system",
  "subtype": "api_error",
  "level": "error",
  "cause": {"code": "ConnectionRefused", "path": "https://api.anthropic.com/...", "errno": 0},
  "retryInMs": 595.8,
  "retryAttempt": 1,
  "maxRetries": 10
}
```

**`compact_boundary`** — Context compaction:

```json
{
  "type": "system",
  "subtype": "compact_boundary",
  "content": "Conversation compacted",
  "logicalParentUuid": "...",
  "isMeta": false
}
```

**`local_command`** — User slash commands:

```json
{
  "type": "system",
  "subtype": "local_command",
  "content": "<command-name>/clear</command-name>\n<command-message>clear</command-message>",
  "level": "info"
}
```

#### `file-history-snapshot` — File Edit Tracking

```json
{
  "type": "file-history-snapshot",
  "messageId": "...",
  "snapshot": {
    "messageId": "...",
    "trackedFileBackups": {
      ".claude/skills/rt-ica/SKILL.md": {
        "backupFileName": "fce4cdb71d269c08@v1",
        "version": 1,
        "backupTime": "2026-01-27T18:04:22.616Z"
      }
    },
    "timestamp": "..."
  },
  "isSnapshotUpdate": false
}
```

#### `queue-operation` — Summary Queue

```json
{
  "type": "queue-operation",
  "operation": "enqueue | dequeue",
  "timestamp": "...",
  "sessionId": "...",
  "content": "[string or json payload]"
}
```

#### `summary` — Session Summaries

```json
{
  "type": "summary",
  "summary": "Credit balance insufficient error loop",
  "leafUuid": "..."
}
```

### 1.5 Tool Call Schema

**Tool use (inside `assistant.message.content[]`):**

```json
{
  "type": "tool_use",
  "id": "toolu_01...",
  "name": "Bash",
  "input": {
    "command": "git status",
    "description": "Show working tree status"
  },
  "caller": {"type": "direct"}
}
```

**Tool inputs by tool name:**

```
Bash:        {"command": "...", "description": "..."}
Read:        {"file_path": "...", "limit": N, "offset": N}
Edit:        {"file_path": "...", "old_string": "...", "new_string": "..."}
Write:       {"file_path": "...", "content": "..."}
Grep:        {"pattern": "...", "path": "...", "output_mode": "content|files_with_matches"}
Glob:        {"pattern": "...", "path": "..."}
Task:        {"description": "...", "subagent_type": "...", "prompt": "...",
               "run_in_background": true, "team_name": "...", "name": "...", "model": "..."}
TaskUpdate:  {"taskId": "...", "status": "...", "owner": "..."}
SendMessage: {"type": "message", "recipient": "...", "content": "...", "summary": "..."}
Skill:       {"skill": "...", "args": "..."}
AskUserQuestion: {"questions": [...]}
```

**Task tool result (for subagents):**

```json
{
  "type": "tool_result",
  "tool_use_id": "toolu_01...",
  "content": [{"type": "text",
    "text": "Async agent launched successfully.\nagentId: a3e2ae9 ..."}]
}
```

The `agentId` in the result links to `{session-dir}/subagents/agent-{agentId}.jsonl`.

### 1.6 Subagent Relationship Model

```
Main session JSONL (orchestrator):
  └─ Task tool_use (id: toolu_01X)
       └─ tool_result → "agentId: abc1234"
            └─ Subagent at: {session-uuid}/subagents/agent-abc1234.jsonl
                 └─ Records with:
                      - isSidechain: true
                      - sessionId: {parent session uuid}
                      - agentId: abc1234
                      - First record is user message = subagent's task prompt
```

Subagent final output is returned via the parent session's tool_result when the Task completes (synchronous) or when `TaskOutput` is polled (asynchronous, `run_in_background: true`).

**Async agent flow:**

1. Orchestrator calls `Task` with `run_in_background: true`
2. Tool result: "Async agent launched successfully. agentId: abc. output_file: /tmp/..."
3. Orchestrator calls `TaskOutput` to poll
4. Eventually `TaskOutput` returns the full subagent output text
5. Also stored in `{session-dir}/tool-results/{task-tool-use-id}.txt`

---

## 2. Signal Catalog

### 2.1 Tool Misuse Detection

**Signal:** `Bash` tool used for `ls`, `cat`, `grep`, `find`, `head`, `tail`, `sed`, `awk` when built-in tools (`Glob`, `Read`, `Grep`) should have been used.

**Data field path:** `assistant.message.content[].input.command` where `assistant.message.content[].name == "Bash"`

**Extraction:** Parse the `command` string for these patterns:

```python
patterns = {
    'should_use_grep': r'\bgrep\b(?!.*\|.*grep)',   # standalone grep
    'should_use_glob': r'\bfind\b\s+\S+\s+-name\b', # find -name
    'should_use_read': r'\bcat\b\s+\S+\.(?:md|py|js|json|yaml)',
    'should_use_read_head': r'\bhead\b\s+-\d+\s+\S+\.\w+',
    'should_use_read_tail': r'\btail\b\s+-\d+\s+\S+\.\w+',
    'should_use_ls': r'^\s*ls\b',  # ls as primary command
}
```

**Corpus-wide counts (this project):**

```
grep:  182 occurrences
ls:    215 occurrences
head:  105 occurrences
find:   42 occurrences
cat:    27 occurrences
tail:   16 occurrences
sed:     6 occurrences
```

**Extraction difficulty:** Trivial — direct string/regex match on `input.command`.

**Caveats:** Must exclude legitimate uses: `git ... | grep`, `uv run ... | head`, pipeline contexts where piping to grep is appropriate. The `description` field on Bash calls sometimes signals intent.

**Example:**

```json
{
  "name": "Bash",
  "input": {
    "command": "grep -r \"Hook Development\" .claude/ ~/.claude/ 2>/dev/null | head -20",
    "description": "Search for hook development references"
  }
}
```

→ This violates the `Grep` tool rule AND pipes to `head` (should use `Grep` with `head_limit`).

---

### 2.2 Error Patterns

**Signal:** Tool result with `is_error: true` or specific error text patterns.

**Data field path:** `user.message.content[].is_error` (where `type == "tool_result"`)

**Error types found:**

```
[6]  Request interrupted by user for tool use
[4]  File has not been read yet (Edit before Read)
[3]  Exit code 1 - pre-commit hook failures
[2]  User denied tool use
[2]  Exit code 2 - file type errors
[1]  String to replace not found in file (Edit failure)
[1]  Unknown skill
[1]  File has been modified since read (stale Edit)
[1]  File content exceeds token limit
[1]  command not found (missing binary)
[1]  File does not exist
[14] billing_error (on assistant records, not tool results)
[587] billing_error sessions (entire sessions)
```

**Extraction difficulty:** Trivial for structured errors, moderate for classifying Bash exit codes.

**Richest pattern for kaizen:** "File has not been read yet" and "String to replace not found" indicate the agent tried to edit without a preceding Read — a trackable anti-pattern.

---

### 2.3 User Frustration Signals

**Signal:** User messages containing correction/frustration language outside of system-generated context.

**Data field path:** `user.message.content[text]` where `user.toolUseResult` is absent/false.

**Patterns to match:**

```python
frustration_signals = [
    r'\[Request interrupted by user\]',        # User hit Ctrl+C
    r'\[Request interrupted by user for tool use\]',  # Denied tool
    r"^(No,|Don't|Stop |Why did you|Why would)",
    r'(wrong|incorrect|that\'s not|you should not)',
]
```

**Corpus stats (88 substantive sessions):**

- User interrupts: 106 occurrences
- Direct corrections ("No,", "Why did you", "wrong"): ~43 real messages

**Key challenge:** Must filter out:

- System-generated messages with `Context: This summary...`
- Session continuation messages (`This session is being continued...`)
- Stop hook feedback messages (`Stop hook feedback:...`)
- Skill injection messages (`Base directory for this skill:...`)

**Real examples found:**

```
"WHy would we ever want to use additional scripts. especially a bash script when this repo is meant to run cross platform?"
"why did your agent planner...end up specifying `command -v claude`...?"
"This means you have not validated your plan."
"No, we need to remove the commands skill..."
"hey, use single line descriptions."
```

**Extraction difficulty:** Moderate — requires filtering system messages. The key discriminator is: real user messages are short (<300 chars typically), directly responsive, and lack XML tags.

---

### 2.4 Workflow Repetition / Common Sequences

**Signal:** Ordered sequences of tool names across sessions, compared to find common patterns.

**Data field path:** `assistant.message.content[].name` (tool_use blocks), in message order.

**Most common 3-tool trigrams (this project):**

```
1066: Bash → Bash → Bash
700:  Read → Read → Read
238:  Edit → Edit → Edit
162:  Bash → Bash → Read
140:  Read → Edit → Edit
137:  SendMessage → SendMessage → SendMessage
136:  Grep → Grep → Grep
119:  Task → Task → Task
114:  Grep → Read → Read
113:  TaskList → TaskList → TaskList
```

**Most common tool usage across all sessions:**

```
Read:       2,248
Bash:       2,245
Edit:         939
Grep:         721
Task:         410
TaskUpdate:   366
Glob:         352
SendMessage:  341
TaskList:     227
Write:        211
```

**Observation:** `Bash` and `Read` are nearly equal — suggesting Bash is used in place of Read in many cases. `Bash` should be ~10-20% of `Read`+`Grep`+`Glob` combined in a well-disciplined workflow.

**Workflow fingerprinting feasibility:** Session-level tool sequences can be extracted and compared with:

- Edit distance (Levenshtein on tool-name arrays)
- Subsequence mining (PrefixSpan)
- TF-IDF on tool-sequence bigrams

**Extraction difficulty:** Trivial to extract sequences; moderate to compare across sessions (needs Python `mlxtend` or `pm4py`).

---

### 2.5 Subagent Delegation Patterns

**Signal:** `Task` tool_use blocks with `subagent_type`, `prompt`, `run_in_background`.

**Data field path:** `assistant.message.content[].input` where `name == "Task"`

**Full input schema:**

```json
{
  "description": "Research LLM observability platforms",
  "subagent_type": "comprehensive-researcher",
  "prompt": "...",
  "run_in_background": true,
  "team_name": "kaizen-scoping",
  "name": "external-researcher",
  "model": "sonnet",
  "max_turns": 10
}
```

**Subagent types used (all sessions, ranked):**

```
127: general-purpose
 57: python3-development:python-cli-architect
 53: claude-context-optimizer
 45: Explore
 15: context-gathering
 11: comprehensive-researcher
 11: plugin-assessor
  9: Plan
  8: python3-development:python-pytest-architect
  7: Bash
  7: holistic-linting:linting-root-cause-resolver
  6: contextual-ai-documentation-optimizer
  6: python3-development:swarm-task-planner
  6: subagent-refactorer
```

**Extraction difficulty:** Trivial. Subagent type, description, and prompt snippet are all directly extractable.

**For kaizen:** Tracking when `general-purpose` is used vs. a more specialized agent (when one existed) reveals delegation anti-patterns.

---

### 2.6 Context Waste Detection

**Signal:** Orchestrator reads a file then immediately delegates a Task where the same file path appears in the Task's prompt.

**Data field paths:**

- `Read` input: `assistant.message.content[].input.file_path` where `name == "Read"`
- `Task` prompt: `assistant.message.content[].input.prompt` where `name == "Task"`

**Detection algorithm:**

```python
# Within a session, track (tool_name, file_path) tuples in order
# Flag: Read(X) followed within 5 turns by Task(prompt contains X)
```

**Corpus result:** 16 instances of context waste in 88 substantive sessions (18% of sessions).

The most common wasted file: `.claude/CLAUDE.md` — orchestrator reads the project instructions then includes the path in a subagent prompt that will read it again.

**Extraction difficulty:** Moderate — requires cross-referencing Read inputs with Task prompts in a sliding window per session.

---

## 3. Cross-Session Analysis Feasibility

### 3.1 What's Possible with Simple Aggregation

These analyses require only regex/string matching on the JSONL corpus:

| Analysis                     | Method                                                 | Output            |
| ---------------------------- | ------------------------------------------------------ | ----------------- |
| Tool misuse rate per session | Count Bash with file-op patterns ÷ total Bash          | Per-session score |
| Error rate per session       | Count `is_error: true` tool results ÷ total tool calls | Per-session score |
| User interrupt rate          | Count interrupt messages ÷ total user turns            | Frustration proxy |
| Billing error sessions       | Scan for `error: "billing_error"`                      | Blocked sessions  |
| Subagent type distribution   | Count `Task.input.subagent_type`                       | Agent usage map   |
| Context waste occurrences    | Read→Task overlap within 5 turns                       | Waste instances   |
| Session length distribution  | Line count of JSONL files                              | Complexity proxy  |
| Token usage per session      | Sum `usage.input_tokens + output_tokens`               | Cost analysis     |
| Cache efficiency             | `cache_read_input_tokens ÷ total_input_tokens`         | Cache hit rate    |

### 3.2 What Requires Deeper Techniques

| Analysis                          | Why Harder                                                 | Approach                          |
| --------------------------------- | ---------------------------------------------------------- | --------------------------------- |
| Workflow clustering               | Need to compare sequences across 720 sessions              | Sequence similarity + k-means     |
| "Same mistake" detection          | Need semantic similarity of error messages                 | Embedding + clustering            |
| Instruction drift detection       | Match user corrections to prior assistant output           | Context-window pair analysis      |
| Inefficiency pattern mining       | Extract frequent tool sequences that correlate with errors | PrefixSpan + rule mining          |
| Missing tool opportunities        | Detect bash pipelines that match Glob/Grep/Read semantics  | AST parsing of bash cmds          |
| User frustration severity scoring | Classify correction messages by severity                   | LLM classifier or keyword scoring |

---

## 4. Recommended Data Pipeline Architecture

### 4.1 Processing Strategy

With 720 sessions and ~57k total records, full corpus processing is fast in Python:

- Estimated corpus size: ~500MB total
- Single-pass scan with Counter: ~30 seconds
- Session-level analysis with dict building: ~2 minutes

**No streaming or distributed processing needed** at this scale.

### 4.2 Data Model

```python
# Core data model for efficient analysis

Session = {
    "session_id": str,       # UUID
    "is_subagent": bool,
    "parent_session_id": str | None,
    "agent_id": str | None,
    "timestamps": {"start": str, "end": str},
    "version": str,
    "git_branch": str,
    "slug": str,
    "team_name": str | None,
    "agent_name": str | None,

    # Derived metrics
    "tool_calls": [{"name": str, "input": dict, "timestamp": str, "turn_idx": int}],
    "tool_errors": [{"tool_name": str, "error_text": str, "timestamp": str}],
    "user_messages": [{"text": str, "is_correction": bool, "timestamp": str}],
    "subagent_delegations": [{"subagent_type": str, "agent_id": str, "prompt_snippet": str}],
    "token_usage": {"input": int, "output": int, "cache_read": int, "cache_created": int},
    "has_compact": bool,
    "api_errors": int,
    "billing_errors": int,
}
```

### 4.3 Recommended Pipeline

```
Phase 1: Enumerate (1 pass, ~30s)
  - Walk /projects/ directories
  - Classify each JSONL as: main-session / subagent / orphan-agent
  - Build parent-child linkage from subagent file paths
  - Output: sessions.json index

Phase 2: Extract (1 pass per session, ~2min)
  - Parse each JSONL line by line
  - Extract: tool_calls, tool_errors, user_messages, token_usage, metadata
  - Apply frustration filters (exclude system-generated messages)
  - Output: sessions_extracted.jsonl (one JSON per session)

Phase 3: Score (in-memory, ~5s)
  - Apply signal detectors to each session
  - Compute: misuse_score, error_rate, frustration_count, waste_count
  - Aggregate: cross-session totals, top-offending sessions
  - Output: kaizen_report.json

Phase 4: Report (template-driven)
  - Top N sessions by each metric
  - Example commands/patterns for each issue type
  - Temporal trends (by date)
```

### 4.4 Key Implementation Notes

1. **Line-by-line parsing** — Do not `json.loads(file.read())`. JSONL files are line-delimited.

2. **Subagent discovery:**

   ```python
   subagent_path = Path(f"{proj_dir}/{session_id}/subagents/agent-{agent_id}.jsonl")
   ```

3. **Tool result linkage:** Match `user.message.content[].tool_use_id` to `assistant.message.content[].id` in previous turn.

4. **Agent progress records** (type=progress, data.type=agent_progress) contain inner assistant turns with tool calls. For subagent tool usage analysis, both:
   - The subagent's own JSONL file (authoritative)
   - The parent's `progress.data.agent_progress` records (streaming duplicate)

   Use only the subagent JSONL to avoid double-counting.

5. **Filtering system messages:** A `user` record is system-generated if:
   - `content` starts with `"Context: This summary"`
   - `content` starts with `"This session is being continued"`
   - `content` contains `<local-command-caveat>`
   - `content` contains `<system-reminder>`
   - `content` contains `Stop hook feedback:`
   - `content` contains `<skill-format>true</skill-format>`

---

## 5. Open Questions

### 5.1 Schema Questions

1. **`agent-*.jsonl` in main dir vs `subagents/` dir** — The 106 top-level `agent-*.jsonl` files appear to be orphaned agent transcripts (agents that completed but whose parent session wasn't written to a subdir). Their `sessionId` field links them to a parent session, but the parent may be in a different project directory. Needs verification.

2. **`tool-results/*.txt` files** — These contain the full output of async Task agents. The naming convention is `{tool_use_id}.txt`. When an async Task completes, is the output also appended to the parent's JSONL as a `user` tool_result? Investigation needed.

3. **`queue-operation` content** — The `content` field is sometimes a JSON string, sometimes a plain string. The enqueue/dequeue pair seems to be the conversation summary generation queue. No analysis value confirmed beyond filtering.

4. **`file-history-snapshot.snapshot.trackedFileBackups`** — Contains file backup locations. Could enable "files changed during session" analysis without parsing Edit tool calls. Potential cross-check signal.

### 5.2 Analysis Questions

5. **Stop hook content** — The `stop_hook_summary` records have `hasOutput: true` but the hook output text itself appears in the next `user` turn as system message text (the "Stop hook feedback:" messages). Need to confirm if hook output is always in the subsequent user turn, or sometimes in a `system-reminder` tag within a later turn.

6. **Temporal correlation** — Can we correlate tool misuse rates with Claude Code version numbers? The version field is available on every record. Versions range from 2.0.50 to 2.1.45 in this corpus.

7. **Which bash commands are genuinely misuse vs. legitimate?** — A bash command `grep -n "pattern" file.py` is clearly misuse (should use Grep tool). But `git log | grep "feat:" | head -10` is a legitimate pipeline. Disambiguation requires NLP or a more detailed regex ruleset.

8. **Cross-project analysis** — There are 40 project directories total. Most have far fewer sessions, but some (like `-home-ubuntulinuxqa2-repos-centrios-stm32`) represent different domains where different tool patterns would be expected. Cross-project comparison needs project-type awareness.

### 5.3 Infrastructure Questions

9. **Real-time vs. batch** — The analysis is designed for batch processing. If real-time kaizen feedback (e.g., hook-based detection) is desired, the architecture changes significantly (single-session streaming analysis).

10. **Embedding similarity** — Some signals (instruction drift, semantic duplicate detection) require text embeddings. What embedding model/service is available in this environment?

---

## 6. Quick Wins (Highest ROI)

Based on corpus analysis, these three signals are **trivial to extract** and **high-value**:

### Signal 1: Bash File-Op Misuse Rate

- **Where:** `assistant.message.content[].input.command` for Bash calls
- **How:** Regex for `grep\b`, `find\b.*-name`, `cat\b.*\.md`, `head\b.*\.\w+`
- **Why:** 593 violations in this project alone

### Signal 2: Edit-Before-Read Error

- **Where:** `user.message.content[].content` for tool_results with `is_error: true`
- **Pattern:** Error text contains "File has not been read yet"
- **Why:** Direct detection of anti-pattern, no ambiguity

### Signal 3: User Interrupt Rate Per Session

- **Where:** `user.message.content` for non-toolUseResult messages containing "[Request interrupted"
- **Why:** Strong frustration proxy, zero false positives

---

_Analysis performed with direct observation of JSONL corpus. No training data assumptions used. All counts and schemas verified against sampled sessions._

---

## 7. Anthropic Official Skills Repository — Token Count Analysis

**Date:** 2026-02-17
**Source:** (Anthropic official skills repository)
**Tool:** `/plugin-creator:lint --tokens-only ../skills/`
**Skills found:** 16 SKILL.md files

### 7.1 Summary Statistics

```text
Total skills:    16
Total tokens:    33,570

Biggest skill:   docx           (4,945 tokens)
Smallest skill:  internal-comms   (326 tokens)

Mean:            2,098.12 tokens
Median:          2,030.50 tokens
Std deviation:   1,407.59 tokens
```

### 7.2 All Skills Sorted by Token Count (Descending)

| Skill Name            | Tokens |
| --------------------- | ------ |
| docx                  | 4,945  |
| algorithmic-art       | 4,150  |
| skill-creator         | 3,705  |
| doc-coauthoring       | 3,289  |
| xlsx                  | 2,814  |
| pptx                  | 2,403  |
| canvas-design         | 2,343  |
| pdf                   | 2,079  |
| slack-gif-creator     | 1,982  |
| mcp-builder           | 1,922  |
| webapp-testing        | 881    |
| frontend-design       | 858    |
| web-artifacts-builder | 702    |
| theme-factory         | 654    |
| brand-guidelines      | 517    |
| internal-comms        | 326    |

### 7.3 Observations

- The distribution is right-skewed: mean (2,098) > median (2,030), with the top 4 skills (docx, algorithmic-art, skill-creator, doc-coauthoring) pulling the mean up.
- The document-format skills (docx, xlsx, pptx, pdf) cluster in the 2,000–5,000 token range, consistent with needing detailed format-specific instructions.
- The bottom 4 skills (webapp-testing, frontend-design, web-artifacts-builder, theme-factory, brand-guidelines, internal-comms) are all under 900 tokens — significantly thinner than the median.
- `internal-comms` at 326 tokens is an outlier on the low end; it is ~15x smaller than `docx`.
- The high standard deviation (1,408 tokens) relative to the mean (2,098) indicates no consistent sizing norm across the official skill set.

_Data collected 2026-02-17. Token counts produced by `plugin_validator.py --tokens-only` run against each SKILL.md file directly._
