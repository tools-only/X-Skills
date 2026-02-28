# hooks-guide Skill Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a cross-platform `hooks-guide` skill to `plugin-creator` that covers hooks for Claude Code, GitHub Copilot, and other AI assistants — with Node.js CJS and Python authoring guides, inline-agent hook conventions, a fetch-and-transform pipeline, and updates to the existing `hook-creator` agent.

**Architecture:** Single skill (`plugins/plugin-creator/skills/hooks-guide/`) with a navigation SKILL.md and per-topic reference files. Three already-fetched raw docs in `.claude/plan/plugin-creator-hooks/` are transformed into AI-facing reference content via `rwr:doc-to-skill`. A bash script orchestrates fetching and transforming platform docs. The existing `hook-creator` agent and `claude-hooks-reference-2026` umbrella are updated to include the new skill.

**Tech Stack:** Bash (fetch script), `CLAUDECODE= claude -p` CLI (rwr:doc-to-skill transform), Markdown (reference files), YAML (SKILL.md frontmatter), `uv run` (plugin validator)

---

## Pre-flight: verify raw docs exist

Before starting, confirm the three already-fetched source docs are present:

```bash
ls .claude/plan/plugin-creator-hooks/
# Expected: hooks-doc.md  sub-agents-doc.md  github-copilot-hooks-doc.md
```

If any are missing, re-fetch with `technical-researcher` agent before proceeding.

---

### Task 1: Scaffold the hooks-guide skill directory

**Files:**

- Create: `plugins/plugin-creator/skills/hooks-guide/SKILL.md`
- Create: `plugins/plugin-creator/skills/hooks-guide/scripts/.gitkeep`
- Create: `plugins/plugin-creator/skills/hooks-guide/references/.gitkeep`

**Step 1: Run init_skill.py to scaffold**

```bash
uv run plugins/plugin-creator/skills/skill-creator/scripts/init_skill.py \
  hooks-guide \
  --path plugins/plugin-creator/skills
```

Expected: `plugins/plugin-creator/skills/hooks-guide/` created with SKILL.md template and example dirs.

**Step 2: Remove generated example files (keep directory structure)**

```bash
rm -f plugins/plugin-creator/skills/hooks-guide/scripts/example_script.py
rm -f plugins/plugin-creator/skills/hooks-guide/references/example_reference.md
rm -f plugins/plugin-creator/skills/hooks-guide/assets/example_asset.txt
rmdir --ignore-fail-on-non-empty plugins/plugin-creator/skills/hooks-guide/assets
```

**Step 3: Verify structure**

```bash
find plugins/plugin-creator/skills/hooks-guide -type f | sort
```

Expected output:

```
plugins/plugin-creator/skills/hooks-guide/SKILL.md
plugins/plugin-creator/skills/hooks-guide/references/
plugins/plugin-creator/skills/hooks-guide/scripts/
```

**Step 4: Commit**

```bash
git add plugins/plugin-creator/skills/hooks-guide/
git commit -m "feat(hooks-guide): scaffold hooks-guide skill directory"
```

---

### Task 2: Transform Claude Code hooks doc → claude-code.md

**Files:**

- Source: `.claude/plan/plugin-creator-hooks/hooks-doc.md` (1744-line raw doc, already fetched)
- Create: `plugins/plugin-creator/skills/hooks-guide/references/claude-code.md`

**Step 1: Delegate to rwr:doc-to-skill via general-purpose agent**

Spawn a `general-purpose` subagent with this prompt:

```
You are running the rwr:doc-to-skill conversion workflow.

Source file: .claude/plan/plugin-creator-hooks/hooks-doc.md
This is the full Claude Code hooks reference (1744 lines), written for human users.

Task: Convert it into an AI-facing reference file at:
  plugins/plugin-creator/skills/hooks-guide/references/claude-code.md

Rules for AI-facing conversion:
- Extract all factual content: events, matchers, JSON schemas, exit codes, field names, examples
- Remove all human UX prose: "you can", "simply", "just", "note that", navigation hints
- Use imperative headings: "Configure PreToolUse" not "How to configure PreToolUse"
- Preserve every code example verbatim
- Group by: Event reference → JSON input schemas → JSON output schemas → Exit codes → Configuration
- Add a Table of Contents at the top (Claude Code peeks at files)
- Target: comprehensive but under 8000 tokens

Write the result to: plugins/plugin-creator/skills/hooks-guide/references/claude-code.md
```

**Step 2: Verify output exists and has content**

```bash
wc -l plugins/plugin-creator/skills/hooks-guide/references/claude-code.md
# Expected: > 100 lines
head -20 plugins/plugin-creator/skills/hooks-guide/references/claude-code.md
# Expected: ToC and frontmatter-style header
```

**Step 3: Commit**

```bash
git add plugins/plugin-creator/skills/hooks-guide/references/claude-code.md
git commit -m "feat(hooks-guide): add claude-code.md reference from hooks.md transform"
```

---

### Task 3: Transform sub-agents doc → inline-agent-hooks.md

**Files:**

- Source: `.claude/plan/plugin-creator-hooks/sub-agents-doc.md` (814-line raw doc, already fetched)
- Create: `plugins/plugin-creator/skills/hooks-guide/references/inline-agent-hooks.md`

**Step 1: Delegate to general-purpose agent**

Spawn a `general-purpose` subagent with this prompt:

```
You are running the rwr:doc-to-skill conversion workflow.

Source file: .claude/plan/plugin-creator-hooks/sub-agents-doc.md
This is the full Claude Code sub-agents reference (814 lines), written for human users.

Task: Extract ONLY the hooks/MCP/skills/memory inline-frontmatter sections and convert
them into an AI-facing reference file at:
  plugins/plugin-creator/skills/hooks-guide/references/inline-agent-hooks.md

Sections to extract and convert (ignore the rest of the doc):
- "Define hooks for subagents" section — hooks in agent frontmatter YAML
- "mcpServers" frontmatter field — inline MCP server definitions
- "skills" frontmatter field — preloading skills into agents
- "memory" frontmatter field — persistent memory scopes
- "SubagentStart" / "SubagentStop" event matchers in settings.json
- "Stop hooks in frontmatter auto-convert to SubagentStop" rule
- "background" and "isolation: worktree" fields

Rules for AI-facing conversion:
- Remove all human UX prose
- Preserve all YAML and JSON examples verbatim
- Use imperative headings
- Add a Table of Contents at top
- Include the complete supported frontmatter fields table (name, description, tools, etc.)

Write the result to: plugins/plugin-creator/skills/hooks-guide/references/inline-agent-hooks.md
```

**Step 2: Verify**

```bash
wc -l plugins/plugin-creator/skills/hooks-guide/references/inline-agent-hooks.md
# Expected: > 80 lines
grep -l "SubagentStop\|mcpServers\|memory:" \
  plugins/plugin-creator/skills/hooks-guide/references/inline-agent-hooks.md
# Expected: file listed (confirms key sections present)
```

**Step 3: Commit**

```bash
git add plugins/plugin-creator/skills/hooks-guide/references/inline-agent-hooks.md
git commit -m "feat(hooks-guide): add inline-agent-hooks.md from sub-agents doc transform"
```

---

### Task 4: Transform GitHub Copilot hooks doc → github-copilot.md

**Files:**

- Source: `.claude/plan/plugin-creator-hooks/github-copilot-hooks-doc.md` (already fetched)
- Create: `plugins/plugin-creator/skills/hooks-guide/references/github-copilot.md`

**Step 1: Delegate to general-purpose agent**

Spawn a `general-purpose` subagent with this prompt:

```
You are running the rwr:doc-to-skill conversion workflow.

Source file: .claude/plan/plugin-creator-hooks/github-copilot-hooks-doc.md
This is the GitHub Copilot coding agent hooks reference, written for human users.

Task: Convert it into an AI-facing reference file at:
  plugins/plugin-creator/skills/hooks-guide/references/github-copilot.md

Rules for AI-facing conversion:
- Extract all factual content: hook events (camelCase), JSON schema, bash/powershell fields,
  cwd, env, timeoutSec, version: 1 requirement, .github/hooks/ location
- Remove all human UX prose, troubleshooting narrative, step-by-step UI instructions
- Preserve all JSON and shell examples verbatim
- Group by: File location → Schema → Events → Exit codes → Debugging
- Add a Table of Contents at top
- Add a "Differences from Claude Code hooks" section at the end noting:
  - camelCase event names vs PascalCase
  - bash/powershell dual keys vs single command
  - .github/hooks/ location vs hooks/hooks.json
  - version: 1 required field
  - No matcher concept (hooks apply to all calls of that type)

Write the result to: plugins/plugin-creator/skills/hooks-guide/references/github-copilot.md
```

**Step 2: Verify**

```bash
grep -c "sessionStart\|preToolUse\|postToolUse" \
  plugins/plugin-creator/skills/hooks-guide/references/github-copilot.md
# Expected: > 0
```

**Step 3: Commit**

```bash
git add plugins/plugin-creator/skills/hooks-guide/references/github-copilot.md
git commit -m "feat(hooks-guide): add github-copilot.md reference from copilot hooks doc transform"
```

---

### Task 5: Author hooks-cjs.md

**Files:**

- Reference: `plugins/plugin-creator/skills/hooks-core-reference/SKILL.md` (read for patterns)
- Reference: `plugins/plugin-creator/skills/hooks-patterns/SKILL.md` (read for templates)
- Create: `plugins/plugin-creator/skills/hooks-guide/references/hooks-cjs.md`

**Step 1: Delegate to contextual-ai-documentation-optimizer agent**

Spawn with this prompt:

```
Write an AI-facing Node.js CommonJS hook authoring guide for Claude Code hooks.

Output file: plugins/plugin-creator/skills/hooks-guide/references/hooks-cjs.md

Source material to read first (do not summarise — extract patterns):
- plugins/plugin-creator/skills/hooks-core-reference/SKILL.md
- plugins/plugin-creator/skills/hooks-patterns/SKILL.md
- plugins/plugin-creator/agents/hook-creator.md (the mandatory constraints section)

The guide must cover:
1. File naming: .cjs extension only (never .js), why (ESM risk)
2. Shebang + 'use strict' header
3. stdin parsing pattern (readFileSync /dev/stdin + try/catch exit 0)
4. execFileSync over execSync — exact pattern with stdio suppression
5. Timeout discipline — local binary: 3000ms, filesystem: 5000ms
6. Exit codes: 0 success, 1 script error, 2 blocking
7. stdout: JSON only via console.log(JSON.stringify(output))
8. stderr: error messages only (shown to Claude on exit 2)
9. Complete templates for each use case:
   - Blocking (PreToolUse, exit 2)
   - Permission decision (permissionDecision JSON)
   - Context injection (SessionStart, additionalContext field)
   - Task verification (Stop/SubagentStop, exit 2 to force continue)
   - Binary availability check (execFileSync pattern)
10. Anti-patterns section with wrong/correct pairs

Format: imperative headings, all code in fenced blocks with language specifiers,
Table of Contents at top. No prose padding.
```

**Step 2: Verify**

```bash
grep -c "execFileSync\|'use strict'\|\.cjs" \
  plugins/plugin-creator/skills/hooks-guide/references/hooks-cjs.md
# Expected: > 3
```

**Step 3: Commit**

```bash
git add plugins/plugin-creator/skills/hooks-guide/references/hooks-cjs.md
git commit -m "feat(hooks-guide): add hooks-cjs.md Node.js CJS authoring guide"
```

---

### Task 6: Author hooks-python.md

**Files:**

- Reference: `.claude/plan/plugin-creator-hooks/hooks-doc.md` (for stdin schema and exit codes)
- Create: `plugins/plugin-creator/skills/hooks-guide/references/hooks-python.md`

**Step 1: Delegate to contextual-ai-documentation-optimizer agent**

Spawn with this prompt:

```
Write an AI-facing Python hook authoring guide for Claude Code hooks.

Output file: plugins/plugin-creator/skills/hooks-guide/references/hooks-python.md

Source to read for stdin schema and exit code spec:
- .claude/plan/plugin-creator-hooks/hooks-doc.md (sections: hook input, exit codes, JSON output)

The guide must cover:
1. Reading stdin: sys.stdin.read() + json.loads() + try/except sys.exit(0)
2. Exit codes: 0 success, 1 script error, 2 blocking
3. stdout: json.dumps() only — never print() for hook responses
4. stderr: sys.stderr.write() for error messages (shown to Claude on exit 2)
5. Subprocess calls: subprocess.run(['binary', 'arg'], capture_output=True, timeout=3) pattern
6. Shebang: #!/usr/bin/env python3
7. Complete templates for each use case:
   - Blocking (PreToolUse, exit 2)
   - Permission decision (permissionDecision dict → json.dumps)
   - Context injection (SessionStart, additionalContext)
   - Task verification (Stop, exit 2 to force continue)
8. When to use Python vs CJS:
   - Python: when integrating with Python tooling (ruff, mypy, pytest), file parsing, regex
   - CJS: default for plugin hooks (avoids Python env dependency)
9. Anti-patterns: using print() for hook output, not catching JSON parse errors

Format: imperative headings, all code fenced with python specifier,
Table of Contents at top. No prose padding.
```

**Step 2: Verify**

```bash
grep -c "sys.stdin\|json.dumps\|sys.exit" \
  plugins/plugin-creator/skills/hooks-guide/references/hooks-python.md
# Expected: > 3
```

**Step 3: Commit**

```bash
git add plugins/plugin-creator/skills/hooks-guide/references/hooks-python.md
git commit -m "feat(hooks-guide): add hooks-python.md Python authoring guide"
```

---

### Task 7: Author common-schema.md and best-practices.md

**Files:**

- Create: `plugins/plugin-creator/skills/hooks-guide/references/common-schema.md`
- Create: `plugins/plugin-creator/skills/hooks-guide/references/best-practices.md`

**Step 1: Delegate both to contextual-ai-documentation-optimizer agent (two separate spawns)**

**For common-schema.md:**

```
Write an AI-facing cross-platform hooks common schema reference.

Output file: plugins/plugin-creator/skills/hooks-guide/references/common-schema.md

Source material:
- .claude/plan/plugin-creator-hooks/hooks-doc.md (Claude Code full schema)
- plugins/plugin-creator/skills/hooks-guide/references/github-copilot.md (already written)

Cover:
1. Concept map: what "hook" means across platforms (event → trigger → command/script)
2. Common event categories present across platforms:
   - Session lifecycle (start/end)
   - Tool/action pre/post
   - User input submission
   - Agent/subagent lifecycle
3. Common fields present across platforms: timeout, working directory, env vars
4. Exit code conventions (0/1/2 for Claude Code; 0/non-zero for Copilot)
5. stdin JSON pattern (Claude Code) vs env vars (some platforms)
6. JSON output control: suppressOutput, hookSpecificOutput (Claude Code only)
7. Comparison table: Claude Code vs GitHub Copilot — event names, schema, location, language

Format: imperative headings, Table of Contents at top, comparison table in markdown.
```

**For best-practices.md:**

```
Write an AI-facing cross-platform hooks best practices guide.

Output file: plugins/plugin-creator/skills/hooks-guide/references/best-practices.md

Source material to read:
- plugins/plugin-creator/agents/hook-creator.md (mandatory constraints + anti-patterns sections)
- .claude/plan/plugin-creator-hooks/hooks-doc.md (async hooks, prompt hooks sections)
- plugins/plugin-creator/skills/hooks-guide/references/github-copilot.md

Cover:
1. Timeout discipline — values by operation type, why network ops are inappropriate
2. Idempotency — hooks may fire multiple times, must be safe to repeat
3. Stdin hygiene — always wrap parse in try/except/catch, exit 0 on bad input
4. stderr discipline — only error messages, never debug noise
5. stdout discipline — JSON only, never raw text in hook responses
6. Binary availability checks before shelling out
7. Blocking vs non-blocking — when exit 2 is appropriate
8. Async hooks (Claude Code) — when to use, timeout implications
9. Prompt hooks (Claude Code) — when to use vs command hooks
10. Testing hooks locally before wiring (pipe stdin, check exit code, validate JSON)
11. Anti-patterns: shell injection, stderr leak, deleting hooks.json, .js extension in ESM projects

Format: imperative headings, Table of Contents at top. Anti-patterns as wrong/correct pairs.
```

**Step 2: Verify both files**

```bash
wc -l plugins/plugin-creator/skills/hooks-guide/references/common-schema.md \
        plugins/plugin-creator/skills/hooks-guide/references/best-practices.md
# Expected: both > 60 lines
```

**Step 3: Commit**

```bash
git add plugins/plugin-creator/skills/hooks-guide/references/common-schema.md \
        plugins/plugin-creator/skills/hooks-guide/references/best-practices.md
git commit -m "feat(hooks-guide): add common-schema.md and best-practices.md"
```

---

### Task 8: Write platform-coverage.md

**Files:**

- Create: `plugins/plugin-creator/skills/hooks-guide/references/platform-coverage.md`

**Step 1: Write the file directly (no transform needed — static registry)**

Content to write:

```markdown
# Platform Coverage Registry

Tracks which AI assistant hook systems are covered, their doc URLs, and fetch status.
Update `Last verified` when re-running fetch-and-transform-hooks-docs.sh.

## Covered platforms

| Platform | Hook concept | Doc URL | Reference file | Last verified |
|----------|-------------|---------|----------------|---------------|
| Claude Code | Yes — hooks.json + settings.json | https://code.claude.com/docs/en/hooks.md | claude-code.md | 2026-02-27 |
| Claude Code (inline agent) | Yes — agent frontmatter | https://code.claude.com/docs/en/sub-agents.md | inline-agent-hooks.md | 2026-02-27 |
| GitHub Copilot | Yes — .github/hooks/ | https://docs.github.com/en/copilot/how-tos/use-copilot-agents/coding-agent/use-hooks.md | github-copilot.md | 2026-02-27 |

## Fetch-attempted platforms (no hooks doc found at time of last run)

| Platform | URL attempted | Result | Notes |
|----------|--------------|--------|-------|
| Cursor | https://docs.cursor.com/context/rules | Pending first run | Rules system, not hooks |
| Windsurf | https://docs.windsurf.com/windsurf/memories | Pending first run | Memories system, not hooks |
| Amp | https://ampcode.com/docs | Pending first run | Unverified |

## Adding a new platform

1. Add row to the table above with URL and "Pending" status
2. Add fetch entry to scripts/fetch-and-transform-hooks-docs.sh
3. Run the script: `bash plugins/plugin-creator/skills/hooks-guide/scripts/fetch-and-transform-hooks-docs.sh`
4. Update Last verified date if reference file was created/updated
```

**Step 2: Commit**

```bash
git add plugins/plugin-creator/skills/hooks-guide/references/platform-coverage.md
git commit -m "feat(hooks-guide): add platform-coverage.md registry"
```

---

### Task 9: Write SKILL.md

**Files:**

- Modify: `plugins/plugin-creator/skills/hooks-guide/SKILL.md`

**Step 1: Delegate to contextual-ai-documentation-optimizer agent**

Spawn with this prompt:

```
Write the SKILL.md for the hooks-guide skill in plugin-creator.

Output file: plugins/plugin-creator/skills/hooks-guide/SKILL.md
(Overwrite the scaffolded template completely)

This skill is a navigation hub. Its job is to route to the correct reference file.
Keep SKILL.md lean — all detail lives in references/.

Frontmatter to use exactly:
---
description: 'Cross-platform hooks reference for AI coding assistants — Claude Code, GitHub Copilot, Cursor, Windsurf, Amp. Covers hook authoring in Node.js CJS and Python, per-platform event schemas, inline-agent hooks and MCP in agent frontmatter, common JSON I/O, exit codes, best practices, and a fetch script to refresh docs from official sources. Use when writing, reviewing, or debugging hooks for any AI assistant.'
allowed-tools: Read, Grep, Glob, Bash, Write, Edit
---

Body structure (write this content):

1. Quick-start decision flowchart (Mermaid flowchart TD):
   - Which platform? → Claude Code / GitHub Copilot / Other
   - Which language? → Node.js CJS / Python
   - Inline agent hooks? → yes/no
   Each leaf node names the reference file to read.

2. Reference file index (brief list — one line each):
   - common-schema.md — shared concepts, cross-platform comparison
   - claude-code.md — Claude Code hooks full reference
   - inline-agent-hooks.md — hooks/mcpServers/skills/memory in agent frontmatter
   - github-copilot.md — GitHub Copilot coding agent hooks
   - hooks-cjs.md — Node.js CJS authoring guide and templates
   - hooks-python.md — Python authoring guide and templates
   - best-practices.md — cross-platform conventions and anti-patterns
   - platform-coverage.md — known platforms, fetch URLs, coverage status

3. Refresh instructions (2-3 lines):
   Run: bash plugins/plugin-creator/skills/hooks-guide/scripts/fetch-and-transform-hooks-docs.sh
   This re-fetches all platform docs and runs rwr:doc-to-skill on each.

4. Sources section with all three fetched doc URLs and access date 2026-02-27.

Rules: imperative headings, no prose padding, no table of contents needed (body is short),
all code fenced with language specifiers, Mermaid diagram for the flowchart.
```

**Step 2: Verify frontmatter is valid**

```bash
uv run plugins/plugin-creator/scripts/plugin_validator.py \
  plugins/plugin-creator/skills/hooks-guide
```

Fix any reported issues before committing.

**Step 3: Commit**

```bash
git add plugins/plugin-creator/skills/hooks-guide/SKILL.md
git commit -m "feat(hooks-guide): write navigation SKILL.md"
```

---

### Task 10: Write fetch-and-transform-hooks-docs.sh

**Files:**

- Create: `plugins/plugin-creator/skills/hooks-guide/scripts/fetch-and-transform-hooks-docs.sh`

**Step 1: Delegate to python-cli-architect agent**

Note: this is a bash script, but python-cli-architect handles shell scripts too. Spawn with:

```
Write a bash script that fetches AI assistant hook documentation and transforms it
into AI-facing reference files using the rwr:doc-to-skill pipeline.

Output file:
  plugins/plugin-creator/skills/hooks-guide/scripts/fetch-and-transform-hooks-docs.sh

The script must:

1. Define an array of platforms with: name, URL, output reference file path

   Platforms:
   - claude-code | https://code.claude.com/docs/en/hooks.md | references/claude-code.md
   - inline-agent | https://code.claude.com/docs/en/sub-agents.md | references/inline-agent-hooks.md
   - github-copilot | https://docs.github.com/en/copilot/how-tos/use-copilot-agents/coding-agent/use-hooks.md | references/github-copilot.md
   - cursor | https://docs.cursor.com/context/rules | references/cursor.md
   - windsurf | https://docs.windsurf.com/windsurf/memories | references/windsurf.md
   - amp | https://ampcode.com/docs | references/amp.md

2. For each platform:
   a. curl the URL to a temp file in /tmp/hooks-fetch-<platform>.md
   b. If curl fails (non-zero exit or empty file): log "SKIP <platform>: fetch failed" and continue
   c. If file is < 500 bytes: log "SKIP <platform>: response too small (likely no hooks content)" and continue
   d. Run: `CLAUDECODE= claude -p "You are running rwr:doc-to-skill. Convert the human-facing documentation
      in the file at /tmp/hooks-fetch-<platform>.md into an AI-facing reference file at
      <output_path>. Rules: remove UX prose, preserve all code examples verbatim, add ToC,
      use imperative headings, group by concept." 2>&1`
   e. If claude exits non-zero: log "FAIL <platform>: transform failed" and continue
   f. Log "OK <platform>: wrote <output_path>" on success

3. After all platforms, update platform-coverage.md:
   - Set Last verified date for successfully updated platforms to today's date (date +%Y-%m-%d)
   - This is a sed-in-place update on the table row

4. Exit 0 if all attempted platforms succeeded, exit 1 if any failed

Script conventions:
- Shebang: #!/usr/bin/env bash
- set -euo pipefail NOT used (graceful partial — don't abort on one platform failure)
- Explicit error handling per platform instead
- SCRIPT_DIR via $(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
- SKILL_DIR=$(dirname "$SCRIPT_DIR") for reference file paths
- Comments explaining each section

After writing, make it executable:
  chmod +x plugins/plugin-creator/skills/hooks-guide/scripts/fetch-and-transform-hooks-docs.sh
```

**Step 2: Verify the script is executable and parses cleanly**

```bash
bash -n plugins/plugin-creator/skills/hooks-guide/scripts/fetch-and-transform-hooks-docs.sh
echo "Syntax OK: $?"
ls -la plugins/plugin-creator/skills/hooks-guide/scripts/fetch-and-transform-hooks-docs.sh
# Expected: -rwxr-xr-x permissions
```

**Step 3: Commit**

```bash
git add plugins/plugin-creator/skills/hooks-guide/scripts/fetch-and-transform-hooks-docs.sh
git commit -m "feat(hooks-guide): add fetch-and-transform-hooks-docs.sh pipeline script"
```

---

### Task 11: Update hook-creator agent

**Files:**

- Modify: `plugins/plugin-creator/agents/hook-creator.md`

**Step 1: Delegate to contextual-ai-documentation-optimizer agent**

Spawn with this prompt:

```
Update the agent file at plugins/plugin-creator/agents/hook-creator.md.

Read the file first, then make these four targeted changes:

CHANGE 1 — frontmatter skills field:
Add plugin-creator:hooks-guide to the skills list.
Current skills list (read from file to get exact current value).
New skills list: existing skills + plugin-creator:hooks-guide

CHANGE 2 — Scope Determination flowchart:
Add a new branch to the existing Mermaid flowchart after the existing four scope options.
New branch: Q1 -->"Scoped to one agent only"--> InlineAgent["Inline agent frontmatter hooks\n→ hooks: field in agent .md file\n→ Scoped to that agent's lifecycle only\n→ See inline-agent-hooks.md"]

CHANGE 3 — Event Selection flowchart:
Add five new events to the existing flowchart as branches from Q1:
- WorktreeCreate — fires when a git worktree is created
- WorktreeRemove — fires when a git worktree is removed
- TeammateIdle — fires when a team agent goes idle
- TaskCompleted — fires when a task is marked complete
- ConfigChange — fires when settings are changed

CHANGE 4 — Sources section:
Update the access dates to 2026-02-27 and add:
- Inline agent hooks: https://code.claude.com/docs/en/sub-agents.md (accessed 2026-02-27)
- Cross-platform guide: plugin-creator:hooks-guide

Do NOT rewrite any other section. Make surgical edits only.
Preserve all existing content outside these four change areas.
```

**Step 2: Verify changes**

```bash
grep "hooks-guide" plugins/plugin-creator/agents/hook-creator.md
# Expected: line in skills: field

grep "InlineAgent\|Scoped to one agent\|inline-agent-hooks" \
  plugins/plugin-creator/agents/hook-creator.md
# Expected: at least one match

grep "WorktreeCreate\|TeammateIdle\|TaskCompleted" \
  plugins/plugin-creator/agents/hook-creator.md
# Expected: at least one match
```

**Step 3: Validate**

```bash
uv run plugins/plugin-creator/scripts/plugin_validator.py \
  plugins/plugin-creator/agents/hook-creator.md
```

Fix any issues before committing.

**Step 4: Commit**

```bash
git add plugins/plugin-creator/agents/hook-creator.md
git commit -m "feat(hook-creator): add hooks-guide skill, inline-agent scope, new events"
```

---

### Task 12: Update claude-hooks-reference-2026 umbrella

**Files:**

- Modify: `plugins/plugin-creator/skills/claude-hooks-reference-2026/SKILL.md`

**Step 1: Read the current umbrella SKILL.md**

```bash
cat plugins/plugin-creator/skills/claude-hooks-reference-2026/SKILL.md
```

**Step 2: Add hooks-guide to the skills list**

The umbrella loads a list of skills (e.g., `plugin-creator:hooks-core-reference`, etc.).
Add `plugin-creator:hooks-guide` to that list using the Edit tool — surgical edit only.

**Step 3: Validate**

```bash
uv run plugins/plugin-creator/scripts/plugin_validator.py \
  plugins/plugin-creator/skills/claude-hooks-reference-2026
```

**Step 4: Commit**

```bash
git add plugins/plugin-creator/skills/claude-hooks-reference-2026/SKILL.md
git commit -m "feat(claude-hooks-reference-2026): add hooks-guide to umbrella skill list"
```

---

### Task 13: Full plugin validation

**Files:** None created — validation only.

**Step 1: Run full plugin validator**

```bash
uv run plugins/plugin-creator/scripts/plugin_validator.py plugins/plugin-creator
```

Expected: all new and modified files pass. No SK006/SK007 token budget violations.

**Step 2: Fix any issues**

If validator reports token budget violations on any reference file, delegate a split to
`plugin-creator:refactor-plugin` agent.

If validator reports frontmatter issues, fix using Edit tool on the specific file.

Re-run validator after each fix until clean.

**Step 3: Commit any fixes**

```bash
git add -A
git commit -m "fix(hooks-guide): address plugin validator issues"
```

---

### Task 14: Lint all modified files

**Files:** None created — lint only.

**Step 1: Run prek on all new/modified files**

```bash
uv run prek run --files \
  plugins/plugin-creator/skills/hooks-guide/SKILL.md \
  plugins/plugin-creator/skills/hooks-guide/references/claude-code.md \
  plugins/plugin-creator/skills/hooks-guide/references/inline-agent-hooks.md \
  plugins/plugin-creator/skills/hooks-guide/references/github-copilot.md \
  plugins/plugin-creator/skills/hooks-guide/references/hooks-cjs.md \
  plugins/plugin-creator/skills/hooks-guide/references/hooks-python.md \
  plugins/plugin-creator/skills/hooks-guide/references/common-schema.md \
  plugins/plugin-creator/skills/hooks-guide/references/best-practices.md \
  plugins/plugin-creator/skills/hooks-guide/references/platform-coverage.md \
  plugins/plugin-creator/agents/hook-creator.md \
  plugins/plugin-creator/skills/claude-hooks-reference-2026/SKILL.md
```

Fix any reported issues. Re-run until clean.

**Step 2: Final commit if fixes were needed**

```bash
git add -A
git commit -m "fix(hooks-guide): lint fixes"
```

---

## Summary of commits expected

```
feat(hooks-guide): scaffold hooks-guide skill directory
feat(hooks-guide): add claude-code.md reference from hooks.md transform
feat(hooks-guide): add inline-agent-hooks.md from sub-agents doc transform
feat(hooks-guide): add github-copilot.md reference from copilot hooks doc transform
feat(hooks-guide): add hooks-cjs.md Node.js CJS authoring guide
feat(hooks-guide): add hooks-python.md Python authoring guide
feat(hooks-guide): add common-schema.md and best-practices.md
feat(hooks-guide): add platform-coverage.md registry
feat(hooks-guide): write navigation SKILL.md
feat(hooks-guide): add fetch-and-transform-hooks-docs.sh pipeline script
feat(hook-creator): add hooks-guide skill, inline-agent scope, new events
feat(claude-hooks-reference-2026): add hooks-guide to umbrella skill list
fix(hooks-guide): address plugin validator issues  [if needed]
fix(hooks-guide): lint fixes  [if needed]
```
