---
name: extract-learnings
description: >
  This skill should be used when the user asks to persist learnings to memory.
  Triggers on "extract learnings", "save this for next time", "remember this
  pattern", "add this to memory", "update CLAUDE.md with what we learned",
  "store this insight".
allowed-tools:
  - Read
  - Write
  - Edit
  - Glob
  - Grep
  - Bash(python3:*)
  - AskUserQuestion
---

# extract-learnings

## Memory Hierarchy

Learnings are placed based on scope, persistence needs, and attention priority. Closer to Layer 0 = loaded every session = more agent attention.

| Layer | File | When Loaded | Purpose |
|-------|------|-------------|---------|
| 0 | `~/.claude/CLAUDE.md` | Every session, all projects | Universal behavioral preferences |
| 1 | `<repo>/CLAUDE.md` | Every session, this project | Project architecture, conventions, gotchas |
| 2 | `~/.claude/projects/.../memory/MEMORY.md` | Every session, agent-managed | Evolving notes, version history, working knowledge |
| 3 | `memory/*.md` topic files | On-demand, when relevant | Detailed reference too long for Layer 2 |
| Meta | Suggest new skill/command | N/A | Repeatable workflow patterns deserve automation |

### Placement Decision Tree

For each candidate learning, ask in order:

1. **Project-independent behavioral preference?** → Layer 0 (`~/.claude/CLAUDE.md`)
   Example: "Always use bun instead of npm", "Never auto-commit without asking"

2. **Project-specific technical knowledge?** → Layer 1 (`<repo>/CLAUDE.md`)
   Example: "FTS cascade: FTS5 → FTS4 → LIKE", "Always bump version in two files"

3. **Concise working note the agent should see every session?** → Layer 2 (`MEMORY.md`)
   Example: "v0.5.1: added Python 3.7 compat", "User prefers conventional commits"

4. **Detailed reference too long for Layer 2?** → Layer 3 (topic file in `memory/`)
   Example: Deep debugging guide, complex architecture notes, long reference tables

5. **Repeatable workflow pattern?** → Meta: suggest creating a skill or command
   Example: "We keep doing X manually → propose a skill for it"

## Context Gathering Tools

Reuse the recall-conversations tools when retrieval from prior sessions is needed.

### recent_chats

```bash
python3 ${CLAUDE_PLUGIN_ROOT}/skills/recall-conversations/scripts/recent_chats.py --n 10
```

Returns markdown-formatted sessions with exchange headers, timestamps, and project paths. Use `--format json` for structured filtering.

### search_conversations

```bash
python3 ${CLAUDE_PLUGIN_ROOT}/skills/recall-conversations/scripts/search_conversations.py --query "keyword"
```

Returns matching sessions ranked by relevance (BM25 when FTS5 available). Use `--format json` for structured result parsing.

For full option catalogs, load `${CLAUDE_PLUGIN_ROOT}/skills/recall-conversations/references/tool-reference.md`. Both scripts default to markdown (token-efficient for synthesis). Use `--format json` when filtering or counting sessions programmatically.

## Workflow

Four stages: **Gather → Analyze → Propose → Execute**. Never skip the proposal stage.

### Stage 1: Context Gathering

Deduce intent from available context. Do not assume where learnings come from.

- **Current conversation signal** ("save this for next time", "remember this pattern"): Source material is already in context. Do not fetch past sessions.

- **Past session signal** ("we figured out that...", "remember when we..."): Use `search_conversations.py` or `recent_chats.py` to retrieve that context. Extract substantive keywords to build the query.

- **Broad extraction request** ("extract learnings from recent sessions"): Use `recent_chats.py --n 10` to gather recent context.

- **Ambiguous intent**: Use `AskUserQuestion` to clarify what to extract and from where. Do not guess.

Exit Stage 1 when at least one source of candidate learnings is in context (current conversation content or retrieved past sessions). If ambiguous intent remains unresolved after one `AskUserQuestion` round, stop and report.

### Stage 2: Analysis & Distillation

Identify candidate learnings from gathered context. For each candidate, determine: the learning itself (condensed to 1-2 sentences), the target layer (via the decision tree above), and the target section within that file.

Limit to 3-7 candidates per invocation. If more surface, rank by impact and present the top set.

### Stage 3: Placement Proposal

This is the critical stage. Never write without explicit approval.

1. **Read current state** of all target files that have candidates aimed at them. Only read files with pending changes.

2. **Check for duplicates** — scan target files for existing content covering the same concept. If already captured (even in different words), skip it and note the duplicate.

3. **Present the proposal** — for each learning, use this format:

   ```
   ### Learning: <one-line summary>
   **Target:** <file path> → <section name>
   **Rationale:** <why this layer, why this section>

   ```diff
   + <the exact line(s) to add>
   ```
   ```

4. **Layer 0 extra gate** — if any learning targets `~/.claude/CLAUDE.md`, add an explicit warning: "This will be added to your global instructions loaded in every session across all projects. Confirm?"

5. **MEMORY.md line check** — read the target MEMORY.md and count lines. If line count exceeds 170, warn that content beyond line 200 is truncated by auto-memory and suggest moving lower-priority sections to topic files (Layer 3) before adding new content.

6. **Get approval** — use `AskUserQuestion` with options:
   - "Approve all" — apply all proposed changes
   - "Approve selectively" — let user pick which learnings to apply
   - "Reject and redirect" — user wants changes placed differently

### Stage 4: Execution

Apply approved edits:

- **Existing files** (Layers 0-2): Use `Edit` to insert at the target section. Preserve existing structure and formatting.
- **New topic files** (Layer 3): Use `Write` to create `memory/<topic>.md`. Link from MEMORY.md if appropriate.
- **Skill suggestions** (Meta): Describe the proposed skill and let the user decide.

Output a summary table:

```
| Learning | Target | Status |
|----------|--------|--------|
| "Always use bun" | ~/.claude/CLAUDE.md | Applied |
| "FTS cascade order" | repo CLAUDE.md | Applied |
| "v0.6.0 changes" | MEMORY.md | Applied |
```

## Content Quality Rules

Every candidate must pass these filters before being proposed:
- Would the agent benefit from knowing this in future sessions?
- Is it condensed to the minimum needed to be useful?
- Is it placed at the right layer (not too broad, not too narrow)?
- Does the target file not already contain this knowledge?

### Examples of passing

- Commands or workflows discovered through trial and error
- Non-obvious gotchas that caused debugging time
- Architectural decisions and their rationale
- Behavioral corrections from the user ("don't do X", "always do Y")
- Configuration quirks not obvious from reading code
- Package/module relationships not obvious from imports
- Version history milestones

### Examples of failing

- Information obviously readable from the code itself
- Generic programming best practices ("write tests", "use meaningful names")
- One-off bug fixes with no recurring pattern
- Verbose explanations — condense to a one-liner or skip
- Content already captured in the target file (semantic dedup)
- Temporary state or session-specific context
- Speculative conclusions not verified against the codebase
