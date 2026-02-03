# fn-1: Memory system for persistent learning

## Overview

Add `.flow/memory/` for project-scoped persistent learnings that survive context compaction. The core insight: **reviewers surface things the implementing model misses**. By capturing NEEDS_WORK feedback, the system learns patterns, framework quirks, and project conventions that models tend to overlook.

## Problem

1. Context compaction loses important decisions and learnings
2. Each session/iteration starts fresh, repeating past mistakes
3. Reviewer feedback is valuable but ephemeral - fixed once, forgotten next time
4. Models consistently miss the same types of things (framework idioms, project conventions, edge cases)

## Design Principles

1. **Simple**: Markdown files, no database, no worker service
2. **Project-scoped**: Lives in `.flow/memory/`, travels with the repo
3. **Reviewer-driven**: Primary signal is NEEDS_WORK feedback from cross-model reviews
4. **Selective capture**: Only record learnings that models tend to miss:
   - Framework/library-specific patterns
   - Project conventions not in CLAUDE.md
   - Common pitfalls for this codebase
   - Architectural decisions with non-obvious rationale
5. **Subagent retrieval**: Memory-scout finds relevant entries without polluting main context

## Directory Structure

```
.flow/
  config.json         # Project settings (memory enabled, etc.)
  memory/
    pitfalls.md       # Things that went wrong, extracted from NEEDS_WORK feedback
    conventions.md    # Project patterns discovered during work
    decisions.md      # Architectural choices with rationale
```

## Configuration

Settings live in `.flow/config.json`:

```json
{
  "memory": {
    "enabled": false
  }
}
```

| Setting | Default | Description |
|---------|---------|-------------|
| `memory.enabled` | `false` | Enable memory system (opt-in for now) |

One switch. When enabled: read + write both active.

**Check in code:**
```python
# In flowctl.py or hooks
config = load_flow_config()
if not config.get("memory", {}).get("enabled", False):
    return  # Memory disabled, skip

# In skills
if memory_enabled:
    Task flow-next:memory-scout(<request>)
```

**Enable for a project:**
```bash
flowctl config set memory.enabled true
# or manually edit .flow/config.json
```

## Backwards Compatibility

Memory is **opt-in** and disabled by default. This ensures:

- Existing plan/work/review workflows unchanged when disabled
- Ralph prompts and receipt flow untouched
- No subagents spawn, no hooks fire when disabled
- Zero overhead for users who don't enable memory

## Integration with existing patterns

### flowctl init
Update `cmd_init()` in flowctl.py to create `config.json` with defaults:

```python
# In cmd_init(), after creating directories
config_path = flow_dir / "config.json"
if not config_path.exists():
    default_config = {
        "memory": {
            "enabled": False
        }
    }
    config_path.write_text(json.dumps(default_config, indent=2))
```

### First run / missing config
Any command that needs config should handle missing file gracefully:

```python
def load_flow_config() -> dict:
    config_path = get_flow_dir() / "config.json"
    if not config_path.exists():
        return {}  # Use defaults
    return json.loads(config_path.read_text())

def get_config(key: str, default=None):
    """Get nested config value like 'memory.enabled'."""
    config = load_flow_config()
    for part in key.split('.'):
        config = config.get(part, {})
        if config == {}:
            return default
    return config if config != {} else default
```

### /flow-next:setup
Setup creates config if missing (same as init), then mentions it:

```markdown
# In workflow.md Step 1 (after creating .flow/)
Check if .flow/config.json exists:
- If missing: create with defaults (memory.enabled: false)
- If exists: leave as-is

# In Step 7 summary
Memory system: disabled by default
Enable with: flowctl config set memory.enabled true
```

Both `flowctl init` and `/flow-next:setup` ensure config exists - whichever runs first creates it.

### Config gating (both read and write)

**Read operations (subagent):**
```markdown
# In /flow-next:plan steps.md Step 1 (parallel with other scouts)
# Add to the parallel Task calls if memory.enabled:
- Task flow-next:memory-scout(<request>)

# In /flow-next:work phases.md Phase 1 (re-anchor)
# Add after flowctl show/cat if memory.enabled:
- Task flow-next:memory-scout(<task-id>: <task-title>)
```

Same subagent, same call pattern, different input context.

**Write operations (hook):**
```python
# In ralph-guard.py PostToolUse handler
if tool_name == "Bash" and "chat-send" in command:
    if not get_config("memory.enabled", False):
        return  # Memory disabled, skip
    # ... parse NEEDS_WORK and append to pitfalls.md
```

One setting controls both read and write.

## What To Capture (and What NOT To)

### Capture (things models miss repeatedly)

| Category | Example | Why Models Miss It |
|----------|---------|-------------------|
| Framework idioms | "React 19: use() instead of useEffect for data" | Training data staleness |
| Project conventions | "This repo uses barrel exports in index.ts" | Not in CLAUDE.md |
| API quirks | "flowctl rp requires --repo-root flag" | Wrapper-specific |
| Build/test gotchas | "Must run npm test from repo root" | Implicit assumptions |
| Edge cases | "Empty arrays need special handling in X" | Discovered via review |

### DON'T Capture

- Generic coding best practices (models know these)
- Obvious bugs (typos, syntax errors)
- One-off mistakes unlikely to recur
- Sensitive data (API keys, credentials)

## Implementation

### Phase 1: Capture (PostToolUse hook)

Extend `ralph-guard.py` to capture learnings from NEEDS_WORK reviews:

```python
# In PostToolUse handler for chat-send
if "NEEDS_WORK" in response or "MAJOR_RETHINK" in response:
    feedback = extract_feedback(response)
    if is_learnable(feedback):  # Filter to actionable patterns
        append_to_memory("pitfalls", {
            "date": today,
            "task": task_id,
            "issue": feedback.issue,
            "fix": feedback.fix,
            "category": classify(feedback)  # framework|convention|api|edge-case
        })
```

**Filtering logic** (`is_learnable`):
- Has specific actionable fix (not vague "improve this")
- References code pattern, API, or convention (not just "wrong output")
- Not a one-off typo or obvious bug

### Phase 2: Retrieval (memory-scout subagent)

One subagent used in both plan and work phases:

```markdown
# agents/memory-scout.md

You search `.flow/memory/` for entries relevant to the current task.

Input: <context> (planning request OR task title/description)

Steps:
1. Read `.flow/memory/pitfalls.md`, `conventions.md`, `decisions.md`
2. Find entries semantically related to the context
3. Return ONLY relevant entries (not the whole memory)

Output format:
## Relevant Memory
- [pitfall] When doing X, remember to Y (from fn-1.1)
- [convention] This repo uses Z pattern for W
- [decision] We chose X over Y because Z
```

**Called from two places (matches existing scout pattern):**
```markdown
# In /flow-next:plan steps.md Step 1 (parallel with repo-scout, practice-scout, docs-scout)
- Task flow-next:memory-scout(<request>)

# In /flow-next:work phases.md Phase 1 (after flowctl show/cat)
- Task flow-next:memory-scout(<task-id>: <task-title>)
```

Same call pattern as other scouts. Subagent returns only relevant entries.

### Phase 3: Configuration scope

**Memory config is separate from Ralph config:**

| Config | Location | Scope |
|--------|----------|-------|
| Memory | `.flow/config.json` | All flow-next usage (manual + Ralph) |
| Ralph | `scripts/ralph/env.config` | Ralph loop only |

Memory works in:
- `/flow-next:plan` (manual)
- `/flow-next:work` (manual)
- Ralph mode (automated)

Document clearly: "Memory is a flow-next feature, not Ralph-specific. Enable via `.flow/config.json`, not `env.config`."

## flowctl Commands

```bash
# Auto-populated by hooks
flowctl memory init                    # Create .flow/memory/ structure

# Manual additions (for decisions, conventions discovered manually)
flowctl memory add --type pitfall "Always use flowctl rp wrappers"
flowctl memory add --type convention "Error handlers go in src/utils/errors.ts"
flowctl memory add --type decision "Chose SQLite over Postgres for simplicity"

# Retrieval
flowctl memory read                    # Dump all memory
flowctl memory read --type pitfalls    # Filter by type
flowctl memory list                    # Show entry count per file
flowctl memory search "pattern"        # Regex search across all memory
```

## Memory File Format

```markdown
# pitfalls.md

## 2026-01-10 fn-1.1 impl-review [framework]
**Issue**: placeholder export still present after adding new export
**Fix**: Remove unused exports before review - reviewers flag dead code
**Context**: src/index.ts had `export const placeholder = 0` from scaffold

## 2026-01-09 fn-2.3 impl-review [api]
**Issue**: flowctl rp builder called without --window flag
**Fix**: Always use `flowctl rp setup-review` which handles window selection
**Context**: Direct builder calls require manual window ID lookup
```

## Curation (preventing bloat)

Start simple, add complexity only if needed:

1. **v1 (MVP)**: Global memory, manual review, no limits
2. **v2 (if bloat)**: Add `[archived]` tag, exclude from default read
3. **v3 (if still bloat)**: LLM summarization to condense old entries
4. **v4 (large projects)**: Epic-scoped memory subdirs

If memory gets large, memory-scout subagent returns only relevant entries (not everything).

## Tasks

### fn-1.1: Config + memory directory structure
- [ ] Add `.flow/config.json` support to flowctl
- [ ] `flowctl config get <key>` / `flowctl config set <key> <value>`
- [ ] Default config with `memory.enabled: false`
- [ ] Add `flowctl memory init` command (creates `.flow/memory/` if enabled)
- [ ] Create memory templates: pitfalls.md, conventions.md, decisions.md
- [ ] Memory commands check `memory.enabled` before running

### fn-1.2: Capture hook
- [ ] Extend ralph-guard.py PostToolUse for chat-send
- [ ] Parse NEEDS_WORK feedback into structured format
- [ ] Implement `is_learnable()` filter
- [ ] Append to `.flow/memory/pitfalls.md`

### fn-1.3: flowctl memory commands
- [ ] `flowctl memory add --type <type> "<content>"`
- [ ] `flowctl memory read [--type <type>]`
- [ ] `flowctl memory list`
- [ ] `flowctl memory search "<pattern>"` (regex)

### fn-1.4: Memory-scout subagent
- [ ] Create `agents/memory-scout.md`
- [ ] Add to plan skill as 4th parallel subagent (if memory.enabled)
- [ ] Add to work skill re-anchor phase (if memory.enabled)
- [ ] Test relevance matching in both contexts

### fn-1.5: Integration testing
- [ ] Test in manual `/flow-next:plan`
- [ ] Test in manual `/flow-next:work`
- [ ] Test in Ralph mode end-to-end

### fn-1.6: Documentation
- [ ] Add memory section to flow-next README
- [ ] Document in ralph.md (note: separate from Ralph config)
- [ ] Add to CLAUDE.md project guide
- [ ] Clarify: memory is flow-next feature, works in manual + Ralph

## Acceptance Criteria

- [ ] NEEDS_WORK feedback auto-appends to pitfalls.md (filtered)
- [ ] `flowctl memory` commands work
- [ ] Memory-scout runs in parallel during planning
- [ ] Re-anchoring includes relevant memory
- [ ] Ralph prompts include memory context
- [ ] Memory-scout returns relevant entries only (handles large memory gracefully)
- [ ] Manual additions via `flowctl memory add` work

## Comparison: flow-next vs claude-mem

| Aspect | claude-mem | flow-next memory |
|--------|------------|------------------|
| Scope | Cross-project | Project-scoped |
| Storage | SQLite + Chroma | `.flow/memory/*.md` |
| Capture | All tool use | NEEDS_WORK feedback only |
| Retrieval | MCP tools + worker | Subagent + flowctl |
| Complexity | High (worker, DB, vectors) | Low (files, grep, subagent) |
| Signal | Everything | Reviewer insights only |

**Our advantage**: Reviewers (different model) catch things the implementing model misses. By capturing only NEEDS_WORK feedback, we get high-signal learnings without noise.

## Quick Commands

```bash
# Validate flowctl
plugins/flow-next/scripts/flowctl validate --all

# Test memory commands (after implementation)
plugins/flow-next/scripts/flowctl memory init
plugins/flow-next/scripts/flowctl memory add --type pitfall "test entry"
plugins/flow-next/scripts/flowctl memory read
```

## References

- plans/flow-next.md - original roadmap mention
- Anthropic long-context guidance on re-anchoring
- claude-mem architecture: /tmp/claude-mem/docs/public/architecture/hooks.mdx
- Ralph E2E notes: plans/ralph-e2e-notes.md (examples of learnable issues)
