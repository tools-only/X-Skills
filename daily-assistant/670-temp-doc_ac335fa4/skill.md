---
argument-hint: "/temp-doc <document>"
description: "Manage temporary documentation created during task execution - prevent git pollution from AI-generated checklists, reports, and investigation notes."
disable-model-invocation: false
---

# Temporary Documentation Lifecycle Command

## Purpose

Manage temporary documentation created during task execution - prevent git pollution from AI-generated checklists, reports, and investigation notes.

## Script Location

`.claude/utilities/find-temp-documentation.py` - Standalone Python script with Typer CLI.

**Primary User**: Claude Code orchestrator (not human users) **Secondary User**: Sub-agents (to update document status)

## Critical Rules

**PROBLEM**: AI agents create temporary docs (bug reports, checklists, investigation notes) that:

- Pollute git history
- Confuse future AI sessions (treats fixed bugs as active)
- Create maintenance burden
- Are never read by humans or AI after task completion

**SOLUTION**: Enforce strict lifecycle for ALL temporary documentation.

## Your Task

Apply temporary documentation lifecycle policy to: $ARGUMENTS

If $ARGUMENTS is empty, audit current directory for temporary docs and recommend cleanup.

## Script Usage

### For Orchestrator (Claude Code Main Instance)

**Context-Efficient Workflow**:

```bash
# Step 1: Check summary (50-100 tokens)
.claude/utilities/find-temp-documentation.py scan --summary

# Step 2: Peek at frontmatter of specific doc (100-150 tokens)
.claude/utilities/find-temp-documentation.py show --index 0 --frontmatter

# Step 3: Read full content only if needed (use Read tool)
Read("path/from/frontmatter")

# Step 4: Before completing task, cleanup completed docs
.claude/utilities/find-temp-documentation.py cleanup --no-dry-run
```

**Common Scenarios**:

```bash
# User asks: "What bugs did you find?"
# → Find bug reports quickly
.claude/utilities/find-temp-documentation.py scan --summary
# See: bugs: 2
.claude/utilities/find-temp-documentation.py show --index 0 --frontmatter
# description: "Bug tracking for config.py. 3 bugs (1 critical)"

# Multiple parallel tasks - organize by task
.claude/utilities/find-temp-documentation.py scan --group-by task

# Check what needs cleanup
.claude/utilities/find-temp-documentation.py cleanup  # dry-run by default

# Actually cleanup (moves to .wastepaper/)
.claude/utilities/find-temp-documentation.py cleanup --no-dry-run
```

### For Sub-Agents

**After Creating Temporary Doc**:

```bash
# Report to orchestrator in final message:
"Created temporary documentation: scripts/pypis_delivery_service/tests/IMPLEMENTATION_BUGS.md
Run: .claude/utilities/find-temp-documentation.py show --index 0 --frontmatter"
```

**Update Status When Finished**:

```bash
# Mark as for-review when sub-agent completes
.claude/utilities/find-temp-documentation.py state implementation-bugs-pypis-config for-review

# Or by index
.claude/utilities/find-temp-documentation.py state --index 0 for-review
```

**Check Current Status**:

```bash
# Get status without Read+Edit overhead
.claude/utilities/find-temp-documentation.py state --index 0
```

## Orchestrator Review Workflow

When sub-agents mark documents as 'for-review', the orchestrator must assess relevance and decide next action:

### 1. Discover Documents Needing Review

```bash
# Find all for-review docs
.claude/utilities/find-temp-documentation.py scan --summary
# Output shows: for-review: 2

# See which docs need review
.claude/utilities/find-temp-documentation.py scan --group-by status
```

### 2. Review Each Document

For each 'for-review' document:

```bash
# Step 1: Peek at frontmatter to understand context
.claude/utilities/find-temp-documentation.py show --index 0 --frontmatter

# Output shows:
# - task-file: Which task this relates to
# - cleanup-trigger: When it should be deleted
# - type: What kind of doc (bugs, investigation, etc.)
# - description: What Claude Code needs to know

# Step 2: Read full content ONLY if relevance unclear
Read("path/from/frontmatter")
```

### 3. Assess Relevance

**Decision Matrix**:

| Cleanup Trigger Met? | Doc Type      | Action                                            |
| -------------------- | ------------- | ------------------------------------------------- |
| ✅ Yes               | bugs          | All bugs fixed? → state completed                 |
| ✅ Yes               | investigation | Insights extracted to code? → state completed     |
| ✅ Yes               | checklist     | All items done? → state completed                 |
| ✅ Yes               | findings      | Already incorporated into task? → state completed |
| ❌ No                | any           | Still relevant? → state active                    |
| ❌ No                | any           | No longer relevant? → state completed             |

**Examples**:

```bash
# Example 1: Bug doc, all bugs fixed
# Frontmatter shows: cleanup-trigger: "When all 3 bugs fixed"
# You've verified all bugs are fixed in code
.claude/utilities/find-temp-documentation.py state --index 0 completed

# Example 2: Investigation doc, insights extracted
# You've added the design decisions to code comments
.claude/utilities/find-temp-documentation.py state investigation-auth completed

# Example 3: Checklist, task still ongoing
# User paused work, checklist still needed
.claude/utilities/find-temp-documentation.py state --index 2 active
```

### 4. Cleanup Completed Documents

After marking documents as 'completed':

```bash
# Verify what will be cleaned
.claude/utilities/find-temp-documentation.py cleanup
# Shows: Would move 2 files to .wastepaper/

# Execute cleanup
.claude/utilities/find-temp-documentation.py cleanup --no-dry-run

# Verify cleanup worked
.claude/utilities/find-temp-documentation.py scan --summary
# Should show: completed: 0
```

### 5. Before Marking Task Complete

**CRITICAL**: Always verify no orphaned temp docs:

```bash
# Check for any docs related to current task
.claude/utilities/find-temp-documentation.py scan --group-by task

# If any docs remain for this task:
# - active → Task not complete, work ongoing
# - for-review → Review and mark completed
# - completed → Run cleanup

# Only mark task complete when:
.claude/utilities/find-temp-documentation.py scan --summary
# Shows: total: 0 (or only docs from other tasks)
```

## Execution Steps

### 1. Identify Temporary Documentation

Scan for files matching these patterns:

- `*BUGS*.md`, `*BUG*.md` (not in docs/)
- `*INVESTIGATION*.md`, `*NOTES*.md`
- `*CHECKLIST*.md`, `*TODO*.md` (outside .claude/tasks/)
- `*FINDINGS*.md`, `*SUMMARY*.md`, `*REPORT*.md`
- Any `.md` file with YAML frontmatter containing `temporary: true`

### 2. Validate Frontmatter

ALL temporary docs MUST have this YAML frontmatter:

```yaml
---
temporary: true
type: checklist|report|summary|investigation|bugs|findings
task: "Brief task description or task file path"
agent: "agent-name or orchestrator"
created: YYYY-MM-DD
cleanup_trigger: "When should this be deleted/archived?"
cleanup_action: delete|archive|extract
status: active|completed|archived
---
```

**Required keys**: `temporary`, `type`, `task`, `agent`, `created`, `cleanup_trigger`, `cleanup_action`, `status`

### 3. Determine Cleanup Action

Based on `cleanup_action` in frontmatter:

#### delete (Default for most temp docs)

- Delete file immediately when status=completed
- Use for: checklists, investigation notes, scratch work
- **99% of temporary docs should use this**

#### archive (Rare - only for historical value)

- Move to `docs/history/YYYY-MM-DD-<name>.md`
- Update frontmatter: `status: archived`
- Use ONLY for: critical bug discoveries that prevent regression
- **Requires user approval - ask first**

#### extract (For docs with reusable insights)

- Extract key insights into code comments or test docstrings
- Delete original file
- Use for: investigation notes with design decisions

### 4. Execute Cleanup

```bash
# For delete action:
rm <temp-doc>.md
git add <temp-doc>.md
git commit -m "docs: remove temporary <type> (task complete)"

# For archive action (RARE - ask user first):
mkdir -p docs/history
mv <temp-doc>.md docs/history/$(date +%Y-%m-%d)-<name>.md
# Update frontmatter status to archived
git add docs/history/
git commit -m "docs: archive temporary <type> for historical reference"

# For extract action:
# 1. Add insights to code/tests
# 2. Delete temp doc
rm <temp-doc>.md
git commit -m "docs: extract insights from <temp-doc> into code comments"
```

## Template for Sub-Agents

When a sub-agent needs to create temporary documentation, use this template:

@include templates/temporary-doc-template.md

## Sub-Agent Instructions

**What Counts as Temporary Documentation?**

Any markdown file you create that is ONLY needed during active task execution and has NO long-term value to the project:

- Bug tracking lists (while fixing bugs)
- Investigation notes (exploring an issue)
- Work-in-progress checklists (tracking sub-tasks)
- Analysis reports (understanding current state)
- Decision logs (choosing between options)
- Test strategy notes (before writing actual tests)

**What is NOT Temporary?**

- User-facing documentation (README, guides, API docs)
- Architecture Decision Records (ADRs) - these are permanent
- Test files themselves (.py test files)
- Code comments and docstrings

**CRITICAL**: When creating temporary documentation (per above definition), you MUST:

1. **Use the template** above with complete YAML frontmatter
2. **Set cleanup_action: delete** unless you have strong justification for archive
3. **Document in final report**:

   ```markdown
   ## Temporary Documentation Created

   File: path/to/TEMP_DOC.md Type: <type> Cleanup Trigger: <when to delete> Cleanup Action: delete|archive|extract Justification: <why this action?>
   ```

4. **Never commit without frontmatter** - orchestrator will reject it
5. **Keep docs minimal** - only what's needed during active task

### Example Sub-Agent Report

```markdown
## Task Complete

### Temporary Documentation Created

**File**: scripts/pypis_delivery_service/IMPLEMENTATION_BUGS.md **Type**: bugs **Cleanup Trigger**: When all bugs fixed and tests pass **Cleanup Action**: delete **Justification**: Bug tracking during active development only. Once fixed, no historical value.

**Orchestrator Action Required**: Delete this file when marking task complete.
```

## Orchestrator Responsibilities

When receiving a sub-agent report mentioning temporary docs:

1. **Add cleanup to todos**:

   ```python
   TodoWrite([
       {"content": "Fix bugs in IMPLEMENTATION_BUGS.md", "status": "pending"},
       {"content": "Delete IMPLEMENTATION_BUGS.md (temp doc cleanup)", "status": "pending"},
   ])
   ```

2. **Validate frontmatter** before accepting temp doc
3. **Execute cleanup** before marking master task complete
4. **Ask user** if archive action recommended (rare)

## Audit Mode

When run without arguments, audit current directory:

```bash
# Find all potential temporary docs
find . -name "*.md" -type f | grep -E "(BUG|INVESTIGATION|CHECKLIST|FINDINGS|SUMMARY|REPORT)" | grep -v "docs/"

# For each file:
# - Check for YAML frontmatter with temporary: true
# - Verify required keys present
# - Check status field
# - Recommend cleanup if status=completed
```

**Output**:

```markdown
## Temporary Documentation Audit

### Files Needing Cleanup (status=completed)

- scripts/pypis_delivery_service/IMPLEMENTATION_BUGS.md (delete)
  - Created: 2025-10-29
  - Task: Rewrite pypis_delivery_service
  - Action: rm scripts/pypis_delivery_service/IMPLEMENTATION_BUGS.md

### Files Missing Frontmatter (VIOLATION)

- scripts/some_notes.md
  - Action: Add frontmatter or delete if truly temporary

### Active Temporary Docs (OK)

- .claude/tasks/investigation_auth_bug.md (status=active)
  - Cleanup trigger: When auth bug fixed
```

## Decision Matrix

| Type          | Default Action    | Archive Only If                                      |
| ------------- | ----------------- | ---------------------------------------------------- |
| bugs          | delete            | Critical bugs showing design flaws worth documenting |
| investigation | extract or delete | Never - put insights in code                         |
| checklist     | delete            | Never                                                |
| findings      | delete            | Never                                                |
| summary       | delete            | Never                                                |
| report        | delete            | Never                                                |

**Rule of thumb**: If you're unsure → **delete**

## Anti-Patterns to Prevent

❌ **Don't do this**:

```markdown
# Implementation Bugs

1. Bug in list iteration
2. Bug in edge case
3. Bug in templates

[No frontmatter, vague descriptions, no cleanup plan]
```

✅ **Do this**:

```yaml
---
temporary: true
type: bugs
task: "Rewrite pypis_delivery_service config system"
agent: "python-pytest-architect"
created: 2025-10-29
cleanup_trigger: "All 3 bugs fixed and tests pass"
cleanup_action: delete
status: active
---
# Implementation Bugs

## Bug #1: List iteration [FIXED: abc123]
## Bug #2: Edge case [IN PROGRESS]
## Bug #3: Templates [PENDING]
```

## Git Integration

**During active development**:

```bash
# OK to commit temp docs with proper frontmatter
git add IMPLEMENTATION_BUGS.md
git commit -m "docs(temp): track implementation bugs during refactoring"
```

**At task completion**:

```bash
# Delete temp docs (default)
git rm IMPLEMENTATION_BUGS.md
git commit -m "docs: remove temporary bug tracking (all resolved)"
```

## Available Script Commands

### `scan` - Progressive Disclosure (Orchestrator Primary Command)

```bash
# Minimal summary (50-100 tokens)
.claude/utilities/find-temp-documentation.py scan --summary

# Full list with indices (200-500 tokens)
.claude/utilities/find-temp-documentation.py scan

# Grouped view
.claude/utilities/find-temp-documentation.py scan --group-by type
```

**Output Format**: YAML (easy for Claude to parse)

### `show` - View Specific Document

```bash
# Frontmatter only (100-150 tokens)
.claude/utilities/find-temp-documentation.py show --index 0 --frontmatter
.claude/utilities/find-temp-documentation.py show implementation-bugs-pypis-config --frontmatter

# Full content
.claude/utilities/find-temp-documentation.py show --index 0
.claude/utilities/find-temp-documentation.py show path/to/doc.md
```

### `state` - Efficient Status Management

```bash
# Get current status
.claude/utilities/find-temp-documentation.py state --index 0

# Set status (no Read+Edit needed)
.claude/utilities/find-temp-documentation.py state --index 0 completed
.claude/utilities/find-temp-documentation.py state implementation-bugs for-review
```

**Status Values**:

- `active` - Currently in use
- `for-review` - Sub-agent finished, orchestrator needs to review
- `completed` - Orchestrator reviewed, ready for cleanup
- `archived` - Historical value (moved to docs/history/)

### `cleanup` - Safe Removal to .wastepaper/

```bash
# Dry-run (see what would be moved)
.claude/utilities/find-temp-documentation.py cleanup

# Actually move files
.claude/utilities/find-temp-documentation.py cleanup --no-dry-run

# Clean specific task
.claude/utilities/find-temp-documentation.py cleanup --task-file .claude/tasks/task.md --no-dry-run
```

**Safety**: Files moved to `.wastepaper/` (gitignored), not deleted. Recoverable if needed.

### `audit` - Find Violations

```bash
.claude/utilities/find-temp-documentation.py audit
```

Finds files matching temporary doc patterns but missing valid YAML frontmatter.

## Verification Before Task Completion

Before marking ANY task complete, verify:

```bash
# Check summary - should show 0 active or for-review docs for this task
.claude/utilities/find-temp-documentation.py scan --summary

# Cleanup completed docs
.claude/utilities/find-temp-documentation.py cleanup --no-dry-run

# Verify cleanup worked
.claude/utilities/find-temp-documentation.py scan --summary
# Should show: total: 0 (or only active docs from other tasks)
```

**Rules**:

- If any docs have `status: active` → TASK NOT COMPLETE YET
- If any docs have `status: for-review` → ORCHESTRATOR MUST REVIEW FIRST
- If any docs have `status: completed` → RUN CLEANUP
- Only mark task complete when all temp docs cleaned up

## Expected Output

When running this command:

```markdown
## Temporary Documentation Lifecycle Results

### Action Taken: <delete|archive|extract|audit>

### Files Processed:

- scripts/pypis_delivery_service/IMPLEMENTATION_BUGS.md
  - Status: completed → deleted
  - Commit: abc123 "docs: remove temporary bug tracking"

### Remaining Active Temp Docs: 0

### Task Complete: Yes
```

## Integration

**Prerequisites**: Sub-agent has completed its task **Follow-up**: Mark master task as complete **Related**: `/cleanup-session`, `/archive-task`

## Example Usage

```bash
# Cleanup specific temp doc
/temp-doc scripts/pypis_delivery_service/IMPLEMENTATION_BUGS.md

# Audit current directory
/temp-doc

# Audit specific directory
/temp-doc scripts/pypis_delivery_service/
```
