# Plugin Validator Findings Investigation

## Summary

- **SK005** — Mixed. The validator logic is sound (it is a warning, not an error), but the required trigger phrase list is too narrow. Several affected descriptions convey intent clearly without matching the exact literal phrases. At least 3 of the 9 files examined have descriptions that communicate trigger context effectively — the rule should accept a wider vocabulary or use semantic matching.
- **SK006** — Real Issue. Token counts above 4400 indicate genuinely large skills. The warning is appropriate given the validator's documented rationale (comparing against Anthropic's official skill sizes). No linter bug.
- **SK007** — Real Issue for `orchestrating-swarms` (12453 tokens). For `claude-hooks-reference-2026` (9054 tokens), it is borderline — the content is intentional reference material, but the threshold is still architecturally valid.
- **LK001** — Mixed. `perl-testing/references/mock-patterns.md`, `brainstorming-skill/references/geeky-gadgets.md`, and `agent-browser` references/templates: **real content issues** (files do not exist). The `.claude/skills/*/` links using `./../knowledge/...` pattern: **real content issues** caused by incorrect path depth — the files exist but the links resolve one level too shallow. The `claude-hooks-reference-2026` link to `../claude-skills-reference-2026/SKILL.md`: **real content issue** — the directory was renamed to `claude-skills-overview-2026`.
- **LK002** — Real Issue enforcing a style preference that matters. Links without `./` prefix resolve correctly in practice, but the validator correctly enforces the stated project convention. The `agentskills` and `agent-browser` links that lack `./` also lack the referenced files entirely, making them LK001 errors too.

**Overall: 2 real issues, 2 mixed, 1 borderline — no outright linter bugs.**

---

## SK005 — Description Missing Trigger Phrases

### Files Examined

1. `plugins/python3-development/skills/development/implement-feature/SKILL.md`
   ```text
   description: Execute a SAM task plan (plan/tasks-*.md) by looping ready tasks,
   delegating each task to its specified agent, and relying on hooks to update task
   timestamps/status.
   ```

2. `plugins/plugin-creator/skills/refactor-plugin/SKILL.md`
   ```text
   description: Start a complete plugin refactoring workflow. Analyzes plugin structure,
   creates refactoring plan with tasks, and guides through execution.
   ```

3. `.claude/skills/work-backlog-item/SKILL.md`
   ```text
   description: "Bridges BACKLOG.md to the SAM planning pipeline — no args shows
   interactive backlog browser with grooming status; with args finds item by title
   substring, auto-grooms if needed, runs RT-ICA to BLOCK on missing inputs before
   SAM planning, invokes add-new-feature, then updates backlog with plan reference.
   STOPS if item already has a Plan field or RT-ICA returns BLOCKED."
   ```

4. `plugins/plugin-creator/skills/arl/SKILL.md`
   ```text
   description: Knowledge reference for Autonomous Refinement Loop research — pattern
   research into prerequisites for autonomous execution without synchronous human
   blocking gates.
   ```

5. `.claude/skills/commit-staged/SKILL.md`
   ```text
   description: Generate descriptive commit messages by analyzing git diffs, very fast
   and context pollution safe. Use on any request to commit staged changes.
   ```

6. `.claude/skills/subagent-contract/SKILL.md`
   ```text
   description: Global contract for all specialist subagents. Enforces role boundaries,
   scope discipline, and DONE/BLOCKED status signaling. Load this skill in any agent
   that should operate as a bounded specialist following supervisor delegation patterns.
   ```

### Validator Logic (file:line)

Source: `plugins/plugin-creator/scripts/plugin_validator.py`, lines 68-76 and 1947-1964.

```python
REQUIRED_TRIGGER_PHRASES = [
    "use when",
    "use this",
    "used when",
    "used by",
    "when ",
    "trigger",
    "activate",
]
```

The check (lines 1950-1952) lowercases the description and looks for any of these exact substrings. It fires a **warning** (not an error). The check applies only to `FileType.SKILL` files (line 1948).

### Verdict: Mixed

The rule design is reasonable — trigger phrases make skills more discoverable by the orchestrator. The severity (warning) is appropriate.

However, the trigger phrase list is too narrow:

- `commit-staged` description contains "Use on any request" — this is a clear trigger phrase but does not match any of the 7 literals. This is a **false positive** from the linter.
- `work-backlog-item` description contains "no args shows..." and "with args finds" — behavioral triggers clearly present, but no literal match.
- `subagent-contract` description contains "Load this skill in any agent" — functional trigger, not captured.
- `refactor-plugin` description says "Start a complete plugin refactoring workflow" — imperative trigger implied, not matched.
- `arl` and `implement-feature` descriptions are genuinely weak on trigger context. Real issues.

The rule misses trigger intent expressed through other vocabulary: "load this skill", "on any request", "use on".

### Evidence

- `commit-staged`: Contains `"Use on any request"` — functionally identical to `"use when"` but not matched by the literal check. False positive.
- `subagent-contract`: Contains `"Load this skill in any agent"` — clear load trigger not in phrase list. False positive.
- `arl`: No trigger phrase present. Real issue — description reads as a content summary, not a usage trigger.
- `implement-feature`: No trigger phrase. Real issue — describes what it does, not when to invoke it.

---

## SK006 — Skill Body Large (Warning)

### Verdict: Real Issue

SK006 is a warning at `TOKEN_WARNING_THRESHOLD = 4400` tokens (body content only, frontmatter excluded). The validator comment states this threshold is derived from Anthropic's official skill sizes.

### Evidence

Source: `plugin_validator.py` lines 48-49, 2087-2098.

```python
TOKEN_WARNING_THRESHOLD = 4400
TOKEN_ERROR_THRESHOLD = 8800
```

Warning message: "This skill is larger than Anthropic's official skills. Review whether content can be moved to references/ or if the skill covers multiple domains that could be separated."

This is actionable architectural guidance. Skills exceeding 4400 body tokens should be reviewed for content that can move to `references/` files or split into sub-skills. Not a linter bug — valid warning.

---

## SK007 — Skill Body Exceeds Token Limit (Error)

### Files Affected

- `.claude/skills/orchestrating-swarms/SKILL.md` — 12453 tokens (42% over 8800 limit)
- `plugins/plugin-creator/skills/claude-hooks-reference-2026/SKILL.md` — 9054 tokens (3% over 8800 limit)

### Verdict: Real Issue

Source: `plugin_validator.py` lines 48-49, 2077-2086.

The error threshold `TOKEN_ERROR_THRESHOLD = 8800` represents a hard architectural limit with the suggestion to split. Both files genuinely exceed it.

`orchestrating-swarms` at 12453 tokens is significantly over — clear split candidate.

`claude-hooks-reference-2026` at 9054 tokens is marginally over, and its content is intentional reference documentation. The threshold is still valid: reference material of this density should be distributed across `references/` files instead of embedded in the SKILL.md body. The validator's suggestion ("Run /plugin-creator:refactor-skill") is appropriate.

Neither finding is a linter bug.

---

## LK001 — Broken Internal Links

### Files Examined

**Group A — Genuinely missing files:**

- `plugins/perl-development/skills/perl-testing/SKILL.md` → `./references/mock-patterns.md`
  - Verified: `plugins/perl-development/skills/perl-testing/references/` contains no files (Glob returned no results).
  - Verdict: **Real content issue.** The file was never created or was deleted.

- `plugins/brainstorming-skill/skills/brainstorming-skill/SKILL.md` → `./references/geeky-gadgets.md`
  - Verified: `plugins/brainstorming-skill/skills/brainstorming-skill/references/` contains no files.
  - Verdict: **Real content issue.** The references directory and files do not exist.

- `.claude/skills/agent-browser/SKILL.md` → 9 broken links (references/*.md and templates/*.sh)
  - Verified: `.claude/skills/agent-browser/references/` and `templates/` contain no files (Glob returned no results).
  - Verdict: **Real content issue.** References and templates directories do not exist.

**Group B — Incorrect path depth from `.claude/skills/*/` to `.claude/knowledge/`:**

- `.claude/skills/verify/SKILL.md` → `./../knowledge/workflow-diagrams/master-workflow.md`
- `.claude/skills/subagent-contract/SKILL.md` → `./../knowledge/workflow-diagrams/multi-agent-orchestration.md`
- `.claude/skills/scientific-thinking/SKILL.md` → `./../knowledge/workflow-diagrams/investigation-workflow.md`
- `.claude/skills/agent-creator/SKILL.md` → `./../knowledge/workflow-diagrams/asset-decision-tree.md`

  Path resolution test (Python):
  ```
  skill_dir = {git_root}/.claude/skills/subagent-contract
  link = ./../knowledge/workflow-diagrams/multi-agent-orchestration.md
  resolved = {git_root}/.claude/skills/knowledge/workflow-diagrams/multi-agent-orchestration.md
  exists = False

  correct = {git_root}/.claude/knowledge/workflow-diagrams/multi-agent-orchestration.md
  correct.exists = True
  ```

  The path `./../` from `.claude/skills/subagent-contract/` navigates to `.claude/skills/` (not `.claude/`). To reach `.claude/knowledge/`, the correct path is `../../knowledge/workflow-diagrams/multi-agent-orchestration.md`.

  All four workflow-diagram target files verified to exist at the correct absolute paths.

  Verdict: **Real content issue — incorrect path depth.** The links use one too few `../` levels. The linter correctly reports broken links. Fix: change `./../knowledge/` to `../../knowledge/` in all four files.

**Group C — Renamed/moved directory:**

- `plugins/plugin-creator/skills/claude-hooks-reference-2026/SKILL.md` → `../claude-skills-reference-2026/SKILL.md`
  - Verified: No directory named `claude-skills-reference-2026` exists in `plugins/plugin-creator/skills/`. The directory `claude-skills-overview-2026` exists.
  - Verdict: **Real content issue.** Directory was renamed from `claude-skills-reference-2026` to `claude-skills-overview-2026`. The link was not updated.

### Validator Logic (file:line)

Source: `plugin_validator.py`, lines 614-629.

```python
skill_dir = path.parent
link_path = (skill_dir / link_url_no_fragment).resolve()
if not link_path.exists():
    errors.append(ValidationIssue(..., code=LK001))
```

The validator resolves links relative to the SKILL.md's parent directory, strips anchor fragments before checking, and reports an error if the resolved path does not exist. Logic is correct.

The `_should_ignore_link` method (lines 662-680) ignores HTTP/HTTPS/FTP external links, anchor-only links (`#`), and absolute paths (`/`). It does NOT ignore `../` relative paths. This is correct behavior — `../` links should be validated.

### Verdict: Real Issue (all cases)

No linter bugs found in LK001. All reported broken links reflect genuine missing or mis-pathed files.

---

## LK002 — Link Missing ./ Prefix

### Files Examined

- `plugins/plugin-creator/skills/agentskills/SKILL.md` — 4 links without `./` prefix:
  - `../claude-skills-overview-2026/SKILL.md` (2 occurrences)
  - `references/forms.md`
  - `references/api.md`
  - `references/redlining.md`
  - `references/best-practices.md`
  - `references/specification.md`
  - `references/integration.md`

- `.claude/skills/agent-browser/SKILL.md` — 9 links without `./` prefix (all references and templates)

### Validator Logic (file:line)

Source: `plugin_validator.py`, lines 596-608.

```python
if not link_url.startswith("./") and not link_url.startswith("../"):
    warnings.append(ValidationIssue(..., code=LK002, severity="warning"))
```

Links using `../` are NOT flagged for LK002 — only links starting with neither `./` nor `../`. This means bare relative links like `references/file.md` trigger LK002 but `../sibling/SKILL.md` does not.

### Verdict: Real Issue (as a style enforcement warning)

The severity is **warning** (not error), which is appropriate. Links without `./` prefix will resolve correctly at the filesystem level (path resolution works the same), so these are not broken links. They violate the project convention documented in CLAUDE.md:

> "Use markdown links with relative paths starting with ./"

For `agentskills`, the `references/*.md` links that lack `./` also lack the underlying files — these are already LK001 errors too.

For `agent-browser`, same situation — the 9 links without `./` prefix also point to non-existent files (already LK001 errors).

The LK002 warning is a style/convention enforcement finding and is doing its job correctly. No linter bug.

---

## Recommendations

### SK005 — Expand the trigger phrase list

**Action: Fix the linter.** Add more natural trigger vocabulary to `REQUIRED_TRIGGER_PHRASES` in `plugin_validator.py`:

```python
REQUIRED_TRIGGER_PHRASES = [
    "use when",
    "use this",
    "use on",        # catches "Use on any request"
    "used when",
    "used by",
    "when ",
    "trigger",
    "activate",
    "load this",     # catches "Load this skill in any agent"
    "load when",
    "invoke",
    "invoked",
]
```

**Also fix the files:** `arl/SKILL.md` and `implement-feature/SKILL.md` descriptions lack genuine trigger context and should have trigger phrases added regardless of the phrase list.

### SK006 — No action needed on the linter

Review affected skills individually. Move content to `references/` or split as appropriate. The warning is correctly calibrated.

### SK007 — Fix the files, not the linter

- `.claude/skills/orchestrating-swarms/SKILL.md`: Split into sub-skills or move content to `references/`. 12453 tokens is well above the threshold.
- `plugins/plugin-creator/skills/claude-hooks-reference-2026/SKILL.md`: Move dense reference tables to `references/` files and link from SKILL.md body. 9054 tokens at 3% over threshold — a modest restructuring would resolve it.

### LK001 — Fix the files

Four distinct fixes needed:

1. **perl-testing**: Create `plugins/perl-development/skills/perl-testing/references/mock-patterns.md` or remove the link.
2. **brainstorming-skill**: Create `plugins/brainstorming-skill/skills/brainstorming-skill/references/geeky-gadgets.md` or remove the link.
3. **agent-browser**: Create the 9 referenced `references/*.md` and `templates/*.sh` files or remove/update the links.
4. **verify, subagent-contract, scientific-thinking, agent-creator**: Change `./../knowledge/workflow-diagrams/` to `../../knowledge/workflow-diagrams/` in all four SKILL.md files.
5. **claude-hooks-reference-2026**: Change `../claude-skills-reference-2026/SKILL.md` to `../claude-skills-overview-2026/SKILL.md`.

### LK002 — Fix the files (style)

Add `./` prefix to bare relative links in `agentskills/SKILL.md` and `agent-browser/SKILL.md`. Note that most of these links also need LK001 fixes (missing files) — the `./` prefix fix alone is insufficient.

---

_Investigation date: 2026-02-18. No files modified during investigation._
