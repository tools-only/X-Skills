# CLAUDE.md Rules Extraction Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a rules-extraction phase to the `optimize-claude-md` skill so that when optimizing a CLAUDE.md file, the `contextual-ai-documentation-optimizer` agent also identifies path/filetype-conditional sections and migrates them atomically into `.claude/rules/*.md` files with correct `paths` frontmatter.

**Architecture:** Two files change. The skill (`SKILL.md`) gains a new CLAUDE.md file-type strategy block and a reference to a new reference file. The reference file (`references/claude-rules-extraction.md`) contains the `.claude/rules/` specification — `paths` frontmatter syntax, glob patterns, detection heuristics, filename conventions, and extraction/replacement rules. The agent reads the reference file on demand when processing a CLAUDE.md target and performs the extraction atomically (writes rules files + modifies CLAUDE.md in one step, then shows diff).

**Tech Stack:** Markdown, YAML frontmatter, Claude Code `.claude/rules/` system (paths-scoped rules)

---

### Task 1: Create the reference file

**Files:**
- Create: `.claude/skills/optimize-claude-md/references/claude-rules-extraction.md`

**Step 1: Write the reference file**

Create `.claude/skills/optimize-claude-md/references/claude-rules-extraction.md` with this exact content:

```markdown
# Claude Rules Extraction Reference

## Table of Contents

1. [What is `.claude/rules/`](#what-is-clauderules)
2. [Rules file format](#rules-file-format)
3. [Path glob syntax](#path-glob-syntax)
4. [Detection heuristics](#detection-heuristics)
5. [Filename conventions](#filename-conventions)
6. [What to extract vs what to keep](#what-to-extract-vs-what-to-keep)
7. [Extraction procedure](#extraction-procedure)
8. [Replacement stub format](#replacement-stub-format)

---

## What is `.claude/rules/`

`.claude/rules/` is a Claude Code feature for modular, topic-specific project instructions. All `.md` files in `.claude/rules/` are automatically loaded as project memory. Files with a `paths` YAML frontmatter field are **conditional** — they only apply when Claude is working with files matching the specified glob patterns.

Canonical location: `.claude/rules/<slug>.md`

User-level rules (apply across all projects): `~/.claude/rules/<slug>.md`

SOURCE: https://code.claude.com/docs/en/memory.md (accessed 2026-02-24)

---

## Rules file format

```markdown
---
paths:
  - "src/**/*.ts"
  - "lib/**/*.ts"
---

# TypeScript Rules

- All functions must have explicit return types
- Use `unknown` instead of `any`
```

**Frontmatter fields:**
- `paths` (optional): YAML list of glob patterns. When present, rules only activate for matching files. When absent, rules apply unconditionally.

**No other frontmatter fields are used** in rules files. Do not add `name`, `description`, `model`, or other skill frontmatter fields.

---

## Path glob syntax

| Pattern | Matches |
|---------|---------|
| `**/*.ts` | All TypeScript files in any directory |
| `src/**/*` | All files under `src/` |
| `*.md` | Markdown files in project root only |
| `src/components/*.tsx` | React components in a specific directory |
| `**/*.{ts,tsx}` | Both `.ts` and `.tsx` files (brace expansion) |
| `{src,lib}/**/*.ts` | TypeScript in either `src/` or `lib/` |
| `.github/workflows/**/*.yml` | All GitHub Actions workflow files |
| `**/scripts/**` | Any file under any `scripts/` directory |
| `pyproject.toml` | Specific file in project root |
| `**/*.py` | All Python files |

Brace expansion is supported. Multiple patterns are OR-combined (file matches if ANY pattern matches).

---

## Detection heuristics

A CLAUDE.md section is a candidate for extraction when it exhibits EITHER signal:

### Signal A — Content language

The section body contains phrases that scope rules to a specific file type, language, tool, path, or context. Examples:

- "When working with `*.py` files"
- "For TypeScript" / "In TypeScript code"
- "Python scripts in `plugins/**/scripts/`"
- "When the file is a `.yml` workflow"
- "In `src/api/`" / "Files under `src/`"
- "When editing GitHub Actions"
- "For scripts with shebangs"
- "`**/*.{ts,tsx}` files"

### Signal B — Heading signals

The section heading names a filetype, language, tool, or path area. Examples that qualify:

- `## Python Development Routing`
- `## CI Workflow Modification Protocol`
- `## GitHub Actions CI Workflow`
- `## TypeScript Standards`
- `## API Development Rules`
- `## Script Invocation`

### Disqualifying check (apply after either signal fires)

Even when a signal fires, DO NOT extract if the section content is **universally applicable** — i.e., the rules apply regardless of which file is being edited. Both signals must agree on scope. If the heading names Python but the rules inside apply everywhere, do not extract.

Examples of sections that should NOT be extracted despite heading signals:
- `## GitHub CLI (gh) Usage` — rules about using the `gh` CLI apply when doing any task, not just when editing `.github/` files
- `## Epistemic Identity` — applies universally
- `## Markdown Formatting Standards` — applies to all files Claude edits

---

## Filename conventions

Derive the rules filename from the section heading slug:

| Section heading | Filename |
|-----------------|----------|
| `## Python Development Routing` | `python-development.md` |
| `## CI Workflow Modification Protocol` | `ci-workflows.md` |
| `## GitHub Actions CI Workflow` | `github-actions.md` |
| `## TypeScript Standards` | `typescript.md` |
| `## Script Invocation` | `script-invocation.md` |
| `## API Development Rules` | `api-development.md` |

Rules:
- Lowercase only
- Hyphens between words, no underscores
- Drop generic suffixes: "Protocol", "Standards", "Rules", "Guidelines" → shorten to topic
- Drop "Modification", "Development", "Invocation" when heading has a clearer noun
- Max 30 characters
- Always `.md` extension

---

## What to extract vs what to keep

### Extract

Sections whose rules only apply when working with specific file types, paths, tools, or languages, and where both the heading and content agree on that scope:

- Python-specific coding rules (apply to `**/*.py`)
- CI/CD workflow rules (apply to `.github/workflows/**/*.yml`)
- TypeScript/JavaScript rules (apply to `**/*.{ts,tsx,js}`)
- API endpoint rules (apply to `src/api/**/*`)
- Test file rules (apply to `**/tests/**/*`, `**/*.test.*`)
- Script/shebang rules (apply to `**/scripts/**`)

### Keep in CLAUDE.md

- Identity and core protocol sections
- Universal rules that apply to all files
- Rules about Claude's behavior, communication style, or epistemic stance
- Rules about git, commits, PRs (these apply regardless of file type)
- Delegation and orchestration patterns
- Rules that reference multiple unrelated file types without a clear primary scope
- Short sections (<5 lines of rules) where extraction adds overhead without benefit

---

## Extraction procedure

When a candidate section is identified:

1. **Determine `paths` glob** — derive from content language and heading. Be specific: `**/*.py` not `**/*`. Use brace expansion for multi-extension matches.

2. **Determine filename** — apply filename conventions above.

3. **Write the rules file** at `.claude/rules/<filename>.md`:
   - YAML frontmatter with `paths` list
   - Section heading as H1 title
   - Full section content (preserve all rules, examples, code blocks)

4. **Replace the section in CLAUDE.md** — substitute the extracted section with a compact cross-reference stub (see Replacement stub format below).

5. **Show diff** of both files after writing.

**ATOMIC**: Steps 3 and 4 happen together. Do not write one without the other.

---

## Replacement stub format

Replace the extracted section in CLAUDE.md with:

```markdown
## <Original Heading>

> Extracted to [`.claude/rules/<filename>.md`](.claude/rules/<filename>.md) — applies when working with `<primary glob>`.
```

Keep the heading so the section is still discoverable by heading scan. The stub is one line. Do not preserve any of the original rule content in CLAUDE.md.

---

## CoVe verification questions for extraction

After performing extraction, the CoVe post-check MUST include these additional questions:

1. **Content integrity**: Does the new rules file contain all rules from the original section verbatim?
2. **Stub accuracy**: Does the CLAUDE.md stub reference the correct filename and glob?
3. **Scope correctness**: Does the `paths` glob actually match the files the rules apply to?
4. **No universal rules extracted**: Were any rules that apply to all files accidentally moved to the rules file?
5. **CLAUDE.md coherence**: Does CLAUDE.md still read coherently after the extraction?
```

**Step 2: Verify file written**

Run: `wc -l .claude/skills/optimize-claude-md/references/claude-rules-extraction.md`
Expected: ~140 lines

**Step 3: Commit**

```bash
git add .claude/skills/optimize-claude-md/references/claude-rules-extraction.md
git commit -m "feat(optimize-claude-md): add claude-rules-extraction reference file"
```

---

### Task 2: Update SKILL.md with rules extraction strategy

**Files:**
- Modify: `.claude/skills/optimize-claude-md/SKILL.md`

The current SKILL.md has a `<file_type_strategies>` block inside the delegation template (lines 50–84 in the existing file under `### Phase 3: Delegate to @contextual-ai-documentation-optimizer`). The `**CLAUDE.md:**` bullet currently reads:

```
**CLAUDE.md:** Front-load identity and constraints; use Mermaid flowcharts for decision logic; compress verbose sections using TRIGGER->PROCEDURE->OUTPUT format; minimize content and run the plugin validator after writing to check token complexity; check for behavioral instructions that could be hooks.
```

This needs to expand into a multi-line strategy block that includes the rules extraction phase.

**Step 1: Read the current SKILL.md**

Read `.claude/skills/optimize-claude-md/SKILL.md` in full to locate the exact text to replace.

**Step 2: Replace the CLAUDE.md strategy bullet**

Find this exact text (it appears inside the `<file_type_strategies>` block in the delegation template):

```
**CLAUDE.md:** Front-load identity and constraints; use Mermaid flowcharts for decision logic; compress verbose sections using TRIGGER->PROCEDURE->OUTPUT format; minimize content and run the plugin validator after writing to check token complexity; check for behavioral instructions that could be hooks.
```

Replace it with:

```
**CLAUDE.md:** Front-load identity and constraints; use Mermaid flowcharts for decision logic; compress verbose sections using TRIGGER->PROCEDURE->OUTPUT format; minimize content and run the plugin validator after writing to check token complexity; check for behavioral instructions that could be hooks.

Additionally, run the **Rules Extraction Phase** for CLAUDE.md targets:

1. Read `./references/claude-rules-extraction.md` for the full extraction spec before proceeding.
2. Scan every section for extraction candidates using both detection signals (content language AND heading signals) defined in the reference.
3. Apply the disqualifying check — skip sections whose content is universally applicable despite a scoped heading.
4. For each confirmed candidate: derive `paths` glob and filename per the reference conventions.
5. Write the `.claude/rules/<filename>.md` file with correct `paths` frontmatter and full extracted content.
6. Replace the extracted section in CLAUDE.md with the compact stub format from the reference.
7. Steps 5 and 6 are ATOMIC — write both files before reporting.
8. Show a unified diff of all changed files after extraction.
9. CoVe post-check MUST include the 5 extraction-specific verification questions from the reference.
```

**Step 3: Also update the delegation template TARGET section**

In the same delegation template block (lines 51–94 of SKILL.md), find:

```
TARGET: {resolved path(s)}
FILE TYPE: {CLAUDE.md | SKILL.md | agent definition | reference file}
```

No change needed here — the file type already covers CLAUDE.md.

Find the CONSTRAINTS section in the delegation template. It currently ends with:

```
- Signal DONE when optimization complete, BLOCKED when missing required inputs
```

Add one line before that final line:

```
- For CLAUDE.md: read `./references/claude-rules-extraction.md` before analyzing; perform rules extraction phase after optimization analysis, before CoVe
```

So the CONSTRAINTS block ends with:

```
- For CLAUDE.md: read `./references/claude-rules-extraction.md` before analyzing; perform rules extraction phase after optimization analysis, before CoVe
- Signal DONE when optimization complete, BLOCKED when missing required inputs
```

**Step 4: Verify the SKILL.md reads coherently**

Read the modified SKILL.md in full. Confirm:
- The `<file_type_strategies>` block contains the expanded CLAUDE.md strategy with extraction phase
- The delegation template CONSTRAINTS block includes the reference file instruction
- No other sections were accidentally modified

**Step 5: Run linting**

```bash
uv run prek run --files .claude/skills/optimize-claude-md/SKILL.md
```

Expected: pass (markdown linting)

**Step 6: Commit**

```bash
git add .claude/skills/optimize-claude-md/SKILL.md
git commit -m "feat(optimize-claude-md): add rules extraction phase to CLAUDE.md strategy"
```

---

### Task 3: Validate skill token complexity

**Files:**
- Read: `.claude/skills/optimize-claude-md/SKILL.md` (validator reads it)

**Step 1: Run the plugin validator**

```bash
uv run plugins/plugin-creator/scripts/plugin_validator.py .claude/skills/optimize-claude-md/
```

Expected: No SK006 or SK007 warnings. If SK006 fires (token warning), the added content may need to be trimmed or moved entirely into the reference file.

**Step 2: If SK006 fires**

Move the numbered extraction phase steps (steps 1–9 added in Task 2) out of the SKILL.md delegation template and into the reference file as a new section `## Integration into delegation template`. Replace steps 1–9 in SKILL.md with a single line:

```
Follow the extraction procedure in `./references/claude-rules-extraction.md`.
```

Re-run validator to confirm warning clears.

**Step 3: Commit if changes made**

```bash
git add .claude/skills/optimize-claude-md/SKILL.md .claude/skills/optimize-claude-md/references/claude-rules-extraction.md
git commit -m "fix(optimize-claude-md): move extraction steps to reference to reduce token count"
```

---

### Task 4: Smoke test on the project CLAUDE.md

**Files:**
- Read: `.claude/CLAUDE.md` (target for smoke test)

**Step 1: Identify extraction candidates manually**

Read `.claude/CLAUDE.md`. Identify which sections you expect the agent to flag as extraction candidates based on the heuristics. Write them down as a checklist before running anything. This is your ground truth for verifying agent output.

Sections likely to qualify (verify against actual CLAUDE.md content):
- `## Python Development Routing` (content references `**/*.py`, delegation to Python agents)
- `## CI Workflow Modification Protocol` (content references `.github/workflows/`)
- `## GitHub Actions CI Workflow` (if present as separate section)
- `## Script Invocation` (content references `**/scripts/**`)

Sections likely NOT to qualify:
- `## Epistemic Identity`
- `## Built-in Tools vs Bash Equivalents`
- `## GitHub CLI (gh) Usage` (applies universally)
- `## Verification and Trust Building`

**Step 2: Invoke the skill**

Run: `/optimize-claude-md .claude/CLAUDE.md`

**Step 3: Verify agent output against checklist**

Confirm:
- Agent read `./references/claude-rules-extraction.md` (should be visible in tool use)
- Agent flagged the expected candidates and skipped the expected non-candidates
- For each extracted section: rules file written to `.claude/rules/`, stub in CLAUDE.md is correct
- Diff shown covers both the new rules files and the modified CLAUDE.md
- CoVe includes the 5 extraction-specific questions

**Step 4: If agent missed candidates or extracted incorrectly**

Do NOT commit the smoke test output. Identify which heuristic failed (Signal A, Signal B, or disqualifying check). Update the reference file to sharpen the heuristic. Re-run. Repeat until output matches checklist.

**Step 5: Revert smoke test changes**

```bash
git checkout -- .claude/CLAUDE.md
git rm -f .claude/rules/*.md 2>/dev/null || true
```

The smoke test is validation only — do not commit CLAUDE.md extraction to this branch.

---

## Definition of Done

- [ ] `.claude/skills/optimize-claude-md/references/claude-rules-extraction.md` exists and contains full spec
- [ ] `.claude/skills/optimize-claude-md/SKILL.md` delegation template includes rules extraction phase for CLAUDE.md
- [ ] Plugin validator reports no SK006/SK007 on the skill
- [ ] Smoke test: agent correctly identifies extraction candidates and non-candidates
- [ ] Smoke test: agent writes valid `.claude/rules/` files with correct `paths` frontmatter
- [ ] Smoke test: agent replaces extracted sections with correct stubs
- [ ] Smoke test: CoVe includes extraction-specific verification questions
- [ ] All commits follow `feat(optimize-claude-md):` / `fix(optimize-claude-md):` scope convention
