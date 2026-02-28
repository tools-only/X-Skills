---
name: improve-skill
description: >
  This skill should be used when the user asks to "improve a skill", "make this skill better",
  "add features to a skill", "this skill is missing something", "upgrade my skill", "what's
  missing from this skill", "improve how the skill works step by step", "the skill doesn't do X",
  "make this more useful", "the skill isn't working well for me", "can this skill do more", or
  wants to improve skill effectiveness rather than structural correctness.
argument-hint: <path-to-skill-directory-or-SKILL.md>
model: opus
allowed-tools: Read, Glob, Edit, Write, WebSearch, WebFetch, Task, AskUserQuestion
---

# Skill Improver

Increase the effectiveness of an existing skill by modeling user intent, testing the skill
against that intent through mental simulation and live doc validation, and proposing ranked
improvements — new features, UX gains, accuracy fixes, efficiency wins.

This skill is about whether the skill *accomplishes what users need*. Structural correctness
(malformed frontmatter, missing fields, violated conventions) is repair-skill's domain — if
obvious structural issues are present, note them briefly, recommend `repair-skill`, and
continue with effectiveness analysis.

## Phase 0: Load the Skill

Read `$ARGUMENTS` as the path to a skill directory or SKILL.md file.

- If a directory: read `SKILL.md`, then note which of `references/`, `scripts/`, `examples/`,
  `assets/` exist
- If a file: read it directly, then discover sibling resource directories

If the path is missing or ambiguous, use AskUserQuestion to resolve before proceeding.

Note any obvious structural issues in one sentence ("this description is first-person —
recommend repair-skill for structural fixes") and move on. Do not run a full structural audit.

Phase 0 is complete when SKILL.md is loaded and the sibling directory inventory is noted.

## Phase 1: Understand User Intent

Before analyzing, establish what the user wants. Use AskUserQuestion:

- "What specifically does this skill not do well?"
- Offer options: a specific gap they've noticed, "I'm not sure — run a full effectiveness
  audit", "it works but I want new capabilities", "the workflow feels clunky"

Regardless of the answer, also infer the skill's purpose from its description and body. State
your understanding of what problem it solves and for whom — one sentence — before proceeding.
This grounds the entire analysis in the correct frame.

If the user named a specific complaint: orient the analysis toward that area and scan for
related issues in the same workflow region. If the user is unsure: run the full Phase 2 audit
across all four sub-analyses.

Phase 1 is complete when user intent is established and the skill's purpose is stated.

## Phase 2: Effectiveness Analysis

Load `${CLAUDE_PLUGIN_ROOT}/skills/improve-skill/references/effectiveness-rubric.md` before
starting. It contains the severity framework, improvement type definitions, and effort/impact
calibration criteria used in Phase 3.

Run all four sub-analyses. For each finding record: what the issue is, why it reduces
effectiveness for the user, and the concrete improvement.

---

### 2a — Mental Simulation

Walk through the skill as Claude executing it with a concrete representative user request.
Choose an input that exercises the main workflow path — not an edge case, but the typical use.

For each phase of the skill, evaluate:

- **Missing info**: What does this step need that hasn't been provided or gathered yet? If the
  user hasn't specified something and the skill doesn't ask, what does Claude have to guess?
- **Divergence points**: Where would two different Claude instances execute this instruction and
  arrive at meaningfully different outputs? These are underspecified steps.
- **Dead ends**: Where does the skill's workflow stop but the user's actual goal isn't yet
  accomplished? What does the user have to do manually after the skill finishes?
- **Friction**: Where does the skill pause the user at a low-value decision point (Claude could
  make a good default), or skip user input at a high-value moment (user has a strong preference)?

Document findings by type: stuck points, divergence points, dead ends, friction points.

---

### 2b — Live Doc Validation

Identify all factual claims in the skill that reference external standards: frontmatter field
names, Claude tool names and behavior, API parameters, CLI flags, third-party service interfaces.

For each claim, verify against current documentation:

```
For Claude-specific claims (frontmatter options, tool names, model IDs):
  Use Task tool with subagent_type=claude-code-guide — faster and more accurate than WebSearch.

For third-party claims (npm packages, APIs, external CLIs):
  Use WebSearch or WebFetch against official documentation.
```

Flag drift between what the skill states and current reality. Severity: high if the claim
produces broken output; medium if it produces outdated guidance; low if it's a naming change
with no behavioral difference.

---

### 2c — Feature Adjacency Scan

Given the skill's purpose, identify capabilities that are absent but would be high-value:

- **Adjacent**: Naturally extends what the skill already does — same domain, one step further.
  The user who just ran this skill would almost certainly want this next.
- **Complementary**: Commonly needed right before or after this skill. The user does this
  manually today.
- **End-to-end gap**: The skill starts a job the user finishes by hand — workflow stops at
  "here's a plan" when the user wanted "and apply it."

For each candidate: estimate implementation effort (one instruction change, new phase, or new
script) and user value (rare edge case, common scenario, or blocks the skill in a key scenario).

---

### 2d — UX Flow Review

Evaluate the skill's interaction design:

- Does the skill surface all necessary questions at the start (before heavy work), or does it
  interrupt mid-workflow with requests the user didn't anticipate?
- Is the "I don't know" path explicit? If a user triggers the skill without a specific
  complaint, does the skill handle that gracefully, or does it assume the user knows?
- Does the output format match how users consume it? A report users read once can be dense
  prose; one they apply iteratively needs more structure.
- Are there steps where the skill makes a consequential decision without user input?

Phase 2 is complete when all four sub-analyses are finished and findings are recorded.

## Phase 3: Improvement Proposal

Load `${CLAUDE_PLUGIN_ROOT}/skills/improve-skill/references/effectiveness-report-template.md`
for the output format before constructing the report. Reference
`${CLAUDE_PLUGIN_ROOT}/skills/improve-skill/examples/sample-analysis.md` to calibrate depth
and specificity if needed.

Present findings grouped by improvement type — users think in terms of outcomes, not audit
dimensions. Each entry must include: the sub-analysis code in brackets, what the gap is, why
it matters to the user, and the specific fix. Calibrate severity using the criteria in
`references/effectiveness-rubric.md`.

Ask: "Apply all improvements? Or select specific ones?"

Phase 3 is complete when the report is delivered and user selection is confirmed.

## Phase 4: Apply Improvements

Apply confirmed items in order: new features → accuracy fixes → UX improvements →
efficiency gains.

For each item:
- State what is being changed and why — reference the effectiveness principle, not "you asked"
- Make the edit or create the file
- Confirm the change integrates cleanly with surrounding content

After applying, briefly explain:
- What was changed and why
- What was added and why
- What was left out and why (effort outweighs benefit, or requires domain knowledge the
  user must supply)
- What was not selected — note it remains available to apply later

**Validation:** After delivering the explanation, re-read the modified SKILL.md in full and
confirm: all selected improvements are present, no surrounding content was inadvertently
altered, phases are still numbered with exit conditions, and all references point to files
that exist.

Phase 4 is complete when all confirmed items are applied, the explanation is delivered, and
the validation pass finds no integration failures.

## Phase 5: Structural Lint

After applying effectiveness improvements, invoke the skill-lint agent for a structural
quality pass:

```
Use Task tool with subagent_type=claude-skills:skill-lint:
"Lint the skill at <path-to-skill-directory>. Auto-apply critical and major fixes, report
minor findings for user decision."
```

Wait for the agent to complete. If it auto-applied structural fixes, note them alongside
the effectiveness changes from Phase 4. If it reports minor findings, present them to the
user.

Phase 5 is complete when the lint agent returns and any user-selected minor fixes are applied.
