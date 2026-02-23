# Effectiveness Rubric

Loaded during Phase 2. Defines improvement types, severity calibration, and
effort/impact criteria used to rank findings in Phase 3.

---

## Improvement Types

### NEW FEATURES
Capabilities the skill doesn't have but should. The user needs to do something manually today,
or the skill stops before the user's actual goal is accomplished.

Examples:
- Skill generates a plan but doesn't apply it — user applies manually
- Skill processes one file at a time but the user almost always has multiple
- Skill never verifies its own factual claims against current docs

### UX IMPROVEMENTS
Things that work correctly but create friction. The outcome is right; the path to it is rough.

Examples:
- Necessary question asked in Phase 3 after heavy work is complete (should be Phase 1)
- Skill produces dense prose where the user will need to copy individual sections
- No explicit path for "I don't know what I want" — skill assumes a specific user goal
- Consequential decision made silently where user would want to choose

### ACCURACY FIXES
Factual claims in the skill that have drifted from current reality. The skill is confidently
wrong about something external — a tool name, an API parameter, a frontmatter field.

Examples:
- Frontmatter field listed that no longer exists or has been renamed
- Model ID that has been deprecated
- Third-party CLI flag that changed between versions

### EFFICIENCY GAINS
The same outcome is achievable with less token cost, fewer tool calls, or fewer user
interactions — without changing the output quality.

Examples:
- Three sequential Read calls that could be one
- A tool call where inline `@path` injection would work
- An agentic reasoning step that always produces the same structure (script candidate)

---

## Severity Calibration

Severity is about impact on the user when the skill is triggered — not how hard the fix is.

### HIGH VALUE
Apply this label when the absence of the improvement means the skill fails the user in a
common scenario. The user either can't complete their goal, gets wrong output they'll trust,
or has to do significant manual work that the skill should have done.

Signs: blocks the main use case, produces wrong output silently, creates a dead end after
a long workflow, applies to every invocation not just edge cases.

### MEDIUM VALUE
Apply when the improvement makes the skill meaningfully better for a common scenario, but
the current behavior produces a usable (if degraded) outcome. The user can work around it.

Signs: adds a capability users would reach for regularly, removes friction from the main
path, fixes docs that are outdated but the user could discover the right value independently.

### LOW VALUE
Apply for polish — the skill works well, this makes it slightly better. Affects edge cases
or infrequent scenarios.

Signs: applies only to a specific edge case, or the fix is cosmetic with no behavioral change.

---

## Effort Estimation

Use effort estimates to help the user decide what to apply first. Three tiers:

**One instruction change** — Edit one sentence or paragraph in SKILL.md. No new phases, no
new files. Examples: adding an explicit full-audit path, clarifying an underspecified step,
updating a stale field name.

**New phase or section** — Add a meaningful new workflow phase with its own entry condition
and exit condition. Requires ~10–30 lines of SKILL.md. Examples: adding live doc validation
as a distinct phase, adding a post-apply explanation phase.

**New script** — Requires creating a file in `scripts/`, wiring it into SKILL.md with
trigger condition, exact invocation, and output handling. Higher effort but produces
deterministic, reliable behavior for consistency-critical steps.

---

## The Stuck / Diverge / Dead-End / Friction Framework (2a)

These four categories from mental simulation map to improvement types:

| Finding type | Maps to | Why |
|---|---|---|
| **Stuck point** | NEW FEATURE or UX IMPROVEMENT | Skill needs more info it never asks for |
| **Divergence point** | UX IMPROVEMENT or EFFICIENCY GAIN | Underspecified step → inconsistent output |
| **Dead end** | NEW FEATURE | Skill stops; user goal is not met |
| **Friction point** | UX IMPROVEMENT | Workflow interrupts user at wrong moment |

---

## Quick Reference: High-Signal Improvements

These patterns appear often and are almost always high-value:

1. **No "I don't know" path** — Any skill that starts by asking for user intent but offers no
   fallback for an uncertain user has a gap. Add an explicit full-audit mode.

2. **Factual claims without verification** — Any skill that references versioned external
   standards (frontmatter fields, API parameters, tool names) without checking them will
   eventually produce wrong output. Live doc validation is almost always worth adding.

3. **Workflow stops at proposal** — If the skill's final output is a plan, checklist, or
   report but the user obviously wants those things *applied*, that's a dead end gap.

4. **Mid-workflow questions** — Any question asked after Phase 1 that could have been asked
   in Phase 1 is friction. If the answer is needed to do the work, it belongs at the start.

5. **Silent consequential decisions** — If the skill makes a significant choice (what model
   to use, what to include in output, what to skip) without surfacing it, that's a divergence
   point. Either make it explicit or give the user control.
