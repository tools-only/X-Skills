# Skill Anatomy Reference

Gold standard for what a well-formed skill looks like at each tier. Load this before
running any audit dimension — it is the rubric every dimension evaluates against.

---

## The Three-Level Loading Model

This model governs every progressive disclosure and verbosity decision in a skill.

| Level | What loads | When |
|-------|------------|------|
| **Metadata** | Frontmatter `description` | Always — every session, whether or not the skill triggers |
| **Skill body** | `SKILL.md` (below the `---`) | Only when the skill triggers |
| **Resources** | `references/`, `examples/`, `scripts/` | Only when Claude explicitly reads them |

**Why this matters for every audit decision:**
- Anything in the description costs tokens on *every session* — justify every word
- Anything in SKILL.md costs tokens on *every invocation* — content that isn't needed for most runs belongs in `references/`
- Anything in `references/` or `scripts/` is free until needed — defer detail here aggressively

---

## Gold Standard Directory Anatomy

### Simple skill (no subdirectories needed)

```
skill-name/
└── SKILL.md    # Single file; body under 200 lines
```

Use when: the skill has one clear workflow, no reusable code, no domain-specific reference data.

### Standard skill

```
skill-name/
├── SKILL.md          # Core instructions; 200–400 lines
└── references/
    └── topic.md      # Detail deferred from SKILL.md
```

Use when: SKILL.md would exceed ~300 lines without deferral, or when the skill has sections
only needed for specific sub-tasks (e.g., advanced options, domain-specific lookup tables).

### Complex skill

```
skill-name/
├── SKILL.md          # Navigation + core flow; under 500 lines
├── scripts/          # Deterministic operations
├── references/       # Documentation Claude reads while working
├── examples/         # Complete, runnable artifacts users can copy
└── assets/           # Output-only files (templates, images, fonts)
```

Use when: skill has reusable deterministic operations, domain-specific reference material,
or produces output users will copy and adapt.

---

## Directory Type Definitions

### `scripts/`

Executable code (Python, Bash). The defining characteristic: scripts can run without being
loaded into context — they are invoked by path, not read into the conversation. Use when:

- The same code block appears more than once across invocations (inline = re-generated each time)
- The operation is fragile or must produce identical output consistently
- A utility would benefit multiple phases of the workflow

Scripts that exist but are not referenced in SKILL.md are dead code — Claude cannot use
them without a pointer. Every script must have an explicit reference with a description of
when and how to invoke it.

### `references/`

Documentation Claude reads while working. The defining characteristic: loaded only when
Claude decides it needs them — kept out of SKILL.md to reduce base context cost. Use when:

- A section of SKILL.md would exceed ~100 lines for a single topic
- The information is only needed for specific sub-tasks, not every invocation
- Domain-specific data (lookup tables, option catalogs, field definitions) is needed
- The skill has multiple variants (e.g., AWS vs GCP vs Azure) and each needs its own file

Include a table of contents at the top of any reference file over 100 lines.

### `examples/`

Complete, runnable artifacts. The defining characteristic: users copy these directly without
modification. Distinct from references (docs you read) and scripts (utilities you invoke).
Use when:

- The skill produces output users will adapt (config files, template code, structured prompts)
- Working examples would disambiguate ambiguous instructions more efficiently than prose
- The skill teaches by demonstration and the demo needs to be intact and executable

### `assets/`

Files used in output but never loaded into Claude's context. The defining characteristic:
consumed by other tools or included in deliverables, not read by the model. Use for:
templates, images, fonts, boilerplate files, starter projects.

---

## Degrees of Freedom

Match instruction specificity to the task's fragility and variability. Over-constraining
flexible tasks breaks generalization; under-constraining fragile tasks breaks reliability.

| Level | When to use | Format |
|-------|-------------|--------|
| **High freedom** | Multiple valid approaches; context determines best path | Text heuristics, principles, criteria |
| **Medium freedom** | Preferred pattern exists; some variation is acceptable | Pseudocode, parameterized scripts |
| **Low freedom** | Fragile operations, exact sequence required, consistency critical | Exact scripts, minimal parameters |

**Audit application:** for each agentic step, ask whether the freedom level matches the
task's actual variability. A step described vaguely that must produce consistent output
is under-constrained (needs to move toward low freedom). A step with exact scripts for
a judgment-heavy task is over-constrained (model's generalization is being wasted).

---

## Naming Conventions

| Element | Rule |
|---------|------|
| Skill name | Lowercase letters, digits, hyphens only; max 64 chars |
| Command name | Verb-led: `fix-issue`, `review-pr`, `deploy-staging` |
| Namespace | Use tool prefix when it aids routing: `gh-address-comments`, `linear-close-issue` |
| Reference files | Lowercase, hyphens, descriptive: `frontmatter-options.md`, `skill-anatomy.md` |
| Script files | Lowercase, hyphens or underscores: `init_skill.py`, `validate-input.sh` |

---

## Body Structure Conventions

### Voice and framing

- **Imperative voice throughout**: "Analyze", "Generate", "Identify" — not "You should analyze"
- **No first-person narrative**: Never "I will", "I am", "I'll then check" — the skill is instructions, not a plan
- **No second-person**: Avoid "you" entirely — second-person reads as addressing the user, not Claude
- **No hedging language**: Remove "you might want to", "consider possibly", "generally speaking"

### Header depth

- H2 (`##`) for major phases or top-level sections
- H3 (`###`) for sub-topics within a phase
- No H4+ (`####`) — content this granular belongs in `references/`, not SKILL.md

### Routing guidance placement

- **Description only** — routing guidance ("when to use this skill") belongs exclusively in frontmatter
- **Never in the body** — the body loads only after the skill has already triggered; routing guidance there is never read by the routing decision and wastes context on every invocation

### Size invariants

| SKILL.md | Interpretation |
|----------|----------------|
| Under 200 lines | Simple skill — likely no `references/` needed |
| 200–400 lines | Standard — check if any sections are invocation-frequency-selective |
| 400–500 lines | Approaching limit — audit for deferrable content |
| Over 500 lines | Requires `references/` deferral; content density is harming maintainability |

---

## Gap Analysis Checklist

Use this to identify what a skill *would benefit from* adding, not just what's wrong.

**Would a `scripts/` directory help?**
- [ ] Is there a code block that appears or would appear more than once?
- [ ] Is there a fragile operation that must produce consistent output?
- [ ] Is there a setup/validation/cleanup step that could be scripted?

**Would a `references/` directory help?**
- [ ] Does SKILL.md exceed 300 lines and have sections only needed for specific sub-tasks?
- [ ] Does the skill have domain-specific lookup data (option catalogs, field tables)?
- [ ] Does the skill have multiple variants (per-provider, per-framework) that could be split?

**Would an `examples/` directory help?**
- [ ] Does the skill produce output users will adapt (configs, templates, prompts)?
- [ ] Are there ambiguous instructions that a complete working example would clarify better than prose?
- [ ] Does the skill teach by demonstration?

**Would additional frontmatter fields help?**
- [ ] Does the skill accept a file path argument but lack `argument-hint`?
- [ ] Does the skill do heavy analysis that would benefit from `model: opus`?
- [ ] Does the skill produce heavy output that would pollute main context (`context: fork`)?
- [ ] Does the skill need user decisions mid-workflow but lacks `AskUserQuestion` in `allowed-tools`?
