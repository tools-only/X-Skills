# Refactor-Skill Alignment Review

**Date:** 2026-02-18
**Skill reviewed:** `/plugin-creator:refactor-skill`
**Target skill analyzed:** `stinkysnake` (5428 tokens, SK006 warning)
**Trigger:** Orchestrator proposed splitting stinkysnake into 3 separate skills; user challenged whether that was the correct action

---

## What the Orchestrator Did

1. Read `stinkysnake/SKILL.md` (833 lines, 5428 tokens — SK006 warning, no SK007 error)
2. Read `refactor-skill/SKILL.md`
3. Proposed splitting stinkysnake into three separate skills:
   - A facade `stinkysnake` meta-skill
   - `stinkysnake-analysis` (Phases 1-3)
   - `stinkysnake-review` (Phases 4-5)
4. When the user challenged this, the orchestrator backtracked and offered `references/` extraction instead

---

## What the Refactor-Skill Actually Says

### The opening definition frames the skill's scope narrowly

From refactor-skill/SKILL.md line 10:

> "At the architecture level, this refactoring decomposes a monolithic LLM skill into smaller, purpose-built, independently invocable skills with well-defined responsibilities, inputs, and outputs."

This definition scopes the skill to a specific action: decomposition into independently invocable skills. The word "monolithic" is key — the skill positions itself as a remedy for skills that are structurally problematic, not skills that are merely large.

### Phase 1 Domain Identification Criteria

From refactor-skill/SKILL.md lines 34-43 (the "Domain Identification Criteria" table):

| Signal | Indicates Separate Skill |
|--------|--------------------------|
| Different tool requirements | `tools` would differ |
| Different invocation triggers | Description keywords diverge |
| Independent use cases | Can be used without the other |
| Different expertise domains | Distinct knowledge areas |
| Section size >200 lines | Too large for single concern |
| Different hook requirements | Lifecycle needs differ |

The skill requires the orchestrator to evaluate these signals **before** proposing a split. The Phase 1 instruction (line 25) is: "ANALYZE the source skill thoroughly" — and the criteria table is meant to guide that analysis to a **finding**, not skip directly to "split it."

### Phase 2 requires presenting a plan before executing

From refactor-skill/SKILL.md line 85:

> "STOP and present plan to user before proceeding."

### Quality standards prohibit over-fragmentation

From refactor-skill/SKILL.md lines 363-366:

> "6. Over-fragment (don't create skills <50 lines)"

Minimum viable skill size is defined at lines 368-375:

> "A skill should have enough substance to be useful alone:
> - At least 50 lines of meaningful content
> - At least 2-3 distinct instructions or rules
> - Clear value proposition in description"

### The SK006 validator suggestion

From `plugin_validator.py` lines 2091-2101, the SK006 suggestion text is:

> "This skill is larger than Anthropic's official skills. Review whether content can be moved to references/ or if the skill covers multiple domains that could be separated"

This text presents two options in priority order:
1. Move content to `references/` (content extraction)
2. Split if multiple domains exist (structural decomposition)

---

## Disconnects Found

### 1. Phase 1 analysis was skipped or superficial

The Domain Identification Criteria table gives six signals to evaluate. Applied to stinkysnake:

- **Different tool requirements**: No. All nine phases use the same Read/Grep/Glob/Bash toolset. Phase 4 delegates to an agent and Phase 8 delegates to another agent, but these are task delegations, not tool divergence within the skill itself.
- **Different invocation triggers**: No. The entire skill is invoked with one trigger: "improve Python code quality for a module." All phases serve that single user intent.
- **Independent use cases**: No. Phases 1 through 9 are a sequential workflow. Phase 4 (plan review) is meaningless without Phase 3 (modernization planning). Phase 8 (test-first) is meaningless without Phase 7 (interface design). None of the proposed splits would be usable alone.
- **Different expertise domains**: Arguable. Static analysis vs. type improvement vs. test-first is a conceptual distinction, but stinkysnake already partitions this into named phases. The phases are steps in one workflow, not separate domains a user would invoke independently.
- **Section size >200 lines**: This signal is present — Phase 3 (Modernization Planning) spans lines 238-406, roughly 170 lines. Phase 2 spans lines 153-234. Some phases are large.
- **Different hook requirements**: No.

Of six signals, only one is genuinely present (section size), and it applies to individual phases rather than the skill as a whole. The orchestrator did not document this analysis before proposing a split.

### 2. The orchestrator applied the split as a default response to SK006

The SK006 warning does not mandate splitting into new skills. The validator's suggestion text (lines 2099-2101 of `plugin_validator.py`) lists `references/` extraction first and skill separation second, conditional on "multiple domains." The orchestrator inverted this order by jumping to decomposition before assessing whether domains were truly separate.

### 3. The proposed split created skills that cannot be used independently

The orchestrator proposed `stinkysnake-analysis` (Phases 1-3) and `stinkysnake-review` (Phases 4-5). Phase 4 (plan review) explicitly requires that Phase 3 has already produced a plan file at `.claude/plans/stinkysnake-plan.md`. The review agent prompt at stinkysnake/SKILL.md line 418 references this file directly. A standalone `stinkysnake-review` skill has no path to produce that artifact — it depends entirely on prior phases. This violates the "Independent use cases" criterion the skill requires before splitting.

### 4. The backtrack also misapplied the skill

When challenged, the orchestrator offered `references/` extraction as an alternative. This is the correct instinct — it matches the SK006 suggestion text. However, the orchestrator should have identified this as the primary candidate during Phase 1 analysis, not discovered it only after user correction. The `references/` path was available from the start.

---

## Root Cause Analysis

The orchestrator deviated for two compounding reasons:

**Reason 1: The refactor-skill conflates "this skill exists" with "splitting is appropriate."**

The skill's opening definition, workflow title ("Refactoring Workflow"), all five phases, and the Splitting Strategies section are written entirely around the assumption that splitting into new skills is the goal. The Phase 1 analysis section contains the Domain Identification Criteria, which could prevent an inappropriate split — but the criteria are presented as a table inside an "ANALYZE" section that reads as preparation for the split, not as a gate that might conclude "no split needed."

There is no explicit branching instruction in Phase 1 of the form: "If none of the domain signals are present, do not proceed to Phase 2. Consider references/ extraction instead." The skill implicitly assumes the analysis will confirm that splitting is warranted.

**Reason 2: The orchestrator pattern-matched SK006 → "run refactor-skill" → "split into skills" without fully applying the domain criteria.**

The SK006 warning says "consider" two options. The orchestrator treated invoking refactor-skill as the decision, when the real decision point is inside Phase 1 of refactor-skill. The orchestrator skipped the discrimination step and went to planning.

---

## Gaps in the Refactor-Skill

### Gap 1: No explicit "no-split" exit path in Phase 1

The skill has no instruction for the case where Phase 1 analysis concludes that domain signals are absent. There is no text equivalent to: "If fewer than 2-3 domain signals are present, the skill is not a candidate for splitting. Recommend `references/` extraction instead and stop." Without this exit path, every invocation of refactor-skill defaults to producing a split plan.

### Gap 2: No decision point distinguishing SK006 (warning) from SK007 (error) scenarios

The ERROR_CODES.md documentation (line 437-454) describes `references/` extraction as the SK006 recommendation, while SK007 (line 473) says "Split skill immediately." The refactor-skill itself contains no reference to these thresholds and no guidance that a warning-level signal warrants a more conservative approach than an error-level signal.

### Gap 3: The "Section size >200 lines" criterion is a weak signal when present alone

A single section over 200 lines is listed as an independent indicator of a separate skill. But a large section in a sequential workflow is better addressed by moving detailed reference material to `references/` than by creating a new independently-invocable skill. The criterion does not distinguish between "large because it covers a distinct domain" and "large because it has detailed reference content that belongs in references/."

### Gap 4: No minimum context check before proceeding

The skill does not include a pre-flight check: "Is this skill genuinely monolithic (covers multiple independent user intents) or merely large (covers one intent with substantial detail)?" This distinction determines whether the refactor-skill is the right tool at all. A large but cohesive skill is not monolithic in the architectural sense the skill's definition uses.

---

## Recommendations

### For improving the refactor-skill

1. **Add a Phase 0: Candidate Assessment** before Phase 1 that evaluates whether the skill is actually a candidate for splitting into new skills. The gate criteria should be explicit: if fewer than 2 of the 6 domain signals are present, output "Not a candidate for skill splitting — recommend references/ extraction" and stop.

2. **Add an explicit `references/` extraction path** as a parallel workflow to the splitting workflow. When the only issue is token count (SK006) and content is cohesive, the correct action is to move detailed reference content (large tables, multi-step examples, comprehensive code samples) into `./references/` files and link them from the main skill. This reduces tokens while preserving the skill's single-purpose design.

3. **Qualify the "Section size >200 lines" criterion** with: "This signal indicates potential for splitting only when the section has a distinct invocation trigger. A large section with detailed reference content is a candidate for `references/` extraction, not skill splitting."

4. **Reference the SK006 vs SK007 threshold distinction** in the skill. Add language such as: "If the skill is at warning level (SK006) but not error level (SK007), exhaust `references/` extraction options before proceeding to skill splitting."

5. **Make Phase 1 a gate, not just preparation.** Add a STOP instruction after the Domain Identification Criteria: "If fewer than 2-3 signals are present, this skill does not meet the criteria for splitting. Document the finding and recommend the appropriate alternative action."

### For improving orchestrator compliance

1. **The Domain Identification Criteria must be evaluated explicitly, not inferred.** When running Phase 1, document each criterion as present or absent with evidence before moving to Phase 2. The absence of most signals should block Phase 2.

2. **The SK006 suggestion text is the primary action, not a secondary fallback.** The validator message presents `references/` extraction first. Follow that ordering before considering skill splitting.

3. **Distinguish "sequential workflow phases" from "independently invocable skills."** A skill with 9 sequential phases where each depends on the previous is not a monolith with separable domains — it is a pipeline. Pipeline steps do not become separate skills; they become phases with their detailed reference content moved to `references/`.

---

## Summary

The orchestrator's primary error was skipping Phase 1's discriminating analysis and treating the refactor-skill's split workflow as the default response to an SK006 warning. The refactor-skill contributed to this by having no explicit exit path for the "not a split candidate" conclusion. Both the orchestrator and the skill need adjustment: the orchestrator needs to apply the domain criteria as a gate rather than preparation, and the skill needs to provide an explicit non-split path with `references/` extraction guidance.
