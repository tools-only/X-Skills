---
description: Facilitated discussion of all 3 design phases, producing an Architecture Decision Record
---

# Review Design & Produce ADR

Facilitated discussion of all design phases from an `/arc` review. Produces an ADR as the final deliverable.

## Usage

`/arc-review [name]`

- `[name]` — The review name used in `/arc [name] [target]`

## Example

```
/arc-review payment-refactor
```

## Prerequisites

Run `/arc [name] [target]` first to generate:
- `docs/design-reviews/[name]/design.md`
- `docs/design-reviews/[name]/refined.md`
- `docs/design-reviews/[name]/critique.md`

## What Happens

Discuss each phase in order, then synthesize into an ADR:

### Phase 1: Initial Design
- Present summary of design.md
- Ask: "What did you think about the initial design?"
- Discuss strengths, weaknesses, concerns

### Phase 2: Refinements
- Present summary of refined.md and key changes
- Ask: "What did you think about the refinements?"
- Discuss whether refinements addressed the right issues

### Phase 3: Critique
- Present critique findings one by one
- For each: get user's decision (Accept/Reject/Defer/Partial)
- Record decisions and rationale

### Final: Produce ADR
- Synthesize all discussion into `docs/design-reviews/[name]/adr.md`

### Continuous Improvement
- Agent identifies where user made corrections
- Surfaces patterns for discussion
- Proposes process improvements

## Your Role

You drive the discussion. I facilitate and record:
- Share your reactions to each phase
- Identify strengths and concerns
- Make decisions on critique findings
- The ADR captures your final position

## Output

```
docs/design-reviews/[name]/
└── adr.md ← Architecture Decision Record (final deliverable)
```

---

# Instructions for Claude

You are facilitating a design review discussion across ALL 3 phases, then producing an ADR.

## Your Job

1. Read ALL documents first:
   - `docs/design-reviews/[name]/design.md`
   - `docs/design-reviews/[name]/refinements.md`
   - `docs/design-reviews/[name]/refined.md`
   - `docs/design-reviews/[name]/critique.md`
2. Discuss each phase IN ORDER with the user
3. WAIT for user response at each step
4. Track where user disagrees or makes corrections
5. After all phases discussed, write adr.md
6. Surface patterns and propose process improvements

## CRITICAL RULES

- **DO NOT** skip phases - discuss all 3 in order
- **DO NOT** rush through - let user share their thoughts
- **DO NOT** decide for the user
- **DO NOT** batch critique findings together
- **DO NOT** proceed without user input at each step

## Flow

### Phase 1: Initial Design

```
You: [Summarize key points from design.md]
You: "What did you think about the initial design? What are its strengths?"
User: [shares thoughts]
You: "Any concerns or things you felt were missing?"
User: [shares concerns]
```

### Phase 2: Refinements

```
You: [Summarize key changes from refinements.md]
You: "What did you think about the refinements? Did they address the right issues?"
User: [shares thoughts]
You: "Anything the refiner missed or over-corrected?"
User: [shares feedback]
```

### Phase 3: Critique Findings

```
You: "Now let's review the critique findings one by one."
You: "Finding 1 (CRITICAL): [description]. What's your decision?"
User: "Accept, because..."
You: [Record] "Finding 2 (HIGH): [description]. What's your decision?"
...continue until all findings addressed...
```

### Final: ADR

```
You: "All phases discussed. Writing the ADR to capture your final decisions."
[Write adr.md]
You: "ADR written to docs/design-reviews/[name]/adr.md"
```

### Continuous Improvement

After writing the ADR, analyze the review and surface patterns:

```
You: "Looking back at our discussion, here's where you had to make corrections or disagreed with the process:"

[List specific corrections, e.g.:]
- "The Architect missed X, which you had to add"
- "The Refiner over-complicated Y, which you simplified"
- "The Critique flagged Z as critical, but you rejected it because..."

You: "Let's discuss how we could build this expertise into the process. Some options:"
- "Add guidance to the Architect about X"
- "Constrain the Refiner to avoid Y pattern"
- "Teach the Critique about context Z"

You: "Which of these would be most valuable? Any other improvements?"
User: [discusses]
[Record agreed improvements in ADR]
```

## ADR Structure (adr.md)

```markdown
# ADR: [name]

**Status:** Proposed
**Date:** [date]

## Context

[Problem being solved, constraints, goals - from design.md + user's discussion]

## Decision

[Final recommended design - refined.md modified by accepted critique findings and user feedback]

## Consequences

### Positive
- [Benefits identified during discussion]

### Negative
- [Risks from accepted critique findings]
- [Trade-offs acknowledged]

### Mitigations
- [How negative consequences will be addressed]

## Alternatives Considered

### Original Design
[Summary from design.md]

**Why refined:** [Key changes and user's reasoning]

### Rejected Approaches
[Critique findings user REJECTED, with their rationale]

## Open Issues

[DEFERRED findings with timeline if provided]

## Process Improvements

Patterns identified during this review:

### Corrections Made
- [Specific correction and which phase it applied to]

### Proposed Improvements
- [Agreed improvement and which agent it applies to]
```

## ADR Content Guidelines

- **Context**: Combine design.md context + user's articulated understanding
- **Decision**: Start with refined.md, incorporate accepted findings + user's modifications
- **Alternatives**: Original design is one; rejected findings represent others
- **Process Improvements**: YOU identify patterns, propose options, user decides
- **Capture user voice**: The ADR should reflect their reasoning, not just your synthesis
