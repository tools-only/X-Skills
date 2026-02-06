---
name: arc-prd
description: "Architecture review for a PRD. Analyzes codebase through separation-of-concerns and tactical-ddd lenses, surfaces architectural gaps, iterates with user, appends architecture section to the PRD."
---

# PRD Architecture Review

Read the PRD at `$ARGUMENTS`. Verify its status is "Awaiting Architecture Review".

You are adding the architecture section that makes this PRD implementable. The PRD defines *what* to build — you define *how to structure it* at the package, feature, and domain model level. Not implementation details. Not function signatures. Structural decisions and domain modeling only.

## Step 0: Load Skills (MANDATORY — do this FIRST)

Before doing anything else, call both of these:
1. Skill(development-skills:separation-of-concerns)
2. Skill(development-skills:tactical-ddd)

Do NOT proceed to Phase 1 until both skills show "Successfully loaded skill".

## Phase 1: Research

1. Read the PRD
2. Read project conventions — search for: `docs/architecture/`, `ARCHITECTURE.md`, ADRs, `package.json` files, existing folder structure at the project root
3. Understand the current codebase structure through Glob and Grep
4. **Capture the current package graph** — run `npx nx show projects` and extract internal dependencies from each package's `package.json` (grep for `"@<org>/` entries). Record every package name and its dependency edges. This is the baseline for the Visual Overview.

## Phase 2: Analysis

Use both skills as **analytical lenses** — not just checklists to fill in, but tools to reveal what the PRD leaves architecturally ambiguous.

For each item in the architecture checklist below, do one of:
- **Propose a decision** — when the PRD + codebase provide enough signal. State the rationale grounded in a specific skill principle.
- **Surface a question** — when the skills reveal a gap, tension, or ambiguity the PRD didn't address. Include proposed options.

Present your analysis as:

```
## Proposed Decisions
[Concrete proposals with rationale citing specific skill principles]

## Open Questions
[Gaps/tensions discovered through skill analysis, each with proposed options]
```

### How to surface good questions

Questions should be grounded in specific skill principles and reference the PRD concretely:

- "PRD describes X as part of feature Y, but separation-of-concerns principle 2 says this is shared if multiple features use it. Does Z also use this? → determines features/ vs platform/domain/"
- "This introduces [concept] but doesn't define invariants. Tactical-DDD principle 7: what must always be true when [state change] happens?"
- "PRD mentions [interface] but doesn't say whether it's a domain concept or infrastructure. This affects where it lives."
- "The PRD has [component A] calling [component B] — is B a use case (menu test) or internal machinery?"
- "PRD introduces [term] without defining it. Is this a new domain concept? An existing one being extended? A generic utility?"

Don't ask questions you can answer from the codebase. Research first.

## Phase 3: Iterate

Discuss with user until all questions are resolved and decisions are agreed. Show, don't tell — use ASCII diagrams, folder structure sketches, before/after comparisons.

## Phase 4: Write

Append the `## Architecture` section to the PRD file (replacing the placeholder). Update the PRD status to "Approved".

### Visual Overview (MANDATORY — write this FIRST)

The architecture section MUST begin with a Visual Overview containing three mermaid diagrams. These diagrams give the reader a 30-second visual summary before detailed decisions.

Use these exact colour definitions across all diagrams:

```
classDef existing fill:#e8e8e8,stroke:#999,color:#333
classDef modified fill:#fff3e0,stroke:#f57c00,color:#333,stroke-width:3px
classDef new fill:#e8f5e9,stroke:#2e7d32,color:#1b5e20,stroke-width:3px
```

Additional colours for diagram 2 only:

```
classDef feature fill:#e8f5e9,stroke:#2e7d32,color:#1b5e20,stroke-width:3px
classDef horizontal fill:#e3f2fd,stroke:#1565c0,color:#0d47a1,stroke-width:3px
classDef contract fill:#fff3e0,stroke:#f57c00,color:#333,stroke-width:3px
```

#### Diagram 1: Package Map

Shows ALL packages in the project and their dependency edges. Built from the nx data captured in Phase 1.

- `graph TD` flowchart
- Every package is a node, classified as `existing`, `modified`, or `new`
- Existing dependency edges use `-->` with `linkStyle default stroke:#999,stroke-width:1px`
- New dependency edges use `==>` with green `linkStyle` overrides: `stroke:#2e7d32,stroke-width:3px`
- New edges are labeled with what crosses the boundary (e.g. `|DraftComponent|`)
- Followed by a separate legend diagram: `graph LR` with `a[Existing]:::existing ~~~ b[Modified]:::modified ~~~ c[New]:::new`

#### Diagram 2: New Feature Detail

🚨 **ONE diagram per new feature.** Do NOT combine multiple features into a single diagram. Each feature gets its own mermaid diagram with its own markdown heading `### New Feature: <name>`.

- `graph TD` flowchart
- Subgraphs represent packages — colour the subgraph border to match its status:
  - New package: `style <id> fill:none,stroke:#2e7d32,stroke-width:3px,color:#1b5e20`
  - Modified package: `style <id> fill:#fff3e0,stroke:#f57c00,stroke-width:2px,color:#333`
  - Existing package: `style <id> fill:#e8e8e8,stroke:#999,stroke-width:2px,color:#333`
- Inside each package subgraph, show what's new/changed:
  - New features (verticals) use `feature` class (green)
  - New shared capabilities (horizontals) use `horizontal` class (blue)
  - Changed types/interfaces use `contract` class (orange)
- Show the consumer (what calls this feature) as an existing node at the top
- Arrow styles: thick green for new usage, thick blue for internal horizontal usage, dotted orange for contract dependencies
- Followed by a separate legend diagram with: `a[New Feature]:::feature ~~~ b[New Shared Capability]:::horizontal ~~~ c[Changed Type/Interface]:::contract`

#### Diagram 3: Domain Model

🚨 **ONE diagram per connected group of domain concepts.** Only group concepts in the same diagram if they have dependency relationships (edges between them). If the PRD introduces domain changes that are unrelated to each other, they MUST be separate diagrams. Each diagram gets its own markdown heading `### Domain Model: <name>`.

**Splitting rule:** After drafting the domain model, check if the diagram contains multiple disconnected subgraphs (groups of nodes with no edges between them). If yes, split into separate diagrams — each with its own title. A single diagram must be a connected graph.

- `graph TD` flowchart (NOT erDiagram — classDef not supported in ER diagrams by most renderers)
- Subgraphs represent packages — same border colour rules as diagram 2
- Nodes are domain concepts (aggregates, entities, value objects), classified as `existing`, `modified`, or `new`
- Edges use semantic relationship labels (e.g. `"traced by"`, `"produces"`, `"configures"`) — NOT database cardinality
- Followed by a separate legend diagram with: `a[Existing]:::existing ~~~ b[Modified]:::modified ~~~ c[New]:::new`

#### External Dependencies

After the three diagrams, add a markdown table:

```
### External Dependencies

| Dependency | Package | Purpose | Status |
|-----------|---------|---------|--------|
```

Status is either `existing` or `NEW`. If no external dependencies, write "None."

### Notation Rules

🚨 **Every visual element must use colour/shape to convey ONE thing only.** Do not encode multiple concepts into the same visual element (e.g. don't put package name and status in the same text label — use node colour for status and text for the name).

🚨 **Containers (subgraphs) represent WHERE code lives.** A reader must be able to look at any node and immediately know which package it belongs to.

🚨 **Colours are consistent across all diagrams.** Grey = existing/unchanged. Orange = modified. Green = new. Blue = shared horizontal capability. Never reuse a colour for a different meaning.

### Annotate deliverables

After writing the architecture section, go back through the PRD's milestone deliverables and add architecture references. Under each deliverable that is affected by an architecture decision, add:

```
Architecture: see §9.X ([decision summary])
```

This creates an explicit link from each deliverable to the architecture decisions that constrain it. Task creation will use these references to inject architecture context into each task.

Examples:
- Under a deliverable about a new extraction engine: `Architecture: see §9.1 (new riviere-connection-detection package), §9.3 (call graph as platform/domain)`
- Under a deliverable introducing an aggregate: `Architecture: see §9.5 (Builder aggregate invariants), §9.6 (upsert merge semantics)`

Every deliverable should have at least one architecture reference. If a deliverable has none, either the architecture section is incomplete or the deliverable is purely documentation.

## Architecture Checklist

The architecture section must answer all of these:

### Structural decisions (from separation-of-concerns)

- Which packages are modified vs created?
- What new features (verticals) are introduced?
- What shared capabilities (horizontals) are needed?
- For each feature: which layers apply? (entrypoint / commands / queries / domain)
- What goes in platform/domain vs platform/infra?
- What external clients are introduced or modified?
- Dependency direction between packages (shown in Package Map diagram)

### Domain model (from tactical-ddd)

- What aggregates/entities are introduced or changed?
- What are the key invariants?
- What use cases exist? (apply the menu test)
- What value objects emerge?
- What state transitions matter?
- What domain language/terminology is introduced?

### Integration

- Key interfaces/contracts between new and existing code

### Flexibility markers

Mark each decision as:
- **Firm** — structural, hard to change later, get it right now
- **Flexible** — can iterate during implementation, direction is set but details may shift

## Rules

1. **Do not define implementation details** — no function signatures, no class hierarchies, no algorithms. Structure and domain model only.
2. **Do not duplicate the PRD** — reference it, don't repeat it.
3. **Do not prescribe the full solution** — leave room for implementation to iterate within the structural boundaries you define.
4. **Ground every decision in a skill principle** — if you can't cite why, it's opinion not architecture.
5. **Show, don't tell** — diagrams, folder trees, and concrete examples over prose.
