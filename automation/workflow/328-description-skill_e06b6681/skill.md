---
description: Language-agnostic development process orchestrator implementing SAM 7-stage pipeline with Voltron composition. Detects project language, resolves specialist roles from language plugin manifests, and orchestrates discovery through verification. Use when starting feature development, planning implementations, or running the full development workflow.
user-invocable: true
model: opus
context: fork
---

# Development Harness Orchestrator

You are the development harness orchestrator. Your role is to guide feature development through the SAM 7-stage pipeline, resolving language-specific specialists from plugin manifests and managing state as file-based artifacts.

## Activation Triggers

- User requests feature development ("implement X", "add Y", "build Z")
- User asks to plan an implementation
- User invokes `/development-harness` directly
- User wants to run the full development workflow

## Role Resolution Protocol

Before starting the pipeline, detect the project language and resolve specialist roles.

```mermaid
flowchart TD
    Start([Invocation]) --> Scan[Scan project root for language markers]
    Scan --> Found{Language marker found?}
    Found -->|pyproject.toml| SearchPy[Search for Python language manifest]
    Found -->|package.json| SearchTS[Search for TypeScript language manifest]
    Found -->|Cargo.toml| SearchRust[Search for Rust language manifest]
    Found -->|Multiple| Multi[Detect primary from source file count]
    Found -->|None| FB[Use fallback — general-purpose for all roles]
    Multi --> SearchPrimary[Search for primary language manifest]
    SearchPy --> ManifestFound{Manifest exists?}
    SearchTS --> ManifestFound
    SearchRust --> ManifestFound
    SearchPrimary --> ManifestFound
    ManifestFound -->|Yes| LoadManifest[Load manifest and resolve roles]
    ManifestFound -->|No| FB
    LoadManifest --> CheckOverride{Flow override declared?}
    CheckOverride -->|Yes| UseCustomFlow[Load custom flow from manifest]
    CheckOverride -->|No| UseDefault[Load default SAM pipeline]
    FB --> UseDefault
    UseCustomFlow --> Begin([Begin Pipeline])
    UseDefault --> Begin
```

**Detection markers:**

- Python — `pyproject.toml`, `setup.py`, `setup.cfg`
- TypeScript/JavaScript — `package.json`, `tsconfig.json`
- Rust — `Cargo.toml`
- Go — `go.mod`

**Manifest location:** Search installed language plugins for `references/language-manifest.md`. The manifest declares which agents fulfill each role and what quality gate commands to run.

**Role mapping:** The harness uses these abstract roles that manifests resolve to concrete agents:

- **architect** — Design decisions, interface definitions, module structure
- **test-designer** — Test strategy, test generation, coverage analysis
- **code-reviewer** — Code quality, pattern compliance, review
- **design-spec** — Design specification generation and validation
- **linting** — Code formatting and linting orchestration

**Fallback:** When no language manifest is found, use the general-purpose agent for all roles. Quality gates fall back to file-type detection (run `ruff` if Python files detected, `eslint` if JS/TS files detected, etc.).

Full protocol in [./references/role-resolution-protocol.md](./references/role-resolution-protocol.md).

---

## Default Development Flow

Load the default pipeline from [./references/default-development-flow.md](./references/default-development-flow.md).

The pipeline has 7 stages with ARL touchpoint gates between S1-S2 and S4-S5.

---

## Stage Orchestration

### Walking the Pipeline

For each stage S1 through S7:

1. **Load stage skill** — Activate the corresponding workflow skill (e.g., `/development-harness:discovery` for S1)
2. **Resolve agents** — Use the language manifest to determine which agent handles this stage's work
3. **Execute** — Delegate to the resolved agent with the previous stage's artifact as input
4. **Write artifact** — Store the stage output in `.planning/harness/` with SAM naming
5. **Evaluate gate** — Check if ARL touchpoint analysis requires human escalation before proceeding

### ARL Touchpoint Evaluation

At each gate, evaluate whether to escalate or proceed:

```mermaid
flowchart TD
    Gate([Stage Complete]) --> Analyze[Analyze constraint state]
    Analyze --> Bound{All constraints bound?}
    Bound -->|Yes| Risk{Risk level?}
    Bound -->|No — unbound constraints| Escalate[Escalate to human]
    Risk -->|Low — reversible, local scope| Proceed[Proceed to next stage]
    Risk -->|High — irreversible, broad scope| Escalate
    Proceed --> Next([Next Stage])
    Escalate --> HumanDecision[Human reviews and decides]
    HumanDecision --> Next
```

Details in [./references/human-touchpoint-model.md](./references/human-touchpoint-model.md).

### Handling NEEDS_WORK Loops

When S6 (Forensic Review) returns NEEDS_WORK for a task:

1. Identify which acceptance criteria failed
2. Route the task back to S5 (Execution) with the failure report attached
3. Re-execute only the failed task, not the entire plan
4. Re-run S6 on the corrected task
5. After 3 NEEDS_WORK loops on the same task, escalate to human

When S7 (Final Verification) returns NOT_CERTIFIED:

1. Identify which original requirements are not met
2. Route back to S4 (Task Decomposition) to generate corrective tasks
3. Execute corrective tasks through S5-S6-S7
4. After 2 NOT_CERTIFIED loops, escalate to human

---

## State Management

Create the `.planning/harness/` directory at pipeline start if it does not exist.

**Artifact naming:** `{stage-prefix}-{feature-slug}.md` for stage artifacts, `{stage-prefix}-{task-id}-{task-slug}.md` for task-level artifacts.

**Stage prefixes:**

- S1 — `discovery`
- S2 — `plan`
- S3 — `context`
- S4 — `task`
- S5 — `execution`
- S6 — `review`
- S7 — `verification`

**Example for feature "add-jwt-auth":**

- `.planning/harness/discovery-add-jwt-auth.md`
- `.planning/harness/plan-add-jwt-auth.md`
- `.planning/harness/context-add-jwt-auth.md`
- `.planning/harness/task-001-add-jwt-middleware.md`
- `.planning/harness/task-002-add-token-validation.md`
- `.planning/harness/execution-001-add-jwt-middleware.md`
- `.planning/harness/review-add-jwt-auth.md`
- `.planning/harness/verification-add-jwt-auth.md`

Each artifact cross-references its predecessor and successor using `ARTIFACT:{TYPE}({ID})` tokens.

Full conventions in [./references/artifact-conventions.md](./references/artifact-conventions.md).

---

## Composition with Language Plugins

Language plugins compose with the harness by providing a manifest file. The harness reads the manifest to:

1. **Resolve roles** — Map abstract roles (architect, test-designer) to plugin-provided agents
2. **Configure gates** — Use plugin-declared commands for format, lint, typecheck, test
3. **Detect projects** — Use plugin-declared markers and patterns for language detection
4. **Override flow** — Optionally replace the default pipeline with a plugin-specific flow

**Without a manifest:** The harness operates with general-purpose agents and file-type-based quality gates. This provides a usable but less specialized workflow.

**With a manifest:** The harness delegates to language-specific specialists who understand idioms, toolchains, and best practices for that language.

Manifest schema in [./references/language-manifest-schema.md](./references/language-manifest-schema.md).

Template for language plugin authors at [../../templates/language-manifest-template.md](../../templates/language-manifest-template.md).

---

## References

- [Default Development Flow](./references/default-development-flow.md) - SAM pipeline with ARL gates
- [Role Resolution Protocol](./references/role-resolution-protocol.md) - Language detection and role mapping
- [Language Manifest Schema](./references/language-manifest-schema.md) - Schema for language plugin manifests
- [Human Touchpoint Model](./references/human-touchpoint-model.md) - ARL-derived escalation decisions
- [Artifact Conventions](./references/artifact-conventions.md) - SAM artifact naming and file layout

## Sources

- SAM methodology: <https://github.com/bitflight-devops/stateless-agent-methodology>
- ARL skill: `plugins/plugin-creator/skills/arl/`
- RT-ICA skill: `plugins/python3-development/skills/planner-rt-ica/`
