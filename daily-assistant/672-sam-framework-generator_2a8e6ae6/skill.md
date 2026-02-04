# SAM Framework Generator

**Status**: Conceptual Design
**Purpose**: Generate project-specific SAM implementations through a meta-SAM workflow
**Date**: 2026-01-26

---

## Core Concept

Use **SAM to create SAM** - apply the Stateless Agent Methodology to generate project-specific SAM workflows.

Instead of manually creating commands, agents, and skills for each project type, the SAM Framework Generator uses the SAM 7-stage pipeline to discover requirements, research patterns, design the workflow, and generate all necessary artifacts.

---

## The Meta-SAM Pattern

```
┌────────────────────────────────────────────────────────────────┐
│              SAM FRAMEWORK GENERATOR (META-SAM)                 │
├────────────────────────────────────────────────────────────────┤
│                                                                 │
│  Input: Project type intent (e.g., "Python CLI tool")          │
│  Output: Complete SAM implementation for that project type      │
│                                                                 │
│  Stage 1: DISCOVERY                                             │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │ Interview user about:                                     │  │
│  │ - What workflows does this project type need?            │  │
│  │ - What are the common operations? (add feature, fix bug) │  │
│  │ - What technologies/patterns are used?                   │  │
│  │ - What constraints exist? (Python version, frameworks)   │  │
│  │ - What verification is required? (tests, linting)        │  │
│  │ - What artifacts should be generated?                    │  │
│  │                                                           │  │
│  │ Output: Project Type Requirements Document                │  │
│  └──────────────────────────────────────────────────────────┘  │
│                              ▼                                  │
│  Stage 2: PLANNING (RT-ICA)                                     │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │ Verify prerequisites:                                     │  │
│  │ - Are there existing SAM templates?                      │  │
│  │ - Are the workflow patterns well-defined?                │  │
│  │ - Do we have reference implementations?                  │  │
│  │                                                           │  │
│  │ BLOCKS if missing SAM component templates                │  │
│  │ Output: SAM Generation Plan                               │  │
│  └──────────────────────────────────────────────────────────┘  │
│                              ▼                                  │
│  Stage 3: CONTEXT INTEGRATION                                   │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │ Research project type patterns:                           │  │
│  │ - Analyze existing similar projects                      │  │
│  │ - Identify common file structures                        │  │
│  │ - Map technology stack patterns                          │  │
│  │ - Find existing SAM implementations to reference         │  │
│  │                                                           │  │
│  │ Output: Contextualized SAM Design                         │  │
│  └──────────────────────────────────────────────────────────┘  │
│                              ▼                                  │
│  Stage 4: TASK DECOMPOSITION                                    │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │ For each workflow identified in discovery:               │  │
│  │ - Design the 7-stage SAM pipeline for that workflow      │  │
│  │ - Identify specialized agents needed                     │  │
│  │ - Define artifact templates                              │  │
│  │ - Specify verification gates                             │  │
│  │                                                           │  │
│  │ Output: SAM Component Generation Tasks                    │  │
│  └──────────────────────────────────────────────────────────┘  │
│                              ▼                                  │
│  Stage 5: EXECUTION                                             │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │ Generate artifacts:                                       │  │
│  │ - .claude/commands/{workflow}.md (orchestrator)          │  │
│  │ - .claude/agents/{role}.yaml (specialized agents)        │  │
│  │ - .claude/sam-templates/{artifact}.md (templates)        │  │
│  │ - .claude/sam-config.yaml (project configuration)        │  │
│  │                                                           │  │
│  │ Output: Generated SAM Implementation                      │  │
│  └──────────────────────────────────────────────────────────┘  │
│                              ▼                                  │
│  Stage 6: FORENSIC REVIEW                                       │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │ Verify generated SAM implementation:                      │  │
│  │ - Commands follow orchestrator discipline                │  │
│  │ - Agents have proper DONE/BLOCKED signaling              │  │
│  │ - Artifact templates are complete                        │  │
│  │ - Workflow stages match SAM principles                   │  │
│  │                                                           │  │
│  │ Output: Validation Report                                 │  │
│  └──────────────────────────────────────────────────────────┘  │
│                              ▼                                  │
│  Stage 7: FINAL VERIFICATION                                    │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │ Confirm generated SAM meets requirements:                 │  │
│  │ - All workflows from discovery have implementations      │  │
│  │ - Agents cover all identified roles                      │  │
│  │ - Templates support required artifacts                   │  │
│  │                                                           │  │
│  │ Output: Project-Specific SAM Ready for Use               │  │
│  └──────────────────────────────────────────────────────────┘  │
│                                                                 │
└────────────────────────────────────────────────────────────────┘
```

---

## Discovery Phase Interview Structure

The Discovery Agent asks structured questions to gather complete requirements:

### Workflow Identification

**Questions:**

- What are the primary workflows for this project type? (e.g., "add feature", "fix bug", "refactor module")
- For each workflow, what triggers it? (user request, automated process)
- What is the typical complexity? (simple, moderate, complex)
- What are the success criteria for each workflow?

### Technology Stack

**Questions:**

- What programming language(s)?
- What frameworks or libraries? (e.g., Typer, Django, FastAPI, Rust clap)
- What build tools? (e.g., uv, cargo, npm)
- What testing frameworks? (pytest, jest, cargo test)
- What linting/formatting tools? (ruff, mypy, clippy, prettier)

### Domain Constraints

**Questions:**

- What architectural patterns are used? (CLI, web app, microservices)
- What are the code organization conventions? (modules, packages, directory structure)
- What data models or schemas are common?
- What external systems are integrated? (APIs, databases, SSH)
- What security considerations apply?

### Artifact Requirements

**Questions:**

- What documents should be generated during planning? (architecture specs, API docs)
- What format should task files use?
- What metadata should be tracked?
- What reports are needed for verification?

### Agent Specializations

**Questions:**

- What roles are needed? (researcher, architect, implementer, tester, reviewer)
- What domain expertise should each agent have?
- What tools should each agent have access to?
- What are the DONE/BLOCKED criteria for each agent?

### Verification Gates

**Questions:**

- What automated checks should run? (tests, linting, security scans)
- What manual review points are needed?
- What are the quality thresholds? (coverage %, linting pass rate)
- What are the acceptance criteria for complete features?

---

## Generated Artifacts

For a project type (e.g., "Python CLI Tool"), the generator produces:

### 1. Workflow Commands

`.claude/commands/add-new-feature.md`

- Orchestrates 7-stage SAM pipeline for adding features
- Delegates to specialized agents
- Uses TodoWrite for progress tracking
- Generates feature-specific artifacts

`.claude/commands/fix-bug.md`

- Orchestrates SAM pipeline for bug fixes
- Includes root cause analysis stage
- Generates bug reports and patches

`.claude/commands/refactor-module.md`

- Orchestrates SAM pipeline for refactoring
- Includes technical debt assessment
- Generates refactoring plans and verification reports

### 2. Specialized Agents

`.claude/agents/python-cli-researcher.yaml`

- Role: Discovery and requirements gathering
- Tools: Read, Grep, Glob, AskUserQuestion
- Expertise: Python CLI patterns, Typer conventions
- Output: Feature requirements document

`.claude/agents/python-cli-architect.yaml`

- Role: Architecture design
- Tools: Read, Write, Grep, Glob
- Expertise: Pydantic models, CLI design, Python 3.11+
- Output: Architecture specification

`.claude/agents/python-task-planner.yaml`

- Role: Task decomposition
- Tools: Read, Write, Grep, Glob
- Expertise: TDD patterns, dependency analysis
- Output: Task file with dependencies

`.claude/agents/python-cli-implementer.yaml`

- Role: Code implementation
- Tools: Read, Write, Edit, Bash
- Expertise: Python 3.11+, Typer, Pydantic, async patterns
- Output: Implementation + self-verification results

`.claude/agents/python-test-architect.yaml`

- Role: Test creation
- Tools: Read, Write, Edit, Bash
- Expertise: pytest, pytest-mock, hypothesis, coverage
- Output: Test suite + coverage report

`.claude/agents/python-forensic-reviewer.yaml`

- Role: Independent verification
- Tools: Read, Bash, Grep
- Expertise: Code review, quality assessment
- Output: Review report with COMPLETE/NEEDS_WORK verdict

### 3. Artifact Templates

`.claude/sam-templates/feature-requirements.md`

```markdown
# Feature Requirements: {feature-name}

## Who
{target users}

## When
{trigger conditions}

## Where
{integration point}

## Inputs
{parameters, files, env vars}

## Outputs
{results, exit codes, side effects}

## Guardrails
{validation, safety checks}

## Conditionals
{edge cases, error handling}

## Naming
{CLI command/option names}

## Acceptance Criteria
- [ ] Criterion 1
- [ ] Criterion 2
- [ ] Criterion 3
```

`.claude/sam-templates/architecture-spec.md`

```markdown
# Architecture Specification: {feature-name}

## Overview
{high-level design}

## CLI Interface
{Typer command structure}

## Component Changes
{modules affected}

## Data Models
{Pydantic models}

## Error Handling
{exception strategy}

## Security Considerations
{security analysis}

## Testing Strategy
{test approach}
```

`.claude/sam-templates/task-file.md`

```markdown
# Implementation Tasks: {feature-name}

## Task 1: {name}
**Status**: pending
**Dependencies**: []
**Priority**: high
**Complexity**: moderate
**Agent**: python-cli-implementer

**Acceptance Criteria**:
- [ ] Criterion 1
- [ ] Criterion 2

**Verification Steps**:
- [ ] Step 1
- [ ] Step 2

**Context Manifest**: {will be added by context-gathering agent}
```

### 4. Project Configuration

`.claude/sam-config.yaml`

```yaml
# SAM Configuration for Python CLI Tool Projects

project_type: python-cli-tool
framework_version: "1.0.0"

# Technology Stack
languages:
  - python: "3.11+"

frameworks:
  - typer: "CLI framework"
  - pydantic: "v2 - Data validation"
  - rich: "Terminal UI"

build_tools:
  - uv: "Package management and script runner"

testing:
  - pytest: "Test framework"
  - pytest-mock: "Mocking"
  - pytest-asyncio: "Async tests"
  - hypothesis: "Property-based testing"
  - coverage: "Code coverage (80% minimum)"

quality:
  - ruff: "Linting and formatting"
  - mypy: "Type checking (strict mode)"

# SAM Workflow Configuration
workflows:
  add-feature:
    command: /add-new-feature
    stages: 7
    agents:
      - python-cli-researcher
      - python-cli-architect
      - python-task-planner
      - python-cli-implementer
      - python-test-architect
      - python-forensic-reviewer
    artifacts:
      - feature-requirements.md
      - architecture-spec.md
      - task-file.md
      - implementation-results.md
      - review-report.md

  fix-bug:
    command: /fix-bug
    stages: 7
    agents:
      - bug-investigator
      - root-cause-analyzer
      - python-cli-implementer
      - python-test-architect
      - python-forensic-reviewer
    artifacts:
      - bug-report.md
      - root-cause-analysis.md
      - fix-plan.md
      - implementation-results.md
      - review-report.md

# Verification Gates
verification:
  tests:
    required: true
    min_coverage: 80
    frameworks: [pytest]

  linting:
    required: true
    tools: [ruff, mypy]
    strict_mode: true

  security:
    required: false
    tools: []

# Domain Constraints
constraints:
  python_version: "3.11+"
  type_hints: required
  async_patterns: preferred
  cli_framework: typer
  data_validation: pydantic-v2
```

---

## SAM Component Templates

These templates are used by the generator to create project-specific agents and commands.

### Agent Template

`.claude/sam-templates/agent-template.yaml`

```yaml
# SAM Agent Template
# Variables: {ROLE_NAME}, {EXPERTISE}, {TOOLS}, {OUTPUT_ARTIFACT}

name: {ROLE_NAME}
description: |
  SAM {STAGE_NAME} stage agent for {PROJECT_TYPE}.
  {ROLE_DESCRIPTION}

model: claude-sonnet-4.5

tools:
  {TOOLS_LIST}

prompt: |
  You are a {ROLE_NAME} agent in a Stateless Agent Methodology (SAM) workflow.

  Your ROLE_TYPE is sub-agent.

  ## Your Responsibilities

  {RESPONSIBILITIES}

  ## Domain Expertise

  {EXPERTISE}

  ## Input Artifacts

  {INPUT_ARTIFACTS}

  ## Output Requirements

  You MUST produce: {OUTPUT_ARTIFACT}

  The artifact MUST include:
  {OUTPUT_REQUIREMENTS}

  ## Completion Signaling

  Return STATUS: DONE when:
  - {DONE_CRITERIA}

  Return STATUS: BLOCKED when:
  - {BLOCKED_CRITERIA}

  ## Verification Protocol

  Before signaling DONE:
  1. {VERIFICATION_STEP_1}
  2. {VERIFICATION_STEP_2}
  3. {VERIFICATION_STEP_3}
```

### Command Template

`.claude/sam-templates/command-template.md`

```markdown
---
description: {WORKFLOW_DESCRIPTION}
argument-hint: {ARGUMENT_HINT}
---

# {WORKFLOW_NAME}

Execute the {WORKFLOW_NAME} workflow using the Stateless Agent Methodology (SAM).

<workflow_request>
$ARGUMENTS
</workflow_request>

---

## Orchestrator Discipline

**CRITICAL**: You are an orchestrator. You coordinate work across specialized agents. You do NOT perform discovery work directly.

- **NEVER** use Read tool to gather information for decision-making
- **NEVER** use Glob tool to search the codebase
- **NEVER** use Grep tool to find patterns
- **ALWAYS** delegate discovery and analysis to specialized agents
- **ALWAYS** use TodoWrite for progress tracking (your state management tool)
- **CAN** read files ONLY for state management (checking task status, not for discovery)

Your role: Route context between user and agents. Define success criteria. Track progress. Trust agent expertise.

---

## Mission

Orchestrate the execution of all {STAGE_COUNT} SAM stages sequentially. Each stage depends on the previous stage's output. After each stage completes, RECEIVE the agent's structured STATUS response before proceeding to the next stage.

---

## Completion Tracking (MANDATORY)

**IMMEDIATELY** after reading this command, you MUST create todos using TodoWrite with this exact checklist:

```

TodoWrite(todos=[
{STAGE_TODOS}
])

```

**RULES**:

1. Create ALL todos BEFORE starting Stage 1
2. Mark each todo `in_progress` BEFORE launching the agent
3. Mark each todo `completed` AFTER agent returns STATUS: DONE
4. If agent returns STATUS: BLOCKED, keep todo as `in_progress` and address the blocker
5. DO NOT display final summary until ALL todos are `completed`

---

{STAGE_DEFINITIONS}

---

## Final Summary

**BEFORE displaying the final summary, YOU MUST:**

1. **VERIFY ALL TODOS COMPLETE**: All todos from Completion Tracking should be `completed`
2. If ANY todo is still `pending` or `in_progress`, DO NOT display the final summary - complete the missing work first
3. Mark "Final: Display completion summary" as `in_progress`
4. DISPLAY this summary structure:

```

================================================================================
{WORKFLOW_NAME} COMPLETE
================================================================================

Request: $ARGUMENTS

## DELIVERABLES

{DELIVERABLES_LIST}

## PHASE SUMMARIES

{PHASE_SUMMARIES}

## NEXT STEPS (for user)

{NEXT_STEPS}

## FILES TO REVIEW

{FILES_TO_REVIEW}

```

5. Mark "Final: Display completion summary" as `completed`

**WORKFLOW COMPLETE**
```

### Stage Template

`.claude/sam-templates/stage-template.md`

```markdown
## Stage {N}: {STAGE_NAME}

**Agent**: `{AGENT_NAME}`
**Purpose**: {STAGE_PURPOSE}

### Delegation

```

Task(
subagent_type="{AGENT_NAME}",
description="{STAGE_DESCRIPTION}",
prompt="""
Your ROLE_TYPE is sub-agent.

OBSERVATIONS:
{OBSERVATIONS}

DEFINITION OF SUCCESS:
{SUCCESS_CRITERIA}

CONTEXT:
{CONTEXT}

YOUR TASK:
{TASK_STEPS}

AVAILABLE RESOURCES:
{RESOURCES}
"""
)

```

### Handle Result

IF STATUS: DONE → {DONE_ACTION}
IF STATUS: BLOCKED → {BLOCKED_ACTION}
```

---

## Implementation Phases

### Phase 1: Core Infrastructure

**Deliverables:**

- SAM generator command: `/sam:generate`
- Discovery agent for requirements gathering
- Project type registry (initial types: python-cli, django-web, rust-binary)
- Template storage in `.claude/sam-templates/`

### Phase 2: Template System

**Deliverables:**

- Agent template system
- Command template system
- Artifact template system
- Variable substitution engine

### Phase 3: Project Type Library

**Deliverables:**

- Python CLI tool templates
- Django web app templates
- Rust binary templates
- React/Next.js frontend templates

### Phase 4: Customization

**Deliverables:**

- Project-specific overrides
- Custom agent injection
- Workflow composition
- Template inheritance

### Phase 5: Validation & Testing

**Deliverables:**

- Generated artifact validation
- SAM principle compliance checks
- Integration testing harness
- Example generated projects

---

## Example Usage

### Generate SAM for Python CLI Project

```bash
# User invokes generator
/sam:generate python-cli-tool

# Discovery agent interviews user
Agent: "What workflows does this CLI tool need?"
User: "Add commands, fix bugs, refactor modules"

Agent: "What CLI framework do you use?"
User: "Typer with Rich for UI"

Agent: "What testing approach?"
User: "pytest with 80% coverage minimum"

# ... more questions ...

# Generator produces artifacts
✓ Generated: .claude/commands/add-new-feature.md
✓ Generated: .claude/commands/fix-bug.md
✓ Generated: .claude/commands/refactor-module.md
✓ Generated: .claude/agents/python-cli-researcher.yaml
✓ Generated: .claude/agents/python-cli-architect.yaml
✓ Generated: .claude/agents/python-task-planner.yaml
✓ Generated: .claude/sam-templates/feature-requirements.md
✓ Generated: .claude/sam-config.yaml

SAM framework ready for: Python CLI Tool
Available workflows: /add-new-feature, /fix-bug, /refactor-module
```

### Use Generated Workflow

```bash
# User invokes generated workflow
/add-new-feature Add a command that lists all out-of-date packages on remote managed hosts over SSH, with search and filtering, and options for selecting what to install or installing all updates. With a toggle for a post-update reboot, or to schedule a reboot for its maintenance window.

# Orchestrator follows SAM 7-stage pipeline
Stage 1: Discovery → python-cli-researcher agent
Stage 2: Planning → validates prerequisites
Stage 3: Context Integration → maps to existing code
Stage 4: Task Decomposition → generates task file
Stage 5: Execution → python-cli-implementer agent
Stage 6: Forensic Review → python-forensic-reviewer agent
Stage 7: Final Verification → confirms feature complete

# Complete with all artifacts generated
✓ Feature requirements documented
✓ Architecture spec created
✓ Task file with dependencies
✓ Implementation verified
✓ Tests passing at 85% coverage
```

---

## Integration with Existing SAM Concepts

### Convergence with OctoCode RDD

The generator can produce workflows that use OctoCode's adversarial validation:

```yaml
# In .claude/sam-config.yaml
verification_style: adversarial

# Generates agents with Generator/Discriminator pattern
agents:
  - python-cli-generator  # Produces code
  - python-cli-verifier   # Finds flaws (adversarial)
```

### Convergence with Get Shit Done (GSD)

The generator can produce workflows that use GSD's wave execution:

```yaml
# In .claude/sam-config.yaml
execution_style: wave-based

# Generates task files with parallelization metadata
task_format: gsd-compatible
```

---

## Naming Considerations

Following the pop culture naming exploration, project types could have themed names:

| Project Type   | SAM Nickname           | Rationale                                |
| -------------- | ---------------------- | ---------------------------------------- |
| Python CLI     | "The Ford Line"        | Assembly line metaphor for CLI pipelines |
| Django Web     | "The Groundhog Loop"   | Iteration cycles for web development     |
| Rust Binary    | "The Apollo Pattern"   | Systems rigor and verification           |
| React Frontend | "The Portal Framework" | Component portals and state management   |
| Data Pipeline  | "The Deming Cycle"     | PDCA quality methodology                 |
| Microservices  | "The Memento Pattern"  | Stateless services with external memory  |

These names appear in:

- Generated command descriptions
- Agent role definitions
- Artifact templates
- Project documentation

---

## Success Metrics

| Metric                         | Target | Measurement                                     |
| ------------------------------ | ------ | ----------------------------------------------- |
| **Generation time**            | <5 min | Time from `/sam:generate` to complete artifacts |
| **Workflow completeness**      | 100%   | All identified workflows have implementations   |
| **Agent specialization**       | 5-7    | Agents per workflow for proper separation       |
| **Artifact template coverage** | 100%   | Templates for all required documents            |
| **SAM principle compliance**   | 100%   | Generated artifacts follow SAM structure        |
| **User customization**         | <10min | Time to customize generated workflow            |

---

## Future Enhancements

### Community Template Registry

Allow sharing of project type templates:

```bash
/sam:install community/nestjs-backend
/sam:install community/svelte-frontend
```

### Cross-Project Learning

Generator learns from successful implementations:

```bash
# After successful feature implementation
/sam:learn-pattern add-new-feature

# Incorporates pattern into future generations
```

### Template Composition

Combine multiple project types:

```bash
/sam:generate python-cli + fastapi-backend
# Generates workflows for both CLI and API development
```

---

## Conclusion

The SAM Framework Generator is a **meta-SAM** that uses the Stateless Agent Methodology to create project-specific SAM implementations.

**Core Insight**: Instead of manually creating SAM workflows for each project type, use SAM's own principles to discover, design, and generate those workflows.

**Key Benefits**:

1. **Consistency**: All generated workflows follow SAM principles
2. **Customization**: Tailored to project-specific technologies and constraints
3. **Evolution**: Easy to update templates and regenerate
4. **Reusability**: Templates can be shared and versioned
5. **Learning**: Captures best practices in template form

**The Recursion**: SAM generates SAM, which generates working software through SAM.

**Next Steps**:

1. Implement Phase 1 (Core Infrastructure)
2. Create initial templates for Python CLI tools
3. Validate with real-world project generation
4. Expand to additional project types
5. Build community template registry
