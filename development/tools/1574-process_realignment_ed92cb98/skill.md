## Local glossary

- **Desired outcome**: the end state/value (what "done" means at the user/business level).
- **Objectives**: component outcomes that collectively achieve the desired outcome (often phrased as "the system/user can...").
- **Acceptance criteria**: testable statements proving each objective is met (observable, verifiable, pass/fail).
- **Skill**: a structured, in-band process wrapper (system-prompt-like structure) that drives the conversation and dispatches work.
- **Skill(fork)**: a forked subprocess skill (isolated context) that performs work and returns a structured summary + artifacts back to the main conversation.

Draft Intended workflow:

1. Skill: Discovery of desired outcome and objectives (user provided a feature, issue, task, enablement, structure change, documentation, pipeline and release flow, versioning, stakeholders, ... or anything else that comes in that is part of the project could be provided as the input and this should be treated in a flexible manner)

   NOTE: By "goal" in earlier conversation, I mean the combination of:

   - Desired outcome (end state/value)
   - Objectives (component outcomes)
   - Acceptance criteria (verification)

   A common structure is:

   - Desired outcome → 1..N objectives
   - Each objective → 1..N acceptance criteria (plus optional "anti-criteria" / non-goals)
     a. Primary Assistant: Context Gathering
     - conversational
     - dispatches researcher agents
       - Skill: Agent with tasks for online research
       - Skill: Agents with tasks for local repository discovery
       - Skill: Agents with git analysis tasks
       - Skill: Agents with GitHub or GitLab remote (or local filepath) repository reference gathering tasks.
       - Skill: Agents that check packages, latest version documentation for any mentioned technologies, i.e. if the user says it for gcc 10.3 then ensure the reference doc or url for gcc 10.3 has been collected as reference material. If the user says we are creating a release pipeline for GitLab, we should know what the current gitlab pipeline schema is and the available predefined variables.
     - Skill: validates findings with user asking contextual questions about if the findings and resources to ensure alignment towards intended goal
       b. review and restate the goal and intent to the user until user approves

2. Skill: Architecture/Design - Identifies what is required to achieve the desired outcome + objectives

   NOTE: Artifacts selected based on:

   - Project complexity (scope, cross-system impact, infrastructure changes)
   - Architectural role ownership (Systems/Solutions/Software perspectives)
   - Explicit user request

   COMPLEXITY ASSESSMENT:

   - **Simple**: Single-module, no external interfaces → ARCH-COMPONENT, API only
   - **Moderate**: Multi-module, internal APIs, config → add ARCH-CONTAINER, ARCH-SOFTWARE, SAD-SOLUTION
   - **Complex**: Cross-system, infrastructure, hardware → add ARCH-CONTEXT, ICD, DEPLOYMENT
   - **Firmware/Embedded**: Hardware interaction → add HAL

   --- Requirements & Scope (ALWAYS) ---
   a. Skill(fork): Project Requirements Documentation

   - Desired outcome, objectives, acceptance criteria
   - Feature scope, boundaries, stakeholders
     `ARTIFACT:PRD(SCOPE:...)`

   b. Skill(fork): Non-Functional Requirements

   - Performance, scalability, reliability, security constraints
   - Technology constraints, platform requirements
     `ARTIFACT:NFR(SCOPE:...)`

   --- Architecture Views (complexity-dependent) ---
   c. Skill(fork): Architecture Artifacts (C4 Model + arc42)
   Decides threshold for each level based on project complexity.
   Input: Discovery summary + PRD + NFR

   i. System Context Diagram (C4 L1 - Systems Architect)
   WHEN: External interfaces, integrations, system boundaries - System boundary, external actors, data flows in/out
   ROLE: Systems Architect perspective
   SOURCE: C4 Model <https://c4model.com/diagrams/system-context>
   `ARTIFACT:ARCH-CONTEXT(SCOPE:...)`

   ii. Container Diagram (C4 L2 - Solutions Architect)
   WHEN: Multi-service, technology selections, application boundaries - Applications, databases, queues, communication protocols - Technology choices per container
   ROLE: Solutions Architect perspective
   SOURCE: C4 Model <https://c4model.com/diagrams/container>
   `ARTIFACT:ARCH-CONTAINER(SCOPE:...)`

   iii. Component Diagram (C4 L3 - Software Architect)
   WHEN: Code changes (DEFAULT for most development work) - Internal components, interfaces, dependencies - Module organization and boundaries
   ROLE: Software/Firmware Architect perspective
   SOURCE: C4 Model <https://c4model.com/diagrams/component>
   `ARTIFACT:ARCH-COMPONENT(SCOPE:...)`

   --- Meta Project Infrastructure (when infra changes) ---
   d. Skill(fork): Meta Project Infrastructure Documentation
   WHEN: Infrastructure/tooling/process changes

   - Linting, packaging, CI/CD, tooling, contributing, debugging
     `ARTIFACT:MPI(CATEGORY:...)`

   --- Detailed Specifications (complexity-dependent) ---
   e. Skill(fork): Detailed Architecture Artifacts

   i. Interface Control Document (Systems Architect)
   WHEN: External system integration, hardware-software boundaries - Protocols, message formats, external API contracts - Hardware interface specifications
   ROLE: Systems Architect perspective
   `ARTIFACT:ICD(INTERFACE:...)`

   ii. API Specification (Software Architect)
   WHEN: New public interfaces, REST/RPC/GraphQL endpoints - Endpoints, request/response, data types, errors - Authentication, rate limits, versioning
   ROLE: Software Architect perspective
   FORMAT: OpenAPI, protobuf, GraphQL schema, etc.
   `ARTIFACT:API-SPEC(MODULE:...)`

   iii. Software Architecture Document (Software Architect)
   WHEN: New features, refactoring, architectural changes - Layers, design patterns, module organization - Internal component relationships
   ROLE: Software Architect perspective
   SOURCE: arc42 Section 5 - Building Block View
   `ARTIFACT:ARCH-SOFTWARE(SCOPE:...)`

   iv. Module Design Document (Software Architect)
   WHEN: Complex modules requiring detailed internal design - Data structures, algorithms, state machines - Error handling, concurrency, resource management
   ROLE: Software Architect perspective
   `ARTIFACT:MODULE-DESIGN(MODULE:...)`

   v. Hardware Abstraction Layer Specification (Firmware Architect)
   WHEN: Firmware/hardware interaction, embedded systems - Register definitions, driver APIs, hardware initialization - Platform-specific abstractions
   ROLE: Firmware Architect perspective
   `ARTIFACT:HAL-SPEC(HARDWARE:...)`

   vi. Deployment View (Systems Architect)
   WHEN: Infrastructure changes, new deployments, scaling - Infrastructure topology, scaling, failover - Software-to-hardware mapping
   ROLE: Systems Architect perspective
   SOURCE: arc42 Section 7, C4 Deployment Diagram
   `ARTIFACT:DEPLOYMENT(SCOPE:...)`

   vii. Solution Architecture Document (Solutions Architect)
   WHEN: Moderate-complex projects, integration patterns - Business capabilities mapped to technical components - Integration patterns, technology stack rationale
   ROLE: Solutions Architect perspective
   SOURCE: arc42 Section 4 - Solution Strategy
   `ARTIFACT:SAD-SOLUTION(SCOPE:...)`

   --- Architecture Decisions (ALWAYS) ---
   f. Skill: Architecture Decision Records

   - Technology choices, pattern selections, trade-offs
   - One ADR per significant architectural decision
     SOURCE: <https://github.com/joelparkerhenderson/architecture-decision-record>
     `ARTIFACT:ADR(DECISION:...)`

   ARTIFACT SELECTION MATRIX:

   | Complexity | Always      | Add Systems Arch  | Add Solutions Arch  | Add Software Arch             | Firmware  |
   | ---------- | ----------- | ----------------- | ------------------- | ----------------------------- | --------- |
   | Simple     | a,b,c.iii,f | -                 | -                   | c.iii, e.ii                   | -         |
   | Moderate   | a,b,c.iii,f | c.i (if external) | c.ii, e.vii         | c.iii, e.ii, e.iii            | -         |
   | Complex    | a,b,c.iii,f | c.i, e.i, e.vi    | c.ii, e.vii         | c.iii, e.ii, e.iii, e.iv      | -         |
   | Firmware   | a,b,c.iii,f | c.i, e.vi         | c.ii (if multi-sys) | c.iii, e.ii, e.iii, e.iv, e.v | e.v (HAL) |

   NOTE: "d" (MPI) created/updated when project infrastructure changes, regardless of complexity.

   ARTIFACT INTERDEPENDENCIES:

   - ARTIFACT:PRD → informs all architecture artifacts
   - ARTIFACT:NFR → constrains Solutions and Systems architecture
   - ARTIFACT:ARCH-CONTEXT → defines scope for ARCH-CONTAINER
   - ARTIFACT:ARCH-CONTAINER → defines boundaries for ARCH-COMPONENT
   - ARTIFACT:ARCH-COMPONENT → drives API-SPEC and MODULE-DESIGN
   - All architecture artifacts → referenced by ADRs for decision rationale
   - All architecture artifacts → consumed by Planning Coordination (step 3)

3. Skill(fork): Planning Coordination - The plan references what has or hasn't yet been created towards the architecture documents and each phase or step in the plan has a link to what part of which document it implements. This takes the change to the architecture artifacts and organises when it can be done in the sequence of existing items within the plan, and ensures that tthe tasks that are not overlapping are done concurrently, and the tasks that are overlapping in their scope or files they touch are done sequentially, and tasks are done in a way that future tasks can benefit from the creation of items from the previous task. Such as, if one task is creating shared utility functions, it would be done before a task that may need to use those utility functions.
4. Implementation (with git worktrees)
   a. Orchestrator agent manages the creation and cleanup of worktrees, and starting the Workers concurrently or in sequence based on the plan (not agents, independant claude instances started via the cli, or sdk)
5. Verification (multi-phase)
   └─→ Input: Implementation + Original Desired Outcome request
   Agent: Verification Agent
   Phases:
   1. Correctness (DoD, acceptance criteria)
   2. Regression (existing systems unaffected)
   3. Design Compliance (follows ADRs, patterns)
   4. Quality Standards (linting, tests, CI/CD)
   5. Documentation (runbooks, API docs complete)
      Output: ARTIFACT:VERIFICATION(SCOPE:...) (e.g. verification-report.md) with issue list
