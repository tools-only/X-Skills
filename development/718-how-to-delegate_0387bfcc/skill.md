---
name: how-to-delegate
description: Scientific delegation framework for orchestrators - provide observations and success criteria while preserving agent autonomy
---

# Delegation Preparation Worksheet

**Workflow Reference**: See [Multi-Agent Orchestration](./../knowledge/workflow-diagrams/multi-agent-orchestration.md) for the complete delegation flow, DONE/BLOCKED signaling protocol, and agent selection guide.

**CRITICAL**: You are an orchestrator. You MUST complete ALL steps in this worksheet before invoking the Task tool. Incomplete preparation leads to failed delegations, wasted agent context, and poor outcomes.

<user_instructions>$ARGUMENTS</user_instructions>

---

## Conditional Execution

**If `<user_instructions>` is empty or contains only `$ARGUMENTS`:**

1. Complete Step 1 (Load Framework) NOW
2. State: "Delegation framework loaded. Awaiting task to delegate."
3. STOP - do not proceed to Step 2 until user provides a task

**If `<user_instructions>` contains a task:**

1. Complete ALL steps (1-10) before delegating
2. Fill in every section with actual content
3. Do NOT skip steps or summarize the worksheet

---

## Core Delegation Principles (Read First)

**Your role as orchestrator:**

- Route observations and context between user and agents
- Define measurable success criteria
- Enable comprehensive agent discovery
- NEVER pre-gather data that agents will collect themselves
- NEVER prescribe HOW agents should work (tools, steps, implementation)
- ALWAYS state WHAT must be achieved and WHY it matters

**Why this matters:**

- Agents have 200k context windows and specialized expertise
- Pre-gathering data causes context rot and duplicates work
- Prescribing HOW limits agents from discovering better solutions
- Observations without assumptions enable scientific method

---

## Step 1: Load Framework

**ACTION**: Activate the agent-orchestration skill NOW:

```text
Skill(command: "agent-orchestration")
```

**Then load domain-specific skills based on task type:**

```mermaid
flowchart TD
    START[Task Received] --> ORCH[Load agent-orchestration]
    ORCH --> CHECK{What domain?}

    CHECK -->|Python code| PY[python3-development]
    CHECK -->|Linting/code quality| LINT[holistic-linting]
    CHECK -->|GitLab CI| GL[gitlab-skill]
    CHECK -->|Git commits| CC[conventional-commits]
    CHECK -->|Package management| UV[uv]
    CHECK -->|Documentation| DOCS[mkdocs]
    CHECK -->|Pre-commit hooks| PRE[pre-commit]

    PY --> PROCEED[Proceed to Step 2]
    LINT --> PROCEED
    GL --> PROCEED
    CC --> PROCEED
    UV --> PROCEED
    DOCS --> PROCEED
    PRE --> PROCEED
```

**Why**: Domain skills contain specialized knowledge agents need. Loading them before delegation ensures agents have access to best practices and project conventions.

---

## Step 2: Identify Task Type

**ACTION**: Select the task type that best matches this delegation:

- [ ] **FOCUSED**: Single file, clear test, known location
- [ ] **INVESTIGATIVE**: Unknown cause, needs research
- [ ] **ARCHITECTURAL**: Multi-component, system-wide

**Why**: Task type determines how much context to provide. Focused tasks need precise location. Investigative tasks need all observations. Architectural tasks need system-wide context.

**TASK SUMMARY** (write one clear sentence):

```text
[Example: "Fix authentication failing for OAuth2 users" or "Investigate why CI pipeline times out on large PRs"]
```

---

## Step 3: Gather Observations

**ACTION**: List ONLY data already in your context. Write factual observations using "observed", "measured", "reported" language.

**CRITICAL**: Do NOT run commands to pre-gather data. Agents gather their own comprehensive data.

**Why**: Pre-gathering wastes your context, duplicates agent work, and causes context rot. Pass-through existing observations; let agents collect fresh data.

**OBSERVATIONS FROM USER:**

```text
[Example: "User reported: 'OAuth login redirects to 404'" or "User stated build fails on Python 3.12"]
```

**OBSERVATIONS FROM PRIOR AGENTS** (if any):

```text
[Example: "context-gathering agent found 3 instances of deprecated auth.login() in src/handlers/"]
```

**ERRORS ALREADY IN CONTEXT** (verbatim, if any):

```text
[Example: Exact error text you already received, not pre-gathered by running commands now]
```

**KNOWN LOCATIONS** (file:line references you already have):

```text
[Example: "src/auth/oauth.py:127 - where user reported issue occurs"]
```

---

## Step 4: Define Success

**ACTION**: Define specific, measurable outcomes and verification methods.

**Why**: Clear success criteria prevent scope creep and ensure agent knows when they're done.

**WHAT must be true when done** (measurable outcome):

```text
[Example: "OAuth login completes successfully for all providers" or "All pytest tests in test_auth.py pass"]
```

**HOW will completion be verified:**

```text
[Example: "Run `pytest test_auth.py -v` - all tests pass" or "Manual test: log in with Google/GitHub/Microsoft accounts"]
```

---

## Step 5: Provide World-Building Context

**ACTION**: Define WHERE to look, WHAT to achieve, and WHY it matters. Focus on context, not implementation.

**Why**: World-building enables agents to understand the problem space and make informed decisions about HOW to solve it.

**WHERE** (problem location, scope boundaries):

```text
[Example: "Authentication module at src/auth/ - OAuth handlers specifically" or "CI pipeline .github/workflows/test.yml"]
```

**WHAT** (identification criteria, acceptance criteria):

```text
[Example: "OAuth redirect must return 200 status with valid session token" or "Pipeline must complete within 10 minutes"]
```

**WHY** (expected outcomes, user requirements):

```text
[Example: "Users cannot log in with enterprise SSO accounts, blocking customer onboarding" or "Slow CI blocks PRs, reducing team velocity"]
```

---

## Step 6: Describe Available Resources

**ACTION**: Describe the ecosystem and available tools. List capabilities, don't prescribe which to use.

**Why**: Agents need to know what's available but should choose tools based on their expertise. Prescribing tools limits discovery.

**PROJECT ECOSYSTEM** (language, package manager, build system):

```text
[Example: "Python project using uv for all operations - activate uv skill" or "Node.js with pnpm workspaces"]
```

**AUTHENTICATED CLIS** (gh, glab, aws, etc.):

```text
[Example: "gh CLI pre-authenticated for GitHub operations" or "glab configured for GitLab access"]
```

**MCP TOOLS AVAILABLE** (check your functions list):

```text
[Example: "Excellent MCP servers installed - check <functions> list and prefer MCP tools (Ref, context7, exa) over built-in alternatives"]
```

**PROJECT-SPECIFIC RESOURCES** (scripts, reports, docs):

```text
[Example: "Validation scripts in ./scripts/ - check README.md" or "Previous fixes documented in .claude/reports/"]
```

---

## Step 7: Select Agent

**ACTION**: Choose the agent type that best matches the task domain and requirements.

**Why**: Specialized agents have domain expertise and optimized workflows for their task types.

**AGENT TYPE** (from available subagent_types):

```text
[Example: "python-cli-architect" or "linting-root-cause-resolver" or "context-gathering"]
```

**RATIONALE** (why this agent matches the task):

```text
[Example: "Task involves Python code changes - python-cli-architect has Python expertise and best practices" or "Need comprehensive context without polluting orchestrator context - context-gathering agent is optimized for this"]
```

---

## Step 8: Pre-Flight Verification

**ACTION**: Review your filled sections against these criteria. Mark each with ✓ or ✗.

**Why**: This checklist catches common anti-patterns that limit agent effectiveness.

| Check                                          | Status | Fix If ✗                                  |
| ---------------------------------------------- | ------ | ----------------------------------------- |
| Will start with `Your ROLE_TYPE is sub-agent.` | \_\_\_ | Add to Step 9 prompt                      |
| Contains only factual observations             | \_\_\_ | Remove "I think", "probably", "likely"    |
| No assumptions stated as facts                 | \_\_\_ | Replace with "Hypothesis to verify:"      |
| Defines WHAT, not HOW                          | \_\_\_ | Remove tool names, steps, line numbers    |
| Lists resources, doesn't prescribe tools       | \_\_\_ | Change "Use X" to "X is available"        |
| No pre-gathered data (agent will gather own)   | \_\_\_ | Remove command outputs you just collected |

---

## Step 9: Construct Task Prompt

**ACTION**: Use your filled-in sections to construct the delegation prompt following this template.

**Why**: This template structure ensures agents receive observations, success criteria, and autonomy to apply their expertise.

**Copy this template and fill in from your worksheet:**

```text
Your ROLE_TYPE is sub-agent.

[Task summary from Step 2]

OBSERVATIONS:
[From Step 3 - verbatim, not paraphrased]

DEFINITION OF SUCCESS:
[From Step 4]

CONTEXT:
[From Step 5 - WHERE, WHAT, WHY]

YOUR TASK:
1. Run SlashCommand /is-it-done to understand completion criteria
2. Perform comprehensive context gathering using available tools, skills, and resources
3. Form hypothesis based on evidence
4. Design and execute experiments
5. Verify findings against authoritative sources
6. Implement solution following best practices
7. Verify /is-it-done criteria satisfied with evidence

AVAILABLE RESOURCES:
[From Step 6 - describe ecosystem, don't prescribe tools]
```

---

## Step 10: Delegate

**ACTION**: Final check, then delegate.

**READY TO DELEGATE?** Mark ✓ or ✗:

- [ ] All verification checks in Step 8 passed
- [ ] Prompt constructed from Step 9 template
- [ ] Agent type selected from Step 7
- [ ] No pre-gathered data included
- [ ] Defines WHAT and WHY, not HOW

**If all ✓, invoke the Task tool:**

```text
Task(
  agent="[agent type from Step 7]",
  prompt="[your constructed prompt from Step 9]"
)
```

**Why final check matters**: One incomplete step can cause agent confusion, wasted context, or failed delegation. Taking 30 seconds to verify saves 10 minutes of back-and-forth.
