# AI Agents and Skills Workshop for Embedded Engineers

**Duration:** 1 hour
**Audience:** Embedded Development Engineers, Firmware QA Engineers
**Prerequisites:** Cursor installed, basic familiarity with terminal commands

---

## Workshop Agenda

| Time | Topic |
|------|-------|
| 0:00-0:08 | Introduction: Why AI Assistants Need Domain Expertise |
| 0:08-0:16 | Skills: Extending Cursor Agent with Specialized Knowledge |
| 0:16-0:24 | Subagents: Creating Focused AI Workers |
| 0:24-0:32 | Orchestration: Skills + Subagents for Multi-Step Work |
| 0:32-0:40 | Hands-On: Build a Firmware Code Review Skill |
| 0:40-0:48 | Hands-On: Build an Embedded Debugging Subagent |
| 0:48-0:55 | Installing Community Skills |
| 0:55-1:00 | Resources and Next Steps |

---

## Part 1: Introduction (10 minutes)

### Why Embedded Engineers Need Custom AI Tools

Cursor’s Agent is powered by strong general-purpose models, but embedded development has unique challenges:

- **Memory constraints** - Code must fit in kilobytes, not gigabytes
- **Timing-critical operations** - Interrupt handlers, real-time requirements
- **Hardware interactions** - Register manipulation, peripheral drivers
- **Cross-compilation** - Building for ARM, RISC-V, or custom architectures
- **Testing difficulties** - Hardware-in-the-loop, JTAG debugging, protocol analysis

**The Problem:** General AI assistants don't understand your specific:
- Register maps and memory layouts
- RTOS APIs (FreeRTOS, Zephyr, ThreadX)
- Coding standards (MISRA-C, CERT-C, company guidelines)
- Build systems (CMake cross-compilation, vendor toolchains)

**The Solution:** Skills and subagents let you teach Cursor Agent your domain expertise.

### Cursor Extensibility Model (Context + Skills + Subagents)

```
┌─────────────────────────────────────────────────────────────┐
│                         Cursor                               │
├─────────────────────────────────────────────────────────────┤
│  Rules              │  Skills           │  Subagents         │
│  (.cursor/rules)    │  (.cursor/skills) │  (.cursor/agents)  │
├─────────────────────────────────────────────────────────────┤
│  Always loaded      │  Loaded on        │  Delegated to     │
│  Every conversation │  demand           │  for specific     │
│                     │                   │  tasks            │
└─────────────────────────────────────────────────────────────┘
```

| Extension Type | When to Use | Example |
|----------------|-------------|---------|
| **Rules** | Persistent instructions for a project | “Use MISRA-C 2012 guidelines” |
| **Skills** | Portable workflows / knowledge packages | Register map docs, review checklists |
| **Subagents** | Focused workers for multi-step tasks or isolation | Codebase exploration, crash analysis |

---

## Part 2: Skills - Extending Cursor Agent with Specialized Knowledge (10 minutes)

### What is a Skill?

A skill is a folder containing a `SKILL.md` (plus optional `references/`, `scripts/`, `assets/`) that teaches Cursor Agent specialized knowledge or workflows. Think of it as an “onboarding guide” that transforms a generalist assistant into a domain expert.

### Skill Directory Structure

```
firmware-review/
├── SKILL.md              # Required: Instructions for Cursor Agent
├── references/           # Optional: Detailed documentation
│   ├── misra-rules.md
│   └── register-map.md
├── scripts/              # Optional: Automation scripts
│   └── check_stack.py
└── assets/               # Optional: Templates, configs
    └── review-template.md
```

### Skill Locations

| Location | Scope | Use Case |
|----------|-------|----------|
| `.cursor/skills/` | Project-level | Skills stored in a repo and shared with the team |
| `.claude/skills/` | Project-level | Compatibility directory (Cursor also discovers these) |
| `.codex/skills/` | Project-level | Compatibility directory (Cursor also discovers these) |
| `~/.cursor/skills/` | User-level | Personal skills available across projects |
| `~/.claude/skills/` | User-level | Compatibility directory |
| `~/.codex/skills/` | User-level | Compatibility directory |

### What Cursor does at startup (important)

When Cursor starts, it **automatically discovers skills** from the directories above and presents them to Agent. Agent can also **apply skills automatically** when it decides they’re relevant, or you can **manually invoke** skills by typing `/` in chat and searching for the skill name.

### Required frontmatter fields (per Cursor docs)

In Cursor, `SKILL.md` YAML frontmatter supports:

- `name` (**required**): lowercase letters/numbers/hyphens only; **must match the parent folder name**
- `description` (**required**): what + when; used for relevance
- Optional fields: `license`, `compatibility`, `metadata`, `disable-model-invocation`

### Disable auto-loading (command-like skills)

By default, skills can be applied automatically when Agent judges them relevant. If you want a skill to behave like an explicit slash command (only included when the user runs it), set:

```yaml
disable-model-invocation: true
```

### Optional skill directories (progressive disclosure)

Cursor skills support these optional folders:

- `scripts/`: executable helpers (shell, Python, etc.)
- `references/`: documentation loaded on demand
- `assets/`: templates/configs/static artifacts

Keep `SKILL.md` focused; push detail into `references/` so it’s only loaded when needed.

### SKILL.md Anatomy

```yaml
---
name: firmware-review
description: 'Review embedded C code for MISRA-C compliance, memory safety, and interrupt handling issues. Use when reviewing firmware code, analyzing drivers, or auditing safety-critical embedded systems.'
---

# Firmware Code Review

## Review Checklist

When reviewing embedded C code, check for:

1. **Memory Safety**
   - Stack usage analysis
   - Buffer boundary validation
   - Pointer arithmetic correctness

2. **Interrupt Safety**
   - Volatile qualifiers on shared variables
   - Critical section protection
   - ISR execution time

3. **MISRA-C Compliance**
   - See [MISRA Rules Reference](./references/misra-rules.md)

## Common Patterns

[Example patterns specific to your codebase...]
```

### Key Frontmatter Fields

| Field | Purpose | Example |
|-------|---------|---------|
| `name` | Skill identifier (must match folder name) | `firmware-review` |
| `description` | Tells Cursor when to apply the skill | Include trigger keywords |
| `disable-model-invocation` | Only include when user explicitly runs it | `true` for “command-like” skills |

### Best Practices for Skills

**1. Write Excellent Descriptions**

```yaml
# BAD - Too vague, won't trigger correctly
description: Helps with firmware

# GOOD - Clear triggers and use cases
description: 'Review embedded C code for MISRA-C compliance, memory safety, and interrupt handling. Use when reviewing STM32 drivers, analyzing FreeRTOS tasks, or auditing safety-critical firmware.'
```

**2. Use Progressive Disclosure**

Keep `SKILL.md` under 500 lines. Move detailed reference material to `references/` files:

```markdown
# SKILL.md - Keep lean
For complete MISRA-C rule details, see [MISRA Reference](./references/misra-rules.md)
```

**3. Include Concrete Examples**

```markdown
## Stack Analysis Example

**Input:** Function with deep call chain
**Output:**
- Maximum stack depth: 2,048 bytes
- Deepest path: main() -> process_can() -> parse_frame() -> update_state()
- Recommendation: Reduce recursion in parse_frame()
```

**4. Match Freedom to Fragility**

| Task Type | Freedom Level | Example |
|-----------|---------------|---------|
| Code analysis | High | Text guidance, Claude decides approach |
| Flash programming | Low | Specific script with exact parameters |
| Build validation | Medium | Pseudocode with configurable options |

---

## Part 3: Subagents - Creating Focused AI Workers (10 minutes)

### What is an Agent?

In Cursor, a **subagent** is a specialized assistant that Agent can delegate to. Subagents run with their own context window and return a summary to the parent Agent.

### Foreground vs background (per Cursor docs)

Subagents can run:

- **Foreground**: parent Agent waits; best when you need the output to proceed
- **Background**: parent Agent continues immediately; best for long-running tasks or parallel workstreams

### Built-in subagents (per Cursor docs)

Cursor includes built-in subagents designed to keep noisy output out of the main conversation:

- **Explore**: searches/analyzes codebases (fast, parallel)
- **Bash**: runs shell command sequences (isolates verbose logs)
- **Browser**: controls browser via MCP tools (filters noisy DOM/screenshot output)

### Agent vs Skill: When to Use Which

| Feature | Skill | Subagent |
|---------|-------|-------|
| **Purpose** | Knowledge/workflows | Focused task execution |
| **Invocation** | User runs it via `/` or Agent applies it | Agent delegates when helpful |
| **Context** | In the main Agent context | Runs in its own context window |
| **Best for** | Repeatable checklists + reference | Long investigations, parallel workstreams |
| **Best for** | Reference material, guidelines | Code review, test generation, debugging |

### Subagent File Structure

Subagents are single markdown files in `.cursor/agents/` or `~/.cursor/agents/`:

```
.cursor/agents/
├── firmware-reviewer.md
├── test-generator.md
└── debug-analyzer.md
```

### Subagent Anatomy

```yaml
---
name: firmware-reviewer
description: 'Review embedded C/C++ code for memory safety, MISRA compliance, and interrupt handling. Use when reviewing drivers, HAL code, or RTOS tasks. Focuses on stack usage, volatile correctness, and timing analysis.'
---

# Firmware Code Reviewer

You are a senior embedded systems engineer specializing in safety-critical firmware review.

## Core Competencies

<competencies>
- Memory safety analysis (stack, heap, static allocation)
- Interrupt and concurrency correctness
- MISRA-C 2012 compliance verification
- Peripheral driver review (DMA, timers, ADC/DAC)
- RTOS task analysis (FreeRTOS, Zephyr)
</competencies>

## Review Workflow

<workflow>
1. Identify file type (driver, application, RTOS task, ISR)
2. Check memory safety patterns
3. Verify interrupt handling correctness
4. Assess MISRA-C compliance
5. Generate prioritized findings
</workflow>

## Output Format

```markdown
# Code Review: [filename]

## Critical Issues
- [Issue with file:line reference]

## Warnings
- [Issue with file:line reference]

## Recommendations
- [Improvement suggestion]

## Positive Findings
- [What was done well]
```
```

### Key Subagent Frontmatter Fields

| Field | Purpose | Values |
|-------|---------|--------|
| `name` | Unique identifier | lowercase, hyphens, max 64 chars |
| `description` | When to delegate to this subagent | Include action verbs and triggers |

Cursor subagents also support (optional) config fields:

- `model`: `fast`, `inherit`, or a specific model ID (defaults to inherit)
- `readonly`: `true` to restrict write permissions
- `is_background`: `true` to run without waiting

### Best Practices for Subagents

**1. Write Trigger-Rich Descriptions**

```yaml
# BAD
description: Reviews code

# GOOD
description: 'Review embedded C/C++ code for memory safety, MISRA compliance, and interrupt handling. Use when reviewing STM32 drivers, FreeRTOS tasks, or safety-critical firmware. Specializes in stack analysis and volatile correctness.'
```

**2. Keep responsibilities narrow**

- Good: “Analyze Cortex-M HardFault register dumps”
- Bad: “Do anything embedded”

**3. Make outputs predictable**

Add an “Output Format” section in the subagent prompt so the parent Agent gets consistent artifacts.

---

## Part 4: Orchestration - Skills as Commands with Agent Delegation (8 minutes)

### The power pattern in Cursor: Skill → (optional) Subagent → Skill

The most powerful Cursor pattern combines all three concepts:

1. **Skill** provides a durable workflow/checklist (e.g., `/add-embedded-feature`)
2. Agent can **delegate to subagents** for long research/analysis steps
3. Agent uses other **skills** as supporting references (e.g., MISRA or debug tools)

```text
User runs /firmware-release
        |
        v
firmware-release skill instructions
  - delegate “build validation” to a subagent
  - delegate “test execution” to a subagent
  - keep the main Agent focused on decisions + summary
```

### How orchestration works in practice

Cursor skills are **instructions**. Orchestration is done by the Agent following those instructions, optionally spawning subagents as needed.

When you want isolation/parallelism, explicitly ask for delegation in the instructions or your prompt, e.g.:

```text
Delegate call-graph tracing and dependency mapping to a subagent.
Return a short summary plus a list of files/entry points discovered.
```

### Example: Firmware Release Orchestrator

Create `.cursor/skills/firmware-release/SKILL.md`:

```yaml
---
name: firmware-release
description: 'Execute complete firmware release workflow: build validation, test execution, binary packaging, and release notes generation. Use when preparing a firmware release, creating release candidates, or packaging binaries for production.'
---

# Firmware Release Orchestrator

Execute the complete firmware release workflow for version $ARGUMENTS.

## Orchestration Rules

<orchestration_rules>
You are an orchestrator. You coordinate work across specialized agents.
You MUST delegate, not execute directly.

| Task | Delegate To | Never Do Directly |
|------|-------------|-------------------|
| Build validation | @build-validator | Run cmake yourself |
| Test execution | @test-runner | Run pytest yourself |
| Binary packaging | release scripts | Package binaries yourself |
| Release notes | @doc-generator | Write notes yourself |
</orchestration_rules>

## Phase 1: Build Validation

Delegate build validation to a subagent:

- Validate the build for all targets
- Verify binary sizes within flash constraints
- Check no undefined symbols
- Validate linker script memory layout

## Phase 2: Test Execution

Delegate test execution to a subagent:

- Unit tests (host simulation)
- Integration tests (hardware-in-loop if available)
- Static analysis (MISRA compliance)
- Report coverage metrics

## Phase 3: Binary Packaging

Run your project’s packaging script (exact command depends on the repo; keep it in the task plan as a concrete step).

## Phase 4: Release Notes

Delegate release note generation to a subagent:

- Parse git log since last tag
- Categorize: features, fixes, breaking changes
- Include binary checksums
- Format as a CHANGELOG entry

## Success Output

When all phases complete, provide:
- Build status for each target
- Test results summary
- Binary locations and checksums
- Release notes preview
```

### Calling Other Skills from Within a Skill

Skills can reference other skills in their instructions (as explicit `/skill-name` invocations):

```markdown
## Workflow

1. First, invoke /misra-review to check code compliance
2. Then, invoke /stack-analyzer to verify memory usage
3. Finally, run /test-coverage to ensure adequate testing

Each skill activation loads specialized knowledge for that domain.
```

### Embedded Engineering Orchestration Ideas

| Orchestrator Skill | Delegates To | Purpose |
|--------------------|--------------|---------|
| `/firmware-audit` | @misra-reviewer, @security-scanner, @stack-analyzer | Complete firmware quality audit |
| `/driver-review` | @peripheral-expert, @dma-validator, @timing-analyzer | Hardware driver code review |
| `/release-prep` | @build-validator, @test-runner, @packager | Release preparation workflow |
| `/debug-session` | @crash-analyzer, @memory-debugger, @trace-decoder | Structured debugging workflow |

### Best Practices for Orchestration

**1. Orchestrators delegate, never execute**

```markdown
# GOOD
Delegate “run the test suite” to a subagent and ask it to return a short, structured report.

# BAD
Dump long raw logs into the main conversation context.
```

**2. Each subagent has focused responsibility**

- Build validation subagent: build + size checks
- Test execution subagent: run tests + summarize failures
- Crash analysis subagent: decode dumps + likely causes

**4. Use artifacts for handoff between phases**

```markdown
Phase 1 output: plan/build-report.md
Phase 2 input: Read plan/build-report.md for context
```

---

## Part 5: Hands-On - Build a Firmware Code Review Skill (8 minutes)

### Exercise: Create a MISRA-C Review Skill

We'll create a skill that helps Claude review code for common embedded C issues.

**Step 1: Create the skill directory**

```bash
mkdir -p .cursor/skills/firmware-review/references
```

**Step 2: Create SKILL.md**

Create `.cursor/skills/firmware-review/SKILL.md`:

```yaml
---
name: firmware-review
description: 'Review embedded C code for memory safety, interrupt correctness, and MISRA-C compliance. Use when reviewing firmware, drivers, HAL code, or RTOS tasks. Triggers on: code review, MISRA, embedded, firmware, interrupt, volatile, stack.'
---

# Firmware Code Review Skill

Review embedded C/C++ code systematically for common issues in safety-critical systems.

## Review Categories

### 1. Memory Safety

Check for:
- **Stack overflow risk** - Deep call chains, large local arrays, recursive functions
- **Buffer overflows** - Array bounds, string operations, memcpy sizes
- **Pointer safety** - Null checks, alignment, arithmetic correctness
- **Static allocation** - Prefer static over dynamic in embedded systems

Red flags:
- `malloc()` / `free()` in firmware (memory fragmentation risk)
- Variable-length arrays (VLA) - stack size unpredictable
- Recursive functions without depth limits

### 2. Interrupt Safety

Check for:
- **Volatile qualifiers** - Variables shared between ISR and main code MUST be volatile
- **Critical sections** - Shared data access must be protected
- **ISR execution time** - Keep ISRs short, defer work to tasks
- **Reentrancy** - ISR-called functions must be reentrant

Pattern to verify:
```c
// CORRECT: volatile for ISR-shared variable
volatile uint32_t tick_count;

void SysTick_Handler(void) {
    tick_count++;  // Written in ISR
}

void delay_ms(uint32_t ms) {
    uint32_t start = tick_count;  // Read in main
    while ((tick_count - start) < ms);
}
```

### 3. MISRA-C Critical Rules

High-priority rules for embedded:

| Rule | Description | Why It Matters |
|------|-------------|----------------|
| 11.3 | No pointer/integer casts | Alignment faults, undefined behavior |
| 14.3 | Control expressions must be boolean | Logic errors |
| 17.7 | Function return values used | Missed error handling |
| 21.3 | No dynamic memory (malloc) | Fragmentation, determinism |

For complete rules, see [MISRA Reference](./references/misra-rules.md)

### 4. Hardware Interaction

Check for:
- **Register access** - Use volatile pointers, proper barriers
- **Bit manipulation** - Clear before set, read-modify-write protection
- **Timing** - Required delays after peripheral configuration

## Output Format

Provide findings in this format:

```markdown
## Critical Issues
- [MEMORY] file.c:123 - Buffer overflow in parse_packet()
- [ISR] driver.c:45 - Missing volatile on shared_flag

## Warnings
- [MISRA-11.3] hal.c:89 - Integer to pointer cast

## Recommendations
- Consider static buffer pool instead of malloc()

## Positive Findings
- Good use of const correctness throughout
```
```

**Step 3: Create a reference file**

Create `.cursor/skills/firmware-review/references/misra-rules.md`:

```markdown
# MISRA-C 2012 Quick Reference

## Memory and Pointers

| Rule | Title | Severity |
|------|-------|----------|
| 11.3 | Cast between pointer and integer | Required |
| 11.4 | Cast to pointer of different type | Advisory |
| 11.5 | Cast from void* to object pointer | Advisory |
| 18.1 | Pointer arithmetic bounds | Required |
| 18.2 | Pointer subtraction | Required |

## Control Flow

| Rule | Title | Severity |
|------|-------|----------|
| 14.3 | Controlling expressions boolean | Required |
| 15.1 | goto shall not be used | Advisory |
| 15.4 | Single loop exit point | Advisory |
| 16.2 | Switch top-level shall be case | Required |

## Functions

| Rule | Title | Severity |
|------|-------|----------|
| 17.2 | Recursion shall not be used | Required |
| 17.7 | Return value shall be used | Required |
| 8.13 | Pointer parameters as const | Advisory |

## Standard Library

| Rule | Title | Severity |
|------|-------|----------|
| 21.3 | Memory functions (malloc) | Required |
| 21.6 | I/O library stdio.h | Required |
| 21.7 | String functions (atof, atoi) | Required |
```

**Step 4: Test the skill**

```bash
# Invoke the skill
/firmware-review

# Or ask Claude naturally
"Review the interrupt handler in src/drivers/uart.c for safety issues"
```

---

## Part 6: Hands-On - Build an Embedded Debugging Subagent (8 minutes)

### Exercise: Create a Debug Analysis Agent

We'll create an agent that helps analyze firmware crashes and debugging output.

**Step 1: Create the agent file**

Create `.cursor/agents/firmware-debugger.md`:

```yaml
---
name: firmware-debugger
description: 'Analyze firmware crashes, hard faults, and debugging output. Use when investigating HardFault, stack traces, memory corruption, or crash dumps. Specializes in ARM Cortex-M fault analysis, FreeRTOS crash debugging, and memory dump interpretation.'
---

# Firmware Debug Analyzer

You are an expert embedded systems debugger specializing in ARM Cortex-M crash analysis and RTOS debugging.

## Core Competencies

<competencies>
- ARM Cortex-M HardFault analysis (CFSR, HFSR, MMFAR, BFAR)
- Stack trace reconstruction from memory dumps
- FreeRTOS task state analysis
- Memory corruption pattern recognition
- Register state interpretation
</competencies>

## Analysis Workflow

<workflow>

### Step 1: Identify Crash Type

Determine fault category from registers:

| Fault Type | HFSR Bit | CFSR Bits | Common Cause |
|------------|----------|-----------|--------------|
| Hard Fault | FORCED | Various | Escalated fault |
| MemManage | - | IACCVIOL, DACCVIOL | MPU violation |
| Bus Fault | - | PRECISERR, IMPRECISERR | Invalid address |
| Usage Fault | - | UNDEFINSTR, INVSTATE | Illegal instruction |

### Step 2: Analyze Fault Address

- **MMFAR** - Memory management fault address
- **BFAR** - Bus fault address
- Check against linker map for symbol resolution

### Step 3: Reconstruct Stack Trace

For Cortex-M exception frame:
```
SP+0x00: R0
SP+0x04: R1
SP+0x08: R2
SP+0x0C: R3
SP+0x10: R12
SP+0x14: LR (return address)
SP+0x18: PC (faulting instruction)
SP+0x1C: xPSR
```

### Step 4: Check Task State (FreeRTOS)

If using FreeRTOS, analyze:
- Current task TCB
- Stack high water mark
- Task state (Ready, Blocked, Suspended)

</workflow>

## Common Fault Patterns

### Pattern: Stack Overflow

**Symptoms:**
- HardFault with FORCED=1
- SP points below task stack base
- Memory corruption in adjacent task stack

**Investigation:**
1. Check task stack size in FreeRTOS config
2. Analyze call depth with `-fstack-usage`
3. Look for large local arrays or recursion

### Pattern: Null Pointer Dereference

**Symptoms:**
- Bus fault with PRECISERR=1
- BFAR near 0x00000000
- Often during peripheral access

**Investigation:**
1. Check uninitialized pointers
2. Verify peripheral clock enabled
3. Check handle initialization order

### Pattern: Unaligned Access

**Symptoms:**
- Usage fault with UNALIGNED=1
- Accessing packed struct members
- Casting byte buffers to larger types

**Investigation:**
1. Check struct packing
2. Look for `(uint32_t*)buffer` casts
3. Verify DMA buffer alignment

## Output Format

Provide analysis in this format:

```markdown
# Crash Analysis Report

## Fault Summary
- **Type:** [HardFault/BusFault/etc.]
- **Faulting Address:** 0x????????
- **Faulting Instruction:** [Symbol + offset]

## Root Cause Analysis
[Detailed explanation of what caused the fault]

## Evidence
- [Register values and interpretation]
- [Stack trace]
- [Memory state]

## Recommended Fix
1. [Specific action to take]
2. [Code changes needed]

## Prevention
- [How to prevent similar issues]
```

## Important Notes

<rules>
- Always correlate addresses with .map file symbols
- Check for interrupt/task context at fault time
- Consider peripheral state (clocks, power)
- Look for race conditions with ISRs
</rules>
```

**Step 2: Test the agent**

```bash
# Claude will delegate to this agent automatically when appropriate
"Analyze this crash dump from the STM32:
HFSR: 0x40000000
CFSR: 0x00000400
PC: 0x080012a4
LR: 0x08001287
SP: 0x20001ff0"
```

---

## Part 7: Installing Community Skills (7 minutes)

### The Awesome Agent Skills Repository

The community maintains a curated list of skills at:
**<https://github.com/VoltAgent/awesome-agent-skills>**

### Skills Relevant for Embedded Engineers

| Skill | Source | Use Case |
|-------|--------|----------|
| `trailofbits/constant-time-analysis` | Trail of Bits | Detect timing side-channels in crypto code |
| `trailofbits/dwarf-expert` | Trail of Bits | DWARF debugging format expertise |
| `trailofbits/static-analysis` | Trail of Bits | Static analysis with CodeQL, Semgrep |
| `trailofbits/property-based-testing` | Trail of Bits | Property-based test generation |
| `obra/test-driven-development` | Community | TDD methodology |
| `obra/systematic-debugging` | Community | Methodical debugging approach |

### Installing a Skill

**Method 1: Clone directly to skills directory**

```bash
# Clone a skill to your user skills directory
cd ~/.cursor/skills
git clone https://github.com/trailofbits/skills/tree/main/plugins/static-analysis

# Or for project-specific
cd your-project/.cursor/skills
git clone <skill-repo-url>
```

**Method 2: Manual download**

```bash
# Download the SKILL.md and reference files
mkdir -p ~/.cursor/skills/static-analysis
curl -o ~/.cursor/skills/static-analysis/SKILL.md \
  https://raw.githubusercontent.com/trailofbits/skills/main/plugins/static-analysis/SKILL.md
```

**Method 3: Import from GitHub via Cursor Settings**

Use Cursor Settings → Rules → Add Rule → Remote Rule (GitHub) and point it at a repository that contains skills. Cursor will sync the content into your project.

### Verifying Installation

```bash
# Type `/` in chat and search for the skill name.
# You can also open Cursor Settings → Rules to view discovered skills.
```

### Viewing skills in Cursor (per Cursor docs)

1. Open **Cursor Settings** (Cmd+Shift+J on Mac, Ctrl+Shift+J on Windows/Linux)
2. Go to **Rules**
3. Discovered skills show up in the **Agent Decides** section

### Importing skills from GitHub (per Cursor docs)

Cursor can import skills via:

1. **Cursor Settings → Rules**
2. In **Project Rules**, click **Add Rule**
3. Select **Remote Rule (Github)**
4. Paste the repository URL (Cursor will sync it)

### Migrating rules/commands to skills (Cursor 2.4+)

Cursor includes a built-in `/migrate-to-skills` skill that converts:

- **Dynamic rules** (Apply Intelligently, no globs, not always-apply) → standard skills
- **Slash commands** → skills with `disable-model-invocation: true`

---

## Part 8: Resources and Next Steps (5 minutes)

### Official Documentation

| Resource | URL |
|----------|-----|
| Cursor Skills | <https://cursor.com/docs/context/skills> |
| Cursor Subagents | <https://cursor.com/docs/context/subagents> |
| Cursor Rules | <https://cursor.com/docs/context/rules> |
| Agent Skills Standard | <https://agentskills.io> |

### Community Resources

| Resource | Description |
|----------|-------------|
| [awesome-agent-skills](https://github.com/VoltAgent/awesome-agent-skills) | Curated list of 172+ skills |
| [Anthropic Skills Repo](https://github.com/anthropics/skills) | Official example skills |
| [Trail of Bits Skills](https://github.com/trailofbits/skills) | Security-focused skills |

### Embedded-Relevant External Skills

Consider adapting these for your needs:

1. **trailofbits/constant-time-analysis** - Critical for secure bootloaders and crypto
2. **trailofbits/dwarf-expert** - Debug symbol analysis
3. **trailofbits/static-analysis** - Integrate with your CI pipeline
4. **obra/systematic-debugging** - Methodical debugging approach

### Ideas for Your Own Skills

| Skill Idea | Description |
|------------|-------------|
| `register-map-reader` | Parse and explain your MCU's register definitions |
| `freertos-analyzer` | Analyze task priorities, stack usage, deadlock risks |
| `can-protocol-helper` | CAN message parsing and DBC file interpretation |
| `memory-layout-validator` | Verify linker script against hardware constraints |
| `timing-analyzer` | Analyze interrupt latency and WCET |

### Ideas for Your Own Agents

| Agent Idea | Description |
|------------|-------------|
| `hal-migrator` | Help migrate between HAL versions or MCU families |
| `test-generator` | Generate unit tests for embedded code |
| `power-analyzer` | Review code for power consumption issues |
| `safety-auditor` | Audit code for IEC 61508 / ISO 26262 compliance |

### Next Steps

1. **Create your first skill** - Start with something you explain repeatedly
2. **Build a review agent** - Automate your code review checklist
3. **Share with your team** - Commit skills to `.cursor/skills/` in your repo
4. **Iterate based on usage** - Refine skills as you discover gaps

---

## Quick Reference Card

### Skill Structure

```
skill-name/
├── SKILL.md              # Required
├── references/           # Optional: detailed docs
├── scripts/              # Optional: automation
└── assets/               # Optional: templates
```

### Skill Frontmatter (Basic)

```yaml
---
name: skill-name
description: 'What + when. Include trigger keywords.'
---
```

### Skill “command-like” behavior (optional)

If you want the skill to only be included when the user explicitly runs it, set:

```yaml
disable-model-invocation: true
```

### Subagent Structure

Single file: `.cursor/agents/agent-name.md`

### Subagent Frontmatter

```yaml
---
name: agent-name
description: 'What + when. Include trigger keywords.'
---
```

### Orchestration Patterns

| Pattern | Syntax | Use |
|---------|--------|-----|
| Skill provides workflow | “Follow the steps in `/add-embedded-feature`” | Repeatable process |
| Agent delegates to subagent | “Delegate this investigation to a subagent” | Isolation + parallelism |
| Skill references files | Link to `references/*.md` | Progressive disclosure |

### Delegation Syntax

```markdown
In a skill or prompt, explicitly ask Cursor Agent to delegate:

- “Delegate call-graph tracing to a subagent; return entry points + file references.”
- “Delegate crash dump interpretation to a subagent; return fault type + evidence.”
```

---

## Workshop Complete

You now know how to:
- Create skills to extend Cursor with domain expertise
- Build subagents for focused task execution
- **Orchestrate multi-step workflows using skills + subagents**
- Install community skills from awesome-agent-skills
- Apply best practices for descriptions and progressive disclosure

**Questions?** Check the resources above or experiment with your own skills and agents.
