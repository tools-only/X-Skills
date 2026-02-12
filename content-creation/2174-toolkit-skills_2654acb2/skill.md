# Skill Routing

You do not need explicit `/command` invocation. Auto-select based on task signals:

| Signal | Skill | Rationale |
|--------|-------|-----------|
| "experience the app", "find and fix all bugs", "QA and fix", discovery-first | /forge | Experience-first discovery → analysis → fix → validate |
| "fix", "broken", "debug", error context | /repair | Debugging router (auto-detects web vs mobile) |
| "fix the app", "debug production", "check staging" | /appfix | Autonomous app debugging |
| "fix mobile", "debug mobile app", "maestro test" | /mobileappfix | Mobile app debugging with Maestro |
| "build", "implement", complex task | /melt | Autonomous execution with verification |
| "clean up", "tech debt", "slop" | /burndown | Debt elimination |
| "lint", "enforce quality", "set up linting", "make repo clean" | /enforce | Deterministic quality enforcement |
| "improve design/UX/perf" | /improve | Recursive improvement loop |
| "analyze", "think deeply", "evaluate" | /heavy | Multi-perspective analysis |
| "create video", "promo video", "hero video", "AI cinema" | /video | AI-cinematic video generation |
| "create a loom", "record walkthrough", "screen recording" | /loom | Browser walkthrough with voiceover |
| "create audiobook", "turn into audio" | /audiobook | Document-to-audiobook generation |
| "generate episode", "educational video" | /episode | Minecraft-style educational videos |
| "write an essay", "essay about" | /essay | Publication-quality essay writing |
| "create plugin", "new plugin", "plugin lifecycle" | /plugin | SELD relay plugin authoring |
| "remember this", "capture this learning", "document solution" | /compound | Cross-session memory capture |
| "system health", "check health", "how is memory doing" | /health | Toolkit health diagnostics |
| "test harness", "test skill", "sandbox test" | /harness-test | Isolated harness testing |
| "prompt engineering", "write a skill", "optimize prompt" | /prompt-engineering-patterns | Context engineering |
| "terraform", "fix infra", "debug terraform" | /tfappfix | Terraform debugging |
| "auto mode", "autonomous mode", "enable hooks", "disable hooks" | /auto | Toggle hook enforcement on/off |
| No clear task / research only | No skill | Just answer directly |

## Parallelization Strategy

| Condition | Strategy |
|-----------|----------|
| Single focused task | Single-agent execution |
| 2 independent work items | Parallel `Task()` calls in a single message |
| 3+ independent work items with coordination needs | `TeamCreate` with shared task list + `SendMessage` |
| Fire-and-forget research or verification | `Task()` subagents (cheaper, simpler) |
| Agents need to share findings mid-work | `TeamCreate` with peer-to-peer `SendMessage` |

**Default to `Task()`.** Escalate to `TeamCreate` when:
- Agents need to share intermediate findings mid-research (e.g., /heavy Deep cross-pollination)
- 3+ execution items need coordination on shared resources (e.g., /melt worker teams)

**Trust skill-specific triage.** If a skill's complexity assessment selects TeamCreate, follow it — the skill has determined that peer communication will produce higher quality results. Do not override with cost concerns.

## Skill Fluidity

Skills are capabilities, not cages. If the task evolves, adapt:

- Building (/melt) and discover a critical bug? Fix it inline using /repair techniques. No formal mode switch needed.
- Debugging (/repair) and find the root cause is tech debt? Apply /burndown patterns to the area.
- Any task turns out to be architecturally complex? Use /heavy analysis for the sub-problem, then continue.
- Task decomposes into 3+ independent work items? Spawn a team (`TeamCreate`) inline. No mode switch needed.

The `autonomous-state.json` mode field drives auto-approval and checkpoint enforcement. It does not constrain your cognitive approach. Use the best technique for each sub-problem regardless of which skill activated the session.

When to formally re-invoke a skill (via Skill tool):
- The ENTIRE task has shifted (not just a sub-problem)
- You need the full activation ceremony of another skill

When to just adapt inline:
- A sub-problem needs a different approach
- You discovered something that changes the next step but not the overall goal
