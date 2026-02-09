# Skill Routing

You do not need explicit `/command` invocation. Auto-select based on task signals:

| Signal | Skill | Rationale |
|--------|-------|-----------|
| "fix", "broken", "debug", error context | /repair | Debugging router (auto-detects web vs mobile) |
| "build", "implement", complex task | /melt | Autonomous execution with verification |
| "clean up", "tech debt", "slop" | /burndown | Debt elimination |
| "improve design/UX/perf" | /improve | Recursive improvement loop |
| "analyze", "think deeply", "evaluate" | /heavy | Multi-perspective analysis |
| "create video", "promo video", "hero video", "AI cinema" | /video | AI-cinematic video generation |
| "good morning", "briefing", "catch me up", "what did I miss" | /briefing | Morning digest across all platforms |
| "check discord", "discord messages", "discord DMs" | /discord-digest | Discord activity summary |
| "check X/twitter", "timeline", "bookmarks", "what's on X" | /x-digest | X timeline + Grok analysis |
| "reply on discord", "respond in channel", "send discord" | /discord-reply | Context-aware Discord reply |
| "tweet", "post on X", "draft thread" | /x-post | Voice-consistent tweet drafting |
| "engage on X", "find conversations", "grow audience" | /x-engage | Strategic X engagement |
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
