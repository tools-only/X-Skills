# Context Compaction Survival

## What Is Compaction?

Long Claude Code sessions automatically compress earlier conversation to stay within the context window. This is invisible — there's no warning. After compaction, detailed reasoning and nuance from earlier in the conversation is lost.

## What Survives vs What's Lost

| Survives | Lost |
|----------|------|
| Files on disk (SESSION.md, CLAUDE.md) | Conversation details and reasoning |
| Git commits and their messages | Why a decision was made (unless written down) |
| Tool output summaries | Nuanced understanding built up over many exchanges |
| Recent messages (post-compaction) | Earlier trial-and-error context |

## Rules

1. **Write discoveries to files immediately** — if you find a gotcha, command, or pattern, add it to CLAUDE.md now, not "when we wrap up"
2. **Use sub-agents for heavy operations** — sub-agent context is separate from yours, so reading 20 files in a sub-agent doesn't consume your main context
3. **Checkpoint before context fills up** — don't wait for compaction to happen; checkpoint at milestones
4. **Git commits are your safety net** — WIP commits preserve exact file state even if SESSION.md is vague
5. **Make SESSION.md self-sufficient** — a fresh session reading only SESSION.md + CLAUDE.md should know exactly what to do next

## When to Checkpoint

- After completing a logical unit of work (a function, a component, a migration)
- Before attempting something risky (major refactor, dependency upgrade, destructive operation)
- When context feels large (many files read, long back-and-forth debugging)
- Before switching to a different area of the codebase

## Anti-Patterns

| Anti-Pattern | Why It Fails |
|-------------|--------------|
| "I'll update SESSION.md at the end" | Compaction happens mid-session — you lose the context to write a good summary |
| Vague resume instructions | Next session wastes time re-discovering what you already knew |
| SESSION.md over 100 lines | Defeats the purpose — becomes another doc to maintain instead of a quick handoff |
| Checkpointing every 5 minutes | Creates noise; checkpoint at milestones, not on a timer |
