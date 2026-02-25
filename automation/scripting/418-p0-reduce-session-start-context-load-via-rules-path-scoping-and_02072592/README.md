# P0 Reduce Session Start Context Load Via Rules Path Scoping And

| Property | Value |
|----------|-------|
| **Name** | P0 Reduce Session Start Context Load Via Rules Path Scoping And |
| **Repository** | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/backlog/p0-reduce-session-start-context-load-via-rules-path-scoping-and.md) (⭐ 20) |
| **Original Path** | `.claude/backlog/p0-reduce-session-start-context-load-via-rules-path-scoping-and.md` |
| **Category** | automation |
| **Subcategory** | scripting |
| **Tags** | automation |
| **Created** | 2026-02-24 |
| **Updated** | 2026-02-24 |
| **File Hash** | `02072592453247ba...` |

## Description

Session context is at 50% before any work begins. Root causes identified: (1) All custom agent descriptions load at session start (~5k tokens for 50+ agents) with no lazy loading mechanism. (2) CLAUDE.md files load unconditionally and in full (~18k tokens for memory files). (3) All plugin skills load descriptions at startup. Research confirmed these optimization primitives exist but are unused:

**Tags:** `automation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/backlog/p0-reduce-session-start-context-load-via-rules-path-scoping-and.md)*
