# P1 Backlog System Redesign Github Issues As Source Of Truth Wit

| Property | Value |
|----------|-------|
| **Name** | P1 Backlog System Redesign Github Issues As Source Of Truth Wit |
| **Repository** | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/backlog/p1-backlog-system-redesign-github-issues-as-source-of-truth-wit.md) (⭐ 20) |
| **Original Path** | `.claude/backlog/p1-backlog-system-redesign-github-issues-as-source-of-truth-wit.md` |
| **Category** | content-creation |
| **Subcategory** | writing |
| **Tags** | content creation |
| **Created** | 2026-02-26 |
| **Updated** | 2026-02-27 |
| **File Hash** | `ba06b0202c95e8f7...` |

## Description

Current architecture has BACKLOG.md as primary source of truth with GitHub Issues as a secondary mirror. This is inverted — only works for one agent in one repo clone. A second session, different machine, or teammate sees stale markdown files. Redesign so GitHub Issues + Projects are the backend (available from anywhere) and local .claude/backlog/ files are a derived read cache rebuilt on demand. Key changes: (1) add creates GitHub Issue first, writes local cache from response, (2) list queries gh issue list with label filters, caches locally for speed, (3) update does gh issue edit first then updates cache, (4) close happens via Fixes #N in PR — GitHub auto-closes on merge, no manual skill workflow needed, (5) status lives in GitHub labels not markdown fields, (6) priority lives in GitHub labels not markdown section headers, (7) plan artifacts attached as issue comments or linked in body, (8) groomed content written into issue body expandable sections. Add pull/push commands: backlog.py pull fetches all open issues and rebuilds local cache, backlog.py push syncs local edits back to GitHub. The Fixes #N convention in commits/PRs handles closing automatically.

**Tags:** `content creation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/backlog/p1-backlog-system-redesign-github-issues-as-source-of-truth-wit.md)*
