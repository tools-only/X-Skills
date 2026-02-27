---
name: Reduce session-start context load via rules path-scoping and disable-model-invocation
description: "Session context is at 50% before any work begins. Root causes identified: (1) All custom agent descriptions load at session start (~5k tokens for 50+ agents) with no lazy loading mechanism. (2) CLAUDE.md files load unconditionally and in full (~18k tokens for memory files). (3) All plugin skills load descriptions at startup. Research confirmed these optimization primitives exist but are unused: 'paths' frontmatter in .claude/rules/*.md loads rule only when working with matching file patterns; 'd"
metadata:
  topic: reduce-session-start-context-load-via-rules-path-scoping-and
  source: Session observation — context at 50% before any work starts
  added: '2026-02-24'
  priority: P0
  type: Chore
  status: open
  issue: '#198'
---

**Research first**: Which .claude/rules/*.md files lack paths scoping? Which skills are candidates for disable-model-invocation? Can domain-router subagents with skills field replace always-on orchestrator skill loading? See .planning/research/context-optimization-research.md
