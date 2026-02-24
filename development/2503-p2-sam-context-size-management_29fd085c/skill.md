---
name: 'SAM: Context Size Management'
description: Define explicit guidance for measuring and managing context size per agent. What's the target token budget? How to detect context pressure?
metadata:
  topic: sam-context-size-management
  source: Gap analysis of SAM framework
  added: '2026-02-01'
  priority: P2
  type: Feature
  status: open
---
**Suggested location**: [`stateless-software-engineering-framework.md`](https://github.com/bitflight-devops/stateless-agent-methodology/blob/main/stateless-software-engineering-framework.md) (section 2.1 or Appendix C)

**Research first**: How do agent frameworks measure context usage? What token counting approaches exist? How does Claude Code handle context limits internally?
