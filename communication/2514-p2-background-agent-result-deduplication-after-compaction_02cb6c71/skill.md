---
name: Background agent result deduplication after compaction
description: Background agents launched in session N may complete after context compaction or session restart. The system delivers their results as `<task-notification>` messages, but there is no mechanism to mark results as already consumed. During the kaizen build, 3 review agents (plugin-validator, 2x skill-reviewer) completed during Phase 6 and their findings were applied in commit `0d61480`. After compaction, all 3 re-delivered their notifications in session 3, requiring manual evaluation each time ('wa
metadata:
  topic: background-agent-result-deduplication-after-compaction
  source: agentskill-kaizen plugin build sessions 2-3 (2026-02-18)
  added: '2026-02-18'
  priority: P2
  type: Feature
  status: open
---
