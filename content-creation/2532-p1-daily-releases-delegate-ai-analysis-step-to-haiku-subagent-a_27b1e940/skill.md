---
name: 'daily-releases: delegate AI analysis step to Haiku subagent and fix duplicate draft releases'
description: "Step 2b in the daily-releases skill reads commits_detailed.txt, changes.diff, etc. into the orchestrator context and runs AI analysis inline. When orchestrator context is low or the data is large, the analysis produces zero file stats and generic 'Feature added by Claude' messages. Fix: delegate Step 2b to a Haiku subagent that receives file paths, reads them fresh, and writes analysis.json directly. Separately: publish_daily_release.py creates duplicate draft releases when a run is interrupted "
metadata:
  topic: daily-releases-delegate-ai-analysis-step-to-haiku-subagent-a
  source: 'Session context: bad release body from orchestrator low-context, duplicate draft from partial re-run'
  added: '2026-02-24'
  priority: P1
  type: Bug
  status: open
  issue: '#196'
---