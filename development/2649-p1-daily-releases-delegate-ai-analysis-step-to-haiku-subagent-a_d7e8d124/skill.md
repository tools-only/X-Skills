---
name: 'daily-releases: delegate AI analysis step to Haiku subagent and fix duplicate draft releases'
description: "Step 2b in the daily-releases skill reads commits_detailed.txt, changes.diff, etc. into the orchestrator context and runs AI analysis inline. When orchestrator context is low or the data is large, the analysis produces zero file stats and generic 'Feature added by Claude' messages. Fix: delegate Step 2b to a Haiku subagent that receives file paths, reads them fresh, and writes analysis.json directly. Separately: publish_daily_release.py creates duplicate draft releases when a run is interrupted "
metadata:
  topic: daily-releases-delegate-ai-analysis-step-to-haiku-subagent-a
  source: 'Session context: bad release body from orchestrator low-context, duplicate draft from partial re-run'
  added: '2026-02-24'
  priority: P1
  type: Bug
  status: resolved
  issue: '#196'
---

## Story

As a **developer**, I want **Step 2b in the daily-releases skill reads commits_detailed** so that **backlog items are tracked in GitHub**.

## Description

Step 2b in the daily-releases skill reads commits_detailed.txt, changes.diff, etc. into the orchestrator context and runs AI analysis inline. When orchestrator context is low or the data is large, the analysis produces zero file stats and generic 'Feature added by Claude' messages. Fix: delegate Step 2b to a Haiku subagent that receives file paths, reads them fresh, and writes analysis.json directly. Separately: publish_daily_release.py creates duplicate draft releases when a run is interrupted mid-rename (release deleted, new one created) and then re-run with the same head_ref. The re-run finds no existing release and calls create_git_release a second time. Fix: add orphan-draft cleanup in publish_daily_release.py before creating a new release. Suggested locations: SKILL.md step 2b, publish_daily_release.py.

## Acceptance Criteria

- [ ] Work matches description
- [ ] Plan or implementation complete

## Context

- **Source**: Session context: bad release body from orchestrator low-context, duplicate draft from partial re-run
- **Priority**: P1
- **Added**: 2026-02-24
- **Research questions**: None