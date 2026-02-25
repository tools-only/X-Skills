---
name: Fix pre-existing linting errors in gitlab_context.py
description: 'plugins/gitlab-skill/skills/gitlab-skill/scripts/gitlab_context.py has 3 linting errors: PLC0415 x2 (deferred imports inside functions) and S607 (partial executable path in subprocess). Fix: move imports to top-level, use shutil.which() for glab path resolution. Same pattern was applied in fetch_gitlab_mr.py migration.'
metadata:
  topic: fix-pre-existing-linting-errors-in-gitlabcontextpy
  source: Discovered during fetch_gitlab_mr.py subprocess-to-python-gitlab migration
  added: '2026-02-24'
  priority: P2
  type: Tech Debt
  status: open
---

