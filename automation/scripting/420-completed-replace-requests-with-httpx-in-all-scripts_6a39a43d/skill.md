---
name: Replace requests with httpx in all scripts
description: Migrated `validate_glfm.py` from `requests` to `httpx`. Added TID251 ruff ban rule for `requests` imports. Removed `types-requests` from dev dependencies. `sync_gitlab_docs.py` was already using httpx.
metadata:
  topic: replace-requests-with-httpx-in-all-scripts
  source: CI/pre-commit inconsistency discovery (2026-02-05)
  added: '2026-02-23'
  priority: completed
  type: Feature
  status: done
---
