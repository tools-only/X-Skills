---
name: 'P1: plugin-validator pre-commit output is too noisy'
description: Running `prek --all-files` produces pages of plugin-validator output — every file in the repo gets validated, and pre-existing warnings/errors (broken links, SK006 complexity, missing trigger phrases) for files unrelated to the current change drown out actionable information. The validator needs a review of its pre-commit integration to reduce noise.
metadata:
  topic: p1-plugin-validator-pre-commit-output-is-too-noisy
  source: '`prek --all-files` run during session 2026-02-15'
  added: '2026-02-15'
  priority: completed
  type: Feature
  status: done
  issue: '#130'
  plan: ''
---