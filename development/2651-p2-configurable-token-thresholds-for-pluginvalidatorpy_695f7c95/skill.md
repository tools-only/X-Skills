---
name: Configurable Token Thresholds for plugin_validator.py
description: '`plugin_validator.py` token complexity thresholds (4400 warn, 8800 error) are hardcoded constants. Allow per-project or per-plugin configuration via a config file (e.g., `.claude-plugin/validator.json` or a `[tool.plugin-validator]` section in `pyproject.toml`) so different projects can set their own limits.'
metadata:
  topic: configurable-token-thresholds-for-pluginvalidatorpy
  source: Session 2026-02-18, token threshold analysis
  added: '2026-02-18'
  priority: P2
  type: Feature
  status: open
  issue: '#119'
---

**Suggested location**: `plugins/plugin-creator/scripts/plugin_validator.py`
