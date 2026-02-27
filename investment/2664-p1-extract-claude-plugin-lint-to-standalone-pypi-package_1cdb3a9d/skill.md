---
name: Extract claude-plugin-lint to standalone PyPI package
description: "Extract and enhance `validate_frontmatter.py` into a standalone open-source project. First dedicated linter for Claude Code plugin frontmatter (SKILL.md, agents/*.md, commands/*.md). Official `claude plugin validate` only checks plugin.json structure.\n**Features to include**:\n- YAML frontmatter schema validation with Pydantic models\n- Auto-fix capabilities (arrays → comma-separated, multiline → single-line)\n- Token-based complexity metrics (tiktoken) instead of line counts\n- Cross-reference validation (agent references non-existent skill)\n- Marketplace readiness scoring\n- Pre-commit hook integration\n- CLI with `--fix` and `--report` modes\n**Current source**: `plugins/plugin-creator/scripts/validate_frontmatter.py`\n**Suggested repo name**: `claude-plugin-lint` or `cc-plugin-validator`\n\n---"
metadata:
  topic: extract-claude-plugin-lint-to-standalone-pypi-package
  source: Gap analysis - no existing Claude Code plugin linters exist
  added: '2026-02-01'
  priority: completed
  type: Feature
  status: done
  issue: '#93'
  plan: ''
---