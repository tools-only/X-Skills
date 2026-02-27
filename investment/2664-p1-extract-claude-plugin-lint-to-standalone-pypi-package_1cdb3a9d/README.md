# P1 Extract Claude Plugin Lint To Standalone Pypi Package

| Property | Value |
|----------|-------|
| **Name** | P1 Extract Claude Plugin Lint To Standalone Pypi Package |
| **Repository** | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/backlog/p1-extract-claude-plugin-lint-to-standalone-pypi-package.md) (⭐ 20) |
| **Original Path** | `.claude/backlog/p1-extract-claude-plugin-lint-to-standalone-pypi-package.md` |
| **Category** | investment |
| **Subcategory** | trading |
| **Tags** | investment |
| **Created** | 2026-02-23 |
| **Updated** | 2026-02-27 |
| **File Hash** | `1cdb3a9dab95ee1f...` |

## Description

Extract and enhance `validate_frontmatter.py` into a standalone open-source project. First dedicated linter for Claude Code plugin frontmatter (SKILL.md, agents/*.md, commands/*.md). Official `claude plugin validate` only checks plugin.json structure.\n**Features to include**:\n- YAML frontmatter schema validation with Pydantic models\n- Auto-fix capabilities (arrays → comma-separated, multiline → single-line)\n- Token-based complexity metrics (tiktoken) instead of line counts\n- Cross-reference validation (agent references non-existent skill)\n- Marketplace readiness scoring\n- Pre-commit hook integration\n- CLI with `--fix` and `--report` modes\n**Current source**: `plugins/plugin-creator/scripts/validate_frontmatter.py`\n**Suggested repo name**: `claude-plugin-lint` or `cc-plugin-validator`\n\n---

**Tags:** `investment`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/backlog/p1-extract-claude-plugin-lint-to-standalone-pypi-package.md)*
