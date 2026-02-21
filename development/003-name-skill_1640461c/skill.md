---
name: lint
description: Run the plugin validator on a skill, agent, or plugin directory rts token complexity, broken links, frontmatter issues, and structural problems. Use when checking skill quality, validating before commit, or diagnosing validator warnings. Pass the path as an argument.
argument-hint: <path-to-skill-or-plugin>
user-invocable: true
---
!`${CLAUDE_PLUGIN_ROOT}/scripts/plugin_validator.py $ARGUMENTS`

@${CLAUDE_PLUGIN_ROOT}/references/ERROR_CODES.md
