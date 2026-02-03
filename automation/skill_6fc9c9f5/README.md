# {{ skill_name }}

| Property | Value |
|----------|-------|
| **Name** | {{ skill_name }} |
| **Repository** | [PrefectHQ/colin](https://raw.githubusercontent.com/PrefectHQ/colin/main/src/colin/blueprints/mcp-guide/models/SKILL.md) (⭐ 93) |
| **Original Path** | `src/colin/blueprints/mcp-guide/models/SKILL.md` |
| **Category** | content-creation |
| **Subcategory** | writing |
| **Tags** | content creation |
| **Created** | 2026-01-23 |
| **Updated** | 2026-01-23 |
| **File Hash** | `6fc9c9f5795c125a...` |

## Description

Server name: {{ server.name }}
Server description: {{ server.instructions or "No description provided" }}
{% endllm %}
{% if vars.generate_scripts | default(true) %}
allowedtools: Read, Bash(python:)
{% endif %}

**Tags:** `content creation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [PrefectHQ/colin](https://raw.githubusercontent.com/PrefectHQ/colin/main/src/colin/blueprints/mcp-guide/models/SKILL.md)*
