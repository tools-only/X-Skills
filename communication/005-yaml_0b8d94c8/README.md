# Alerting Rules

| Property | Value |
|----------|-------|
| **Name** | Alerting Rules |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/monitoring-expert/references/alerting-rules.md) (⭐ 216) |
| **Original Path** | `skills/monitoring-expert/references/alerting-rules.md` |
| **Category** | communication |
| **Subcategory** | messaging |
| **Tags** | communication |
| **Created** | 2025-12-14 |
| **Updated** | 2026-01-29 |
| **File Hash** | `0b8d94c81c09ce12...` |

## Description

alert: ServiceDown
        expr: up == 0
        for: 1m
        labels:
          severity: critical
        annotations:
          summary: Service {{ $labels.instance }} is down

**Tags:** `communication`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/monitoring-expert/references/alerting-rules.md)*
