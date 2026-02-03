# app-blue

| Property | Value |
|----------|-------|
| **Name** | app-blue |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/devops-engineer/references/deployment-strategies.md) (⭐ 216) |
| **Original Path** | `skills/devops-engineer/references/deployment-strategies.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2025-12-14 |
| **Updated** | 2026-01-29 |
| **File Hash** | `d09fd909526bac88...` |

## Description

yaml
apiVersion: apps/v1
kind: Deployment
spec:
  strategy:
    type: RollingUpdate
    rollingUpdate:
      maxSurge: 25%         Max pods above desired
      maxUnavailable: 25%   Max pods unavailable

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/devops-engineer/references/deployment-strategies.md)*
