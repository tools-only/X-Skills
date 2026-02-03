# databases.mycompany.io

| Property | Value |
|----------|-------|
| **Name** | databases.mycompany.io |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/kubernetes-specialist/references/custom-operators.md) (⭐ 216) |
| **Original Path** | `skills/kubernetes-specialist/references/custom-operators.md` |
| **Category** | development |
| **Subcategory** | devops |
| **Tags** | development |
| **Created** | 2025-12-21 |
| **Updated** | 2026-01-29 |
| **File Hash** | `381e5cdc26f12cb0...` |

## Description

yaml
apiVersion: mycompany.io/v1
kind: Database
metadata:
  name: ordersdb
  namespace: production
spec:
  engine: postgres
  version: "15.4"
  storage: 100Gi
  replicas: 3

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/kubernetes-specialist/references/custom-operators.md)*
