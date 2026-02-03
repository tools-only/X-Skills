# app-config

| Property | Value |
|----------|-------|
| **Name** | app-config |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/kubernetes-specialist/references/configuration.md) (⭐ 216) |
| **Original Path** | `skills/kubernetes-specialist/references/configuration.md` |
| **Category** | development |
| **Subcategory** | devops |
| **Tags** | development |
| **Created** | 2025-12-15 |
| **Updated** | 2026-01-29 |
| **File Hash** | `0677fa7f70d2001c...` |

## Description

yaml
apiVersion: v1
kind: ConfigMap
metadata:
  name: appconfig
  namespace: production
data:
   Simple keyvalue pairs
  database.host: "postgresservice.database.svc.cluster.local"
  database.port: "5432"
  database.name: "appdb"

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/kubernetes-specialist/references/configuration.md)*
