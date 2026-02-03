# app

| Property | Value |
|----------|-------|
| **Name** | app |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/devops-engineer/references/kubernetes.md) (⭐ 216) |
| **Original Path** | `skills/devops-engineer/references/kubernetes.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2025-12-14 |
| **Updated** | 2026-01-29 |
| **File Hash** | `7da17b4e40c05d8c...` |

## Description

yaml
apiVersion: v1
kind: ConfigMap
metadata:
  name: appconfig
data:
  LOG_LEVEL: "info"
  API_TIMEOUT: "30s"

apiVersion: v1
kind: Secret
metadata:
  name: appsecrets
type: Opaque
stringData:
  databaseurl: "postgres://user:pass@host:5432/db"

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/devops-engineer/references/kubernetes.md)*
