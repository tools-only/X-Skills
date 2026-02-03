# web-app-service

| Property | Value |
|----------|-------|
| **Name** | web-app-service |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/kubernetes-specialist/references/networking.md) (⭐ 216) |
| **Original Path** | `skills/kubernetes-specialist/references/networking.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2025-12-15 |
| **Updated** | 2026-01-29 |
| **File Hash** | `12ff07c336439777...` |

## Description

yaml
apiVersion: v1
kind: Service
metadata:
  name: postgresheadless
  namespace: database
spec:
  clusterIP: None   Headless
  selector:
    app: postgres
  ports:
   name: postgres
    port: 5432
    targetPort: 5432

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/kubernetes-specialist/references/networking.md)*
