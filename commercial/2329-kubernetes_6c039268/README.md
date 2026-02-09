# Kubernetes

| Property | Value |
|----------|-------|
| **Name** | Kubernetes |
| **Repository** | [julianobarbosa/claude-code-skills](https://raw.githubusercontent.com/julianobarbosa/claude-code-skills/main/skills/managing-infra-skill/KUBERNETES.md) (⭐ 14) |
| **Original Path** | `skills/managing-infra-skill/KUBERNETES.md` |
| **Category** | commercial |
| **Subcategory** | ecommerce |
| **Tags** | commercial |
| **Created** | 2025-12-16 |
| **Updated** | 2025-12-16 |
| **File Hash** | `6c0392687f8aa335...` |

## Description

yaml
apiVersion: v1
kind: Service
metadata:
  name: app
spec:
  selector:
    app.kubernetes.io/name: app
  ports:
     port: 80
      targetPort: 8080
  type: ClusterIP

**Tags:** `commercial`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [julianobarbosa/claude-code-skills](https://raw.githubusercontent.com/julianobarbosa/claude-code-skills/main/skills/managing-infra-skill/KUBERNETES.md)*
