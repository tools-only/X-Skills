# secret-reader

| Property | Value |
|----------|-------|
| **Name** | secret-reader |
| **Repository** | [mrgoonie/claudekit-skills](https://raw.githubusercontent.com/mrgoonie/claudekit-skills/main/.claude/skills/devops/references/kubernetes-security-advanced.md) (⭐ 1.5k) |
| **Original Path** | `.claude/skills/devops/references/kubernetes-security-advanced.md` |
| **Category** | daily-assistant |
| **Subcategory** | tasks |
| **Tags** | daily assistant |
| **Created** | 2026-01-17 |
| **Updated** | 2026-01-17 |
| **File Hash** | `26e8cbd516658233...` |

## Description

yaml
apiVersion: rbac.authorization.k8s.io/v1
kind: ClusterRole
metadata:
  name: secretreader
rules:
 apiGroups: [""]
  resources: ["secrets"]
  verbs: ["get"]
  resourceNames: ["appcredentials"]   Restrict to specific

**Tags:** `daily assistant`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [mrgoonie/claudekit-skills](https://raw.githubusercontent.com/mrgoonie/claudekit-skills/main/.claude/skills/devops/references/kubernetes-security-advanced.md)*
