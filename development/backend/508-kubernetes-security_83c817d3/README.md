# pod-reader

| Property | Value |
|----------|-------|
| **Name** | pod-reader |
| **Repository** | [mrgoonie/claudekit-skills](https://raw.githubusercontent.com/mrgoonie/claudekit-skills/main/.claude/skills/devops/references/kubernetes-security.md) (⭐ 1.5k) |
| **Original Path** | `.claude/skills/devops/references/kubernetes-security.md` |
| **Category** | daily-assistant |
| **Subcategory** | tasks |
| **Tags** | daily assistant |
| **Created** | 2026-01-17 |
| **Updated** | 2026-01-17 |
| **File Hash** | `83c817d3dae55f4e...` |

## Description

yaml
spec:
  securityContext:
    runAsNonRoot: true
    runAsUser: 1000
    seccompProfile:
      type: RuntimeDefault
  containers:
   name: app
    securityContext:
      allowPrivilegeEscalation: false
      readOnlyRootFilesystem: true
      capabilities:
        drop: ["ALL"]

**Tags:** `daily assistant`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [mrgoonie/claudekit-skills](https://raw.githubusercontent.com/mrgoonie/claudekit-skills/main/.claude/skills/devops/references/kubernetes-security.md)*
