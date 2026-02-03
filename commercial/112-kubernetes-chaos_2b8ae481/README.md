# Kubernetes Chaos

| Property | Value |
|----------|-------|
| **Name** | Kubernetes Chaos |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/chaos-engineer/references/kubernetes-chaos.md) (⭐ 216) |
| **Original Path** | `skills/chaos-engineer/references/kubernetes-chaos.md` |
| **Category** | commercial |
| **Subcategory** | ecommerce |
| **Tags** | commercial |
| **Created** | 2025-12-15 |
| **Updated** | 2026-01-29 |
| **File Hash** | `2b8ae4817420a156...` |

## Description

yaml
apiVersion: litmuschaos.io/v1alpha1
kind: ChaosEngine
metadata:
  name: nginxchaos
  namespace: default
spec:
   Application information
  appinfo:
    appns: 'default'
    applabel: 'app=nginx'
    appkind: 'deployment'

**Tags:** `commercial`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/chaos-engineer/references/kubernetes-chaos.md)*
