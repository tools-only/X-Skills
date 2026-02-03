# {{ include

| Property | Value |
|----------|-------|
| **Name** | {{ include |
| **Repository** | [mrgoonie/claudekit-skills](https://raw.githubusercontent.com/mrgoonie/claudekit-skills/main/.claude/skills/devops/references/kubernetes-helm-advanced.md) (⭐ 1.5k) |
| **Original Path** | `.claude/skills/devops/references/kubernetes-helm-advanced.md` |
| **Category** | data-analysis |
| **Subcategory** | visualization |
| **Tags** | data analysis |
| **Created** | 2026-01-17 |
| **Updated** | 2026-01-17 |
| **File Hash** | `3fbb59829fbdc9cb...` |

## Description

{{ define "mychart.labels" }}
app.kubernetes.io/name: {{ .Chart.Name }}
app.kubernetes.io/instance: {{ .Release.Name }}
{{ end }}

**Tags:** `data analysis`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [mrgoonie/claudekit-skills](https://raw.githubusercontent.com/mrgoonie/claudekit-skills/main/.claude/skills/devops/references/kubernetes-helm-advanced.md)*
