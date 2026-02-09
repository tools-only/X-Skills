# Architecture

| Property | Value |
|----------|-------|
| **Name** | Architecture |
| **Repository** | [julianobarbosa/claude-code-skills](https://raw.githubusercontent.com/julianobarbosa/claude-code-skills/main/skills/keyvault-csi-driver-skill/references/architecture.md) (⭐ 14) |
| **Original Path** | `skills/keyvault-csi-driver-skill/references/architecture.md` |
| **Category** | development |
| **Subcategory** | devops |
| **Tags** | development |
| **Created** | 2025-12-05 |
| **Updated** | 2025-12-05 |
| **File Hash** | `398d56e57d9d0842...` |

## Description

yaml
apiVersion: secretsstore.csi.xk8s.io/v1
kind: SecretProviderClass
metadata:
  name: example              Referenced by pod volume
  namespace: default         Must match pod namespace
spec:
  provider: azure            Provider type

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [julianobarbosa/claude-code-skills](https://raw.githubusercontent.com/julianobarbosa/claude-code-skills/main/skills/keyvault-csi-driver-skill/references/architecture.md)*
