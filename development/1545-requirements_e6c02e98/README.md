# Requirements

| Property | Value |
|----------|-------|
| **Name** | Requirements |
| **Repository** | [rocklambros/zerg](https://raw.githubusercontent.com/rocklambros/zerg/main/.gsd/specs/fix-codeql-scanning-alerts/requirements.md) (⭐ 17) |
| **Original Path** | `.gsd/specs/fix-codeql-scanning-alerts/requirements.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-15 |
| **Updated** | 2026-02-15 |
| **File Hash** | `e6c02e98fec87c57...` |

## Description

Additionally, the /z:plan command has a recurring bug where it fails to generate its primary output (requirements.md) due to an overbroad WORKFLOW BOUNDARY guard that prohibits the Write tool entirely — including for spec files which are the command's intended output.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [rocklambros/zerg](https://raw.githubusercontent.com/rocklambros/zerg/main/.gsd/specs/fix-codeql-scanning-alerts/requirements.md)*
