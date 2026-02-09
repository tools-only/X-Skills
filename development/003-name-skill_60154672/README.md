# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/openevidence-pack/skills/openevidence-core-workflow-a/SKILL.md) (⭐ 1.3k) |
| **Original Path** | `plugins/saas-packs/openevidence-pack/skills/openevidence-core-workflow-a/SKILL.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-01-26 |
| **Updated** | 2026-01-26 |
| **File Hash** | `6015467227f1f838...` |

## Description

interface ClinicalQueryRequest {
  question: string;
  specialty: string;
  urgency: 'stat' | 'urgent' | 'routine';
  patientContext?: {
    age?: number;
    sex?: 'male' | 'female';
    conditions?: string[];
    medications?: string[];
  };
}

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/openevidence-pack/skills/openevidence-core-workflow-a/SKILL.md)*
