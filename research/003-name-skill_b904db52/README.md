# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/openevidence-pack/skills/openevidence-core-workflow-b/SKILL.md) (⭐ 1.3k) |
| **Original Path** | `plugins/saas-packs/openevidence-pack/skills/openevidence-core-workflow-b/SKILL.md` |
| **Category** | research |
| **Subcategory** | academic |
| **Tags** | research |
| **Created** | 2026-01-26 |
| **Updated** | 2026-01-26 |
| **File Hash** | `b904db52586c40aa...` |

## Description

interface DeepConsultRequest {
  question: string;
  specialty: string;
  researchFocus?: 'treatment' | 'diagnosis' | 'prognosis' | 'epidemiology';
  timeframe?: {
    startYear?: number;
    endYear?: number;
  };
  preferredSources?: string[];
  excludeSources?: string[];
}

**Tags:** `research`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/openevidence-pack/skills/openevidence-core-workflow-b/SKILL.md)*
