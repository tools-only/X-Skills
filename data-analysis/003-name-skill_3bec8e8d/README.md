# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/copy-number/cnv-annotation/SKILL.md) (⭐ 191) |
| **Original Path** | `copy-number/cnv-annotation/SKILL.md` |
| **Category** | data-analysis |
| **Subcategory** | processing |
| **Tags** | data analysis |
| **Created** | 2026-01-17 |
| **Updated** | 2026-01-23 |
| **File Hash** | `3bec8e8de86f2c07...` |

## Description

bash
 Convert CNV segments to BED
awk 'NR>1 {print $1"\t"$2"\t"$3"\t"$5"\t"$6}' sample.cns > sample.cnv.bed

**Tags:** `data analysis`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/copy-number/cnv-annotation/SKILL.md)*
