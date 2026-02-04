# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/metabolomics/statistical-analysis/SKILL.md) (⭐ 191) |
| **Original Path** | `metabolomics/statistical-analysis/SKILL.md` |
| **Category** | data-analysis |
| **Subcategory** | visualization |
| **Tags** | data analysis |
| **Created** | 2026-01-20 |
| **Updated** | 2026-01-24 |
| **File Hash** | `f12a64e3bb169097...` |

## Description

r
 Calculate fold change between groups
calculate_fc < function(data, groups) {
    group_means < aggregate(data, by = list(groups), FUN = mean, na.rm = TRUE)
    rownames(group_means) < group_means$Group.1
    group_means < group_means[, 1]

**Tags:** `data analysis`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/metabolomics/statistical-analysis/SKILL.md)*
