# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/flow-cytometry/fcs-handling/SKILL.md) (⭐ 191) |
| **Original Path** | `flow-cytometry/fcs-handling/SKILL.md` |
| **Category** | data-analysis |
| **Subcategory** | visualization |
| **Tags** | data analysis |
| **Created** | 2026-01-20 |
| **Updated** | 2026-01-24 |
| **File Hash** | `ce59515c56e426da...` |

## Description

r
 Read multiple files into flowSet
files < list.files('data', pattern = '\\.fcs$', full.names = TRUE)
fs < read.flowSet(files, transformation = FALSE, truncate_max_range = FALSE)

**Tags:** `data analysis`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/flow-cytometry/fcs-handling/SKILL.md)*
