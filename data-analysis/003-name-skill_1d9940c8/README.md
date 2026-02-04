# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/metagenomics/abundance-estimation/SKILL.md) (⭐ 191) |
| **Original Path** | `metagenomics/abundance-estimation/SKILL.md` |
| **Category** | data-analysis |
| **Subcategory** | processing |
| **Tags** | data analysis |
| **Created** | 2026-01-16 |
| **Updated** | 2026-01-20 |
| **File Hash** | `1d9940c88d1d62e7...` |

## Description

bash
 Run Bracken on Kraken2 report
bracken d /path/to/kraken2_db \
    i kraken_report.txt \
    o bracken_output.txt \
    r 150 \                        Read length (100, 150, 200, 250, 300)
    l S                            Taxonomic level

**Tags:** `data analysis`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/metagenomics/abundance-estimation/SKILL.md)*
