# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/clinical-databases/hla-typing/SKILL.md) (⭐ 191) |
| **Original Path** | `clinical-databases/hla-typing/SKILL.md` |
| **Category** | data-analysis |
| **Subcategory** | processing |
| **Tags** | data analysis |
| **Created** | 2026-01-29 |
| **Updated** | 2026-01-29 |
| **File Hash** | `4adc4c4c672f5ea4...` |

## Description

bash
 Extract HLA reads from BAM
samtools view h input.bam chr6:2800000034000000 | \
    samtools fastq 1 hla_R1.fq 2 hla_R2.fq

**Tags:** `data analysis`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/clinical-databases/hla-typing/SKILL.md)*
