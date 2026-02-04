# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/atac-seq/atac-peak-calling/SKILL.md) (⭐ 191) |
| **Original Path** | `atac-seq/atac-peak-calling/SKILL.md` |
| **Category** | data-analysis |
| **Subcategory** | visualization |
| **Tags** | data analysis |
| **Created** | 2026-01-17 |
| **Updated** | 2026-01-20 |
| **File Hash** | `63aae66778daa524...` |

## Description

bash
 Standard ATACseq peak calling
macs3 callpeak \
    t sample.bam \
    f BAMPE \
    g hs \
    n sample \
    outdir peaks/ \
    q 0.05 \
    nomodel \
    shift 75 \
    extsize 150 \
    keepdup all \
    B

**Tags:** `data analysis`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/atac-seq/atac-peak-calling/SKILL.md)*
