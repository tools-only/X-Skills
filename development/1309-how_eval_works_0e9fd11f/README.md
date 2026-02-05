# How Eval Works

| Property | Value |
|----------|-------|
| **Name** | How Eval Works |
| **Repository** | [cisco-ai-defense/skill-scanner](https://raw.githubusercontent.com/cisco-ai-defense/skill-scanner/main/evals/HOW_EVAL_WORKS.md) (⭐ 324) |
| **Original Path** | `evals/HOW_EVAL_WORKS.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-01-29 |
| **Updated** | 2026-01-29 |
| **File Hash** | `0e9fd11fde1899a4...` |

## Description

1. Static Analyzer Only (includes YARA):
   bash
   python evals/eval_runner.py testskillsdir evals/skills
   
    Uses StaticAnalyzer which includes YARA rule matching
    Fast, no API calls
    Tests patternbased detection

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [cisco-ai-defense/skill-scanner](https://raw.githubusercontent.com/cisco-ai-defense/skill-scanner/main/evals/HOW_EVAL_WORKS.md)*
