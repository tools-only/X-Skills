# Data Import and Export Reference

| Property | Value |
|----------|-------|
| **Name** | Data Import and Export Reference |
| **Repository** | [K-Dense-AI/claude-scientific-skills](https://raw.githubusercontent.com/K-Dense-AI/claude-scientific-skills/main/scientific-skills/matlab/references/data-import-export.md) (🔥 7.7k) |
| **Original Path** | `scientific-skills/matlab/references/data-import-export.md` |
| **Category** | content-creation |
| **Subcategory** | writing |
| **Tags** | content creation |
| **Created** | 2026-01-05 |
| **Updated** | 2026-01-05 |
| **File Hash** | `6d598e56aa30e143...` |

## Description

% With options
T = readtable('data.csv', 'ReadVariableNames', true);
T = readtable('data.csv', 'Delimiter', ',');
T = readtable('data.csv', 'NumHeaderLines', 2);
M = readmatrix('data.csv', 'Range', 'B2:D100');

**Tags:** `content creation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [K-Dense-AI/claude-scientific-skills](https://raw.githubusercontent.com/K-Dense-AI/claude-scientific-skills/main/scientific-skills/matlab/references/data-import-export.md)*
