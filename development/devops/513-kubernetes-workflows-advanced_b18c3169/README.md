# Build and Deploy

| Property | Value |
|----------|-------|
| **Name** | Build and Deploy |
| **Repository** | [mrgoonie/claudekit-skills](https://raw.githubusercontent.com/mrgoonie/claudekit-skills/main/.claude/skills/devops/references/kubernetes-workflows-advanced.md) (⭐ 1.5k) |
| **Original Path** | `.claude/skills/devops/references/kubernetes-workflows-advanced.md` |
| **Category** | development |
| **Subcategory** | devops |
| **Tags** | development |
| **Created** | 2026-01-17 |
| **Updated** | 2026-01-17 |
| **File Hash** | `b18c3169301374f8...` |

## Description

jobs:
  build:
    runson: ubuntulatest
    steps:
     uses: actions/checkout@v3
     run: docker build . t $REGISTRY/$IMAGE:${{ github.sha }}
     run: docker push $REGISTRY/$IMAGE:${{ github.sha }}

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [mrgoonie/claudekit-skills](https://raw.githubusercontent.com/mrgoonie/claudekit-skills/main/.claude/skills/devops/references/kubernetes-workflows-advanced.md)*
