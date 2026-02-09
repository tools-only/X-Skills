# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/linear-pack/skills/linear-debug-bundle/SKILL.md) (⭐ 1.3k) |
| **Original Path** | `plugins/saas-packs/linear-pack/skills/linear-debug-bundle/SKILL.md` |
| **Category** | development |
| **Subcategory** | tools |
| **Tags** | development |
| **Created** | 2026-01-09 |
| **Updated** | 2026-01-09 |
| **File Hash** | `2ada47fa47135283...` |

## Description

interface DebugOptions {
  logRequests?: boolean;
  logResponses?: boolean;
  logErrors?: boolean;
  onRequest?: (query: string, variables: unknown) => void;
  onResponse?: (data: unknown, duration: number) => void;
  onError?: (error: Error, duration: number) => void;
}

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/linear-pack/skills/linear-debug-bundle/SKILL.md)*
