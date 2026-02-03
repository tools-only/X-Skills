# Semantically Dead Code

| Property | Value |
|----------|-------|
| **Name** | Semantically Dead Code |
| **Repository** | [alchemiststudiosDOTai/tunacode](https://raw.githubusercontent.com/alchemiststudiosDOTai/tunacode/master/.claude/qa/semantically-dead-code.md) (⭐ 112) |
| **Original Path** | `.claude/qa/semantically-dead-code.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-01-15 |
| **Updated** | 2026-01-15 |
| **File Hash** | `824ee2a7e8f91705...` |

## Description

glob.py had a use_gitignore parameter that did nothing. The function _load_gitignore_patterns() was called, populated a global _gitignore_patterns, but that global was never read. The actual filtering used DEFAULT_EXCLUDE_DIRS instead.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [alchemiststudiosDOTai/tunacode](https://raw.githubusercontent.com/alchemiststudiosDOTai/tunacode/master/.claude/qa/semantically-dead-code.md)*
