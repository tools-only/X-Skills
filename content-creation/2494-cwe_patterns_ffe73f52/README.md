# Cwe Patterns

| Property | Value |
|----------|-------|
| **Name** | Cwe Patterns |
| **Repository** | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/static-vulnerability-detector/references/cwe_patterns.md) (⭐ 10) |
| **Original Path** | `skills/static-vulnerability-detector/references/cwe_patterns.md` |
| **Category** | content-creation |
| **Subcategory** | writing |
| **Tags** | content creation |
| **Created** | 2026-02-20 |
| **Updated** | 2026-02-20 |
| **File Hash** | `ffe73f5229c6397a...` |

## Description

// Vulnerable: Offbyone
for(int i = 0; i <= 10; i++) {  // VULNERABLE: should be i < 10
    buffer[i] = data[i];
}

Remediation:
 Use safe functions: strncpy, snprintf
 Check bounds before writing
 Use modern C++ containers (std::string, std::vector)

**Tags:** `content creation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/static-vulnerability-detector/references/cwe_patterns.md)*
