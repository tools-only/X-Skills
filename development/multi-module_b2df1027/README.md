# Multi-Module Migration Patterns

| Property | Value |
|----------|-------|
| **Name** | Multi-Module Migration Patterns |
| **Repository** | [adityamparikh/maven-to-gradle](https://raw.githubusercontent.com/adityamparikh/maven-to-gradle/main/references/multi-module.md) (⭐ 10) |
| **Original Path** | `references/multi-module.md` |
| **Category** | daily-assistant |
| **Subcategory** | tasks |
| **Tags** | daily assistant |
| **Created** | 2026-02-02 |
| **Updated** | 2026-02-02 |
| **File Hash** | `b2df1027e015ad6a...` |

## Description

Maven multimodule:

parent/
├── pom.xml            (packaging: pom, <modules>)
├── moduleapi/
│   └── pom.xml        (<parent> → parent)
├── modulecore/
│   └── pom.xml
└── moduleweb/
    └── pom.xml

**Tags:** `daily assistant`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [adityamparikh/maven-to-gradle](https://raw.githubusercontent.com/adityamparikh/maven-to-gradle/main/references/multi-module.md)*
