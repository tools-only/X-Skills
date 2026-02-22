# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/time-aware-dependency-cve-scanner/SKILL.md) (⭐ 10) |
| **Original Path** | `skills/time-aware-dependency-cve-scanner/SKILL.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-20 |
| **Updated** | 2026-02-20 |
| **File Hash** | `635bd5c10da006f7...` |

## Description

Scan repositories for newly disclosed CVEs in dependencies after a specific cutoff date. Takes a repository path, cutoff date (YYYY-MM-DD), and optional parameters for transitive dependencies. Parses dependency manifests (package.json, pom.xml, requirements.txt, go.mod, Cargo.toml) and lockfiles to extract exact versions. Queries vulnerability databases (OSV.dev, NVD, GitHub Advisory) to identify CVEs disclosed strictly after the cutoff date. Distinguishes between newly disclosed CVEs and previously known CVEs. Use when: (1) Performing security audits to find new vulnerabilities since last review, (2) Checking if new CVEs affect a historical codebase version, (3) Generating compliance reports showing vulnerability status at specific dates, (4) Tracking security posture changes over time. Supports npm, Maven, pip, Go modules, Cargo, and other major ecosystems.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/time-aware-dependency-cve-scanner/SKILL.md)*
