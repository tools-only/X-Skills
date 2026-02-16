# Design

| Property | Value |
|----------|-------|
| **Name** | Design |
| **Repository** | [rocklambros/zerg](https://raw.githubusercontent.com/rocklambros/zerg/main/.gsd/specs/security-hardening/design.md) (⭐ 17) |
| **Original Path** | `.gsd/specs/security-hardening/design.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-05 |
| **Updated** | 2026-02-05 |
| **File Hash** | `52059e1480a4eba8...` |

## Description

2. Security Scan Flow:
    run_security_scan(path) called
    os.walk with followlinks=False iterates directory
    Each file path resolved and validated against scan boundary
    Symlinks outside boundary logged and skipped
    Results aggregated and returned

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [rocklambros/zerg](https://raw.githubusercontent.com/rocklambros/zerg/main/.gsd/specs/security-hardening/design.md)*
