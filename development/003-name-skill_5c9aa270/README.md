# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [existential-birds/beagle](https://raw.githubusercontent.com/existential-birds/beagle/main/skills/prometheus-go-code-review/SKILL.md) (⭐ 17) |
| **Original Path** | `skills/prometheus-go-code-review/SKILL.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2025-12-21 |
| **Updated** | 2025-12-21 |
| **File Hash** | `5c9aa270d0f5209e...` |

## Description

go
// BAD  unique per user/request
counter := promauto.NewCounterVec(
    prometheus.CounterOpts{Name: "requests_total"},
    []string{"user_id", "path"},  // millions of series!
)
counter.WithLabelValues(userID, request.URL.Path).Inc()

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [existential-birds/beagle](https://raw.githubusercontent.com/existential-birds/beagle/main/skills/prometheus-go-code-review/SKILL.md)*
