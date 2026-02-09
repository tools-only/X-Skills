# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/langchain-pack/skills/langchain-performance-tuning/SKILL.md) (⭐ 1.3k) |
| **Original Path** | `plugins/saas-packs/langchain-pack/skills/langchain-performance-tuning/SKILL.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-01-09 |
| **Updated** | 2026-01-09 |
| **File Hash** | `3ab4d693ac179c29...` |

## Description

def benchmark(func: Callable, iterations: int = 10):
    """Benchmark a function's performance."""
    times = []
    for _ in range(iterations):
        start = time.perf_counter()
        func()
        elapsed = time.perf_counter()  start
        times.append(elapsed)

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/langchain-pack/skills/langchain-performance-tuning/SKILL.md)*
