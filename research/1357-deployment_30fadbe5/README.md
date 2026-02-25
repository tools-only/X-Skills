# Deployment

| Property | Value |
|----------|-------|
| **Name** | Deployment |
| **Repository** | [tomascupr/sandstorm](https://raw.githubusercontent.com/tomascupr/sandstorm/main/docs/deployment.md) (⭐ 399) |
| **Original Path** | `docs/deployment.md` |
| **Category** | research |
| **Subcategory** | data-gathering |
| **Tags** | research |
| **Created** | 2026-02-15 |
| **Updated** | 2026-02-15 |
| **File Hash** | `30fadbe5e49dcc8d...` |

## Description

Sandstorm is a stateless FastAPI app. Each request creates an independent E2B sandbox, runs the agent, and tears it down. No shared state, no sticky sessions, no coordination between requests. This means deploying for concurrent agent runs is trivial  just add workers.

**Tags:** `research`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [tomascupr/sandstorm](https://raw.githubusercontent.com/tomascupr/sandstorm/main/docs/deployment.md)*
