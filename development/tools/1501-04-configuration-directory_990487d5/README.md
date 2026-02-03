# 04 Configuration Directory

| Property | Value |
|----------|-------|
| **Name** | 04 Configuration Directory |
| **Repository** | [alchemiststudiosDOTai/tunacode](https://raw.githubusercontent.com/alchemiststudiosDOTai/tunacode/master/docs/codebase-map/structure/04-configuration-directory.md) (⭐ 112) |
| **Original Path** | `docs/codebase-map/structure/04-configuration-directory.md` |
| **Category** | development |
| **Subcategory** | tools |
| **Tags** | development |
| **Created** | 2026-01-05 |
| **Updated** | 2026-01-25 |
| **File Hash** | `990487d52f7b5e63...` |

## Description

Base URL Resolution (in agent_config.py):
1. Perprovider config: settings.providers.<provider>.base_url
2. Registry default: models_registry.json → provider api field
3. OpenAI escape hatch: OPENAI_BASE_URL env var (nonAnthropic only)

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [alchemiststudiosDOTai/tunacode](https://raw.githubusercontent.com/alchemiststudiosDOTai/tunacode/master/docs/codebase-map/structure/04-configuration-directory.md)*
