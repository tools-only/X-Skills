# Spec

| Property | Value |
|----------|-------|
| **Name** | Spec |
| **Repository** | [lfnovo/esperanto](https://raw.githubusercontent.com/lfnovo/esperanto/main/specs/openai-compatible-env-vars/spec.md) (⭐ 142) |
| **Original Path** | `specs/openai-compatible-env-vars/spec.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2025-10-19 |
| **Updated** | 2025-10-19 |
| **File Hash** | `27767d0d3aa1420a...` |

## Description

We have several openai compatible providers:

@src/esperanto/providers/stt/openai_compatible.py
@src/esperanto/providers/tts/openai_compatible.py
@src/esperanto/providers/llm/openai_compatible.py
@src/esperanto/providers/embedding/openai_compatible.py

I noticed an issue with the way we configured this. All of them are using the same BASE_URL and API_KEY environment variables. So, if the user has different services for each of them, it won

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [lfnovo/esperanto](https://raw.githubusercontent.com/lfnovo/esperanto/main/specs/openai-compatible-env-vars/spec.md)*
