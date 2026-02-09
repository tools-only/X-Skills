# Model Information Api

| Property | Value |
|----------|-------|
| **Name** | Model Information Api |
| **Repository** | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/openrouter-pack/skills/openrouter-model-catalog/references/model-information-api.md) (⭐ 1.3k) |
| **Original Path** | `plugins/saas-packs/openrouter-pack/skills/openrouter-model-catalog/references/model-information-api.md` |
| **Category** | automation |
| **Subcategory** | scripting |
| **Tags** | automation |
| **Created** | 2026-01-06 |
| **Updated** | 2026-01-06 |
| **File Hash** | `46c2efb88920da44...` |

## Description

model = get_model_info("anthropic/claude3.5sonnet")
print(f"Context: {model['context_length']}")
print(f"Prompt cost: ${model['pricing']['prompt']}/token")
print(f"Completion cost: ${model['pricing']['completion']}/token")

**Tags:** `automation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/openrouter-pack/skills/openrouter-model-catalog/references/model-information-api.md)*
