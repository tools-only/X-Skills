# 0001 Plugins

| Property | Value |
|----------|-------|
| **Name** | 0001 Plugins |
| **Repository** | [strands-agents/docs](https://raw.githubusercontent.com/strands-agents/docs/main/designs/0001-plugins.md) (⭐ 167) |
| **Original Path** | `designs/0001-plugins.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-13 |
| **Updated** | 2026-02-13 |
| **File Hash** | `a1d2acd875bad1c0...` |

## Description

python
agent = Agent(
    model=BedrockModel(),             Model provider
    tools=[calculator, search],       Tool functions
    system_prompt="You are...",       Instructions
    messages=[],                      Conversation history
    hooks=[],                         Hooks
)

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [strands-agents/docs](https://raw.githubusercontent.com/strands-agents/docs/main/designs/0001-plugins.md)*
