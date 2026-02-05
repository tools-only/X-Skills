# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [existential-birds/beagle](https://raw.githubusercontent.com/existential-birds/beagle/main/skills/pydantic-ai-common-pitfalls/SKILL.md) (⭐ 17) |
| **Original Path** | `skills/pydantic-ai-common-pitfalls/SKILL.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2025-12-21 |
| **Updated** | 2025-12-31 |
| **File Hash** | `503dfc8e23436af7...` |

## Description

python
 ERROR: RunContext not allowed in tool_plain
@agent.tool_plain
async def bad_tool(ctx: RunContext[MyDeps]) > str:
    return "oops"
 UserError: RunContext annotations can only be used with tools that take context

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [existential-birds/beagle](https://raw.githubusercontent.com/existential-birds/beagle/main/skills/pydantic-ai-common-pitfalls/SKILL.md)*
