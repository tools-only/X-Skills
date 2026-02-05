# Code Injection

| Property | Value |
|----------|-------|
| **Name** | Code Injection |
| **Repository** | [existential-birds/beagle](https://raw.githubusercontent.com/existential-birds/beagle/main/skills/elixir-security-review/references/code-injection.md) (⭐ 17) |
| **Original Path** | `skills/elixir-security-review/references/code-injection.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-05 |
| **Updated** | 2026-02-05 |
| **File Hash** | `afeb1df497d0fe21...` |

## Description

elixir
 CRITICAL VULNERABILITY
def calculate(user_expression) do
  {result, _} = Code.eval_string(user_expression)
  result
end

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [existential-birds/beagle](https://raw.githubusercontent.com/existential-birds/beagle/main/skills/elixir-security-review/references/code-injection.md)*
