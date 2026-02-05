# Code Injection

| Property | Value |
|----------|-------|
| **Name** | Code Injection |
| **Repository** | [existential-birds/beagle](https://raw.githubusercontent.com/existential-birds/beagle/main/plugins/beagle-elixir/skills/elixir-security-review/references/code-injection.md) (⭐ 17) |
| **Original Path** | `plugins/beagle-elixir/skills/elixir-security-review/references/code-injection.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-05 |
| **Updated** | 2026-02-05 |
| **File Hash** | `0417e51d48af2ae2...` |

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
*Source: [existential-birds/beagle](https://raw.githubusercontent.com/existential-birds/beagle/main/plugins/beagle-elixir/skills/elixir-security-review/references/code-injection.md)*
