# Query Patterns

| Property | Value |
|----------|-------|
| **Name** | Query Patterns |
| **Repository** | [fcakyon/claude-codex-settings](https://raw.githubusercontent.com/fcakyon/claude-codex-settings/main/plugins/supabase-tools/skills/supabase-usage/references/query-patterns.md) (⭐ 401) |
| **Original Path** | `plugins/supabase-tools/skills/supabase-usage/references/query-patterns.md` |
| **Category** | automation |
| **Subcategory** | scripting |
| **Tags** | automation |
| **Created** | 2025-12-25 |
| **Updated** | 2025-12-25 |
| **File Hash** | `c27804d3c5e1d7a3...` |

## Description

javascript
// Equality
const { data } = await supabase
  .from("users")
  .select("")
  .eq("status", "active");

**Tags:** `automation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [fcakyon/claude-codex-settings](https://raw.githubusercontent.com/fcakyon/claude-codex-settings/main/plugins/supabase-tools/skills/supabase-usage/references/query-patterns.md)*
