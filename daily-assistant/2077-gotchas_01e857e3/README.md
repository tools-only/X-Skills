# Gotchas

| Property | Value |
|----------|-------|
| **Name** | Gotchas |
| **Repository** | [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/agents-sdk/gotchas.md) (🔥 19.8k) |
| **Original Path** | `cli-tool/components/skills/development/cloudflare-deploy/references/agents-sdk/gotchas.md` |
| **Category** | daily-assistant |
| **Subcategory** | notes |
| **Tags** | daily assistant |
| **Created** | 2026-02-08 |
| **Updated** | 2026-02-08 |
| **File Hash** | `01e857e317f073cf...` |

## Description

Cause: Mutating state directly or not calling setState() after modifications  
Solution: Always use setState() with immutable updates:
ts
// ❌ this.state.count++
// ✅ this.setState({...this.state, count: this.state.count + 1})

**Tags:** `daily assistant`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/agents-sdk/gotchas.md)*
