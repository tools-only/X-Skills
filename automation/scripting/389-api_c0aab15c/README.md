# Api

| Property | Value |
|----------|-------|
| **Name** | Api |
| **Repository** | [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/vectorize/api.md) (🔥 19.8k) |
| **Original Path** | `cli-tool/components/skills/development/cloudflare-deploy/references/vectorize/api.md` |
| **Category** | automation |
| **Subcategory** | scripting |
| **Tags** | automation |
| **Created** | 2026-02-08 |
| **Updated** | 2026-02-08 |
| **File Hash** | `c0aab15c2d7fc72b...` |

## Description

typescript
interface VectorizeVector {
  id: string;                    // Max 64 bytes
  values: number[];              // Must match index dimensions
  namespace?: string;            // Optional partition (max 64 bytes)
  metadata?: Record<string, any>; // Max 10 KiB
}

**Tags:** `automation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/vectorize/api.md)*
