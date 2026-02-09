# Signature Verification

| Property | Value |
|----------|-------|
| **Name** | Signature Verification |
| **Repository** | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/supabase-pack/skills/supabase-webhooks-events/references/signature-verification.md) (⭐ 1.3k) |
| **Original Path** | `plugins/saas-packs/supabase-pack/skills/supabase-webhooks-events/references/signature-verification.md` |
| **Category** | investment |
| **Subcategory** | crypto |
| **Tags** | investment |
| **Created** | 2026-01-06 |
| **Updated** | 2026-01-06 |
| **File Hash** | `6ca612e1e1d1a775...` |

## Description

typescript
function verifySupabaseSignature(
  payload: Buffer,
  signature: string,
  timestamp: string
): boolean {
  const secret = process.env.SUPABASE_WEBHOOK_SECRET!;

**Tags:** `investment`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/supabase-pack/skills/supabase-webhooks-events/references/signature-verification.md)*
