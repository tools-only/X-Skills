# Signature Verification

| Property | Value |
|----------|-------|
| **Name** | Signature Verification |
| **Repository** | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/vercel-pack/skills/vercel-webhooks-events/references/signature-verification.md) (⭐ 1.3k) |
| **Original Path** | `plugins/saas-packs/vercel-pack/skills/vercel-webhooks-events/references/signature-verification.md` |
| **Category** | investment |
| **Subcategory** | crypto |
| **Tags** | investment |
| **Created** | 2026-01-06 |
| **Updated** | 2026-01-06 |
| **File Hash** | `ea3287657f54ca63...` |

## Description

typescript
function verifyVercelSignature(
  payload: Buffer,
  signature: string,
  timestamp: string
): boolean {
  const secret = process.env.VERCEL_WEBHOOK_SECRET!;

**Tags:** `investment`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/vercel-pack/skills/vercel-webhooks-events/references/signature-verification.md)*
