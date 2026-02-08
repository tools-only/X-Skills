# Transactions

| Property | Value |
|----------|-------|
| **Name** | Transactions |
| **Repository** | [existential-birds/beagle](https://raw.githubusercontent.com/existential-birds/beagle/main/plugins/beagle-go/skills/go-data-persistence/references/transactions.md) (⭐ 20) |
| **Original Path** | `plugins/beagle-go/skills/go-data-persistence/references/transactions.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-07 |
| **Updated** | 2026-02-07 |
| **File Hash** | `e7b4eec11834128c...` |

## Description

Transactions should be managed at the service layer, not the store (repository) layer. The service layer knows which operations must be atomic. Individual store methods should not start their own transactions because:

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [existential-birds/beagle](https://raw.githubusercontent.com/existential-birds/beagle/main/plugins/beagle-go/skills/go-data-persistence/references/transactions.md)*
