# Group Action Oauth Schema Merging Fix

| Property | Value |
|----------|-------|
| **Name** | Group Action Oauth Schema Merging Fix |
| **Repository** | [microsoft/simplechat](https://raw.githubusercontent.com/microsoft/simplechat/main/docs/explanation/fixes/v0.237.008/GROUP_ACTION_OAUTH_SCHEMA_MERGING_FIX.md) (⭐ 112) |
| **Original Path** | `docs/explanation/fixes/v0.237.008/GROUP_ACTION_OAUTH_SCHEMA_MERGING_FIX.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-09 |
| **Updated** | 2026-02-09 |
| **File Hash** | `919f798a32fda554...` |

## Description

Critical Discovery: When comparing global vs group action data:
 Global action (working): additionalFields: {auth_method: 'bearer', base_url: '...', ...}
 Group action (failing): additionalFields: {} ← Empty object!

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [microsoft/simplechat](https://raw.githubusercontent.com/microsoft/simplechat/main/docs/explanation/fixes/v0.237.008/GROUP_ACTION_OAUTH_SCHEMA_MERGING_FIX.md)*
