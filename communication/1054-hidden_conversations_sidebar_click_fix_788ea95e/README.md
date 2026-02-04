# Hidden Conversations Sidebar Click Fix

| Property | Value |
|----------|-------|
| **Name** | Hidden Conversations Sidebar Click Fix |
| **Repository** | [microsoft/simplechat](https://raw.githubusercontent.com/microsoft/simplechat/main/docs/explanation/fixes/v0.235.001/HIDDEN_CONVERSATIONS_SIDEBAR_CLICK_FIX.md) (⭐ 110) |
| **Original Path** | `docs/explanation/fixes/v0.235.001/HIDDEN_CONVERSATIONS_SIDEBAR_CLICK_FIX.md` |
| **Category** | communication |
| **Subcategory** | messaging |
| **Tags** | communication |
| **Created** | 2026-01-13 |
| **Updated** | 2026-01-13 |
| **File Hash** | `788ea95e68c59b69...` |

## Description

javascript
// If this conversation is hidden, ensure the main conversation list also shows hidden conversations
if (convo.is_hidden && window.chatConversations && window.chatConversations.setShowHiddenConversations) {
  window.chatConversations.setShowHiddenConversations(true);
}

**Tags:** `communication`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [microsoft/simplechat](https://raw.githubusercontent.com/microsoft/simplechat/main/docs/explanation/fixes/v0.235.001/HIDDEN_CONVERSATIONS_SIDEBAR_CLICK_FIX.md)*
