# Firestore Security Agent

| Property | Value |
|----------|-------|
| **Name** | Firestore Security Agent |
| **Repository** | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/community/jeremy-firestore/agents/firestore-security-agent.md) (⭐ 1.3k) |
| **Original Path** | `plugins/community/jeremy-firestore/agents/firestore-security-agent.md` |
| **Category** | development |
| **Subcategory** | tools |
| **Tags** | development |
| **Created** | 2025-11-11 |
| **Updated** | 2025-12-26 |
| **File Hash** | `fc094b1094c633bf...` |

## Description

javascript
rules_version = '2';
service cloud.firestore {
  match /databases/{database}/documents {
    // Users can only access their own documents
    match /users/{userId} {
      allow read, write: if request.auth != null && request.auth.uid == userId;
    }
  }
}

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/community/jeremy-firestore/agents/firestore-security-agent.md)*
