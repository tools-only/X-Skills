# Pdf Embedding Token Tracking Fix

| Property | Value |
|----------|-------|
| **Name** | Pdf Embedding Token Tracking Fix |
| **Repository** | [microsoft/simplechat](https://raw.githubusercontent.com/microsoft/simplechat/main/docs/explanation/fixes/v0.235.001/PDF_EMBEDDING_TOKEN_TRACKING_FIX.md) (⭐ 110) |
| **Original Path** | `docs/explanation/fixes/v0.235.001/PDF_EMBEDDING_TOKEN_TRACKING_FIX.md` |
| **Category** | content-creation |
| **Subcategory** | media |
| **Tags** | content creation |
| **Created** | 2026-01-13 |
| **Updated** | 2026-01-13 |
| **File Hash** | `ab2e8d7864a1a086...` |

## Description

1. Added token tracking initialization in process_di_document()
python
def process_di_document(...):
      Token tracking initialization 
    total_embedding_tokens = 0
    embedding_model_name = None

**Tags:** `content creation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [microsoft/simplechat](https://raw.githubusercontent.com/microsoft/simplechat/main/docs/explanation/fixes/v0.235.001/PDF_EMBEDDING_TOKEN_TRACKING_FIX.md)*
