# Terraform Provider Configuration

| Property | Value |
|----------|-------|
| **Name** | Terraform Provider Configuration |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/terraform-engineer/references/providers.md) (⭐ 216) |
| **Original Path** | `skills/terraform-engineer/references/providers.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2025-12-15 |
| **Updated** | 2026-01-29 |
| **File Hash** | `0eff7b26b2d91ac9...` |

## Description

Basic Configuration
hcl
terraform {
  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = "~> 5.0"
    }
  }
}

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/terraform-engineer/references/providers.md)*
