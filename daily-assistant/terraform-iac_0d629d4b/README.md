# Terraform Infrastructure as Code

| Property | Value |
|----------|-------|
| **Name** | Terraform Infrastructure as Code |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/devops-engineer/references/terraform-iac.md) (⭐ 216) |
| **Original Path** | `skills/devops-engineer/references/terraform-iac.md` |
| **Category** | daily-assistant |
| **Subcategory** | tasks |
| **Tags** | daily assistant |
| **Created** | 2025-12-14 |
| **Updated** | 2026-01-29 |
| **File Hash** | `0d629d4b9826cbd9...` |

## Description

hcl
terraform {
  required_providers {
    aws = { source = "hashicorp/aws", version = "~> 5.0" }
  }
  backend "s3" {
    bucket = "terraformstate"
    key    = "app/terraform.tfstate"
    region = "useast1"
  }
}

**Tags:** `daily assistant`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/devops-engineer/references/terraform-iac.md)*
