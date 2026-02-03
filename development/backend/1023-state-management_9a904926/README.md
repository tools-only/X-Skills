# Terraform State Management

| Property | Value |
|----------|-------|
| **Name** | Terraform State Management |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/terraform-engineer/references/state-management.md) (⭐ 216) |
| **Original Path** | `skills/terraform-engineer/references/state-management.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2025-12-15 |
| **Updated** | 2026-01-29 |
| **File Hash** | `9a904926d847c67f...` |

## Description

Backend Configuration
hcl
 backend.tf
terraform {
  backend "s3" {
    bucket         = "myterraformstate"
    key            = "production/vpc/terraform.tfstate"
    region         = "useast1"
    encrypt        = true
    dynamodb_table = "terraformstatelock"

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/terraform-engineer/references/state-management.md)*
