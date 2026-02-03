# Terraform Module Patterns

| Property | Value |
|----------|-------|
| **Name** | Terraform Module Patterns |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/terraform-engineer/references/module-patterns.md) (⭐ 216) |
| **Original Path** | `skills/terraform-engineer/references/module-patterns.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2025-12-15 |
| **Updated** | 2026-01-29 |
| **File Hash** | `8bbad929701c35c8...` |

## Description

main.tf
hcl
resource "aws_vpc" "this" {
  cidr_block           = var.cidr_block
  enable_dns_hostnames = var.enable_dns_hostnames
  enable_dns_support   = var.enable_dns_support

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/terraform-engineer/references/module-patterns.md)*
