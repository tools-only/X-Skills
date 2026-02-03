# Terraform Best Practices

| Property | Value |
|----------|-------|
| **Name** | Terraform Best Practices |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/terraform-engineer/references/best-practices.md) (⭐ 216) |
| **Original Path** | `skills/terraform-engineer/references/best-practices.md` |
| **Category** | development |
| **Subcategory** | devops |
| **Tags** | development |
| **Created** | 2025-12-15 |
| **Updated** | 2026-01-29 |
| **File Hash** | `4d62edf91d441d05...` |

## Description

Use Modules for Reusability
hcl
 Bad  Repeated code
resource "aws_vpc" "app1" {
  cidr_block = "10.0.0.0/16"
  enable_dns_hostnames = true
  tags = { Name = "app1vpc", Environment = "prod" }
}

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/terraform-engineer/references/best-practices.md)*
